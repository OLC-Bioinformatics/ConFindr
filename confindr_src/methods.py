#!/usr/bin/env python3
"""
Core ConFindr methods.

This module implements the low-level helpers used by ConFindr's
contamination-detection pipeline (I/O, BAM parsing, statistics and
reporting helpers).
"""

# Standard library imports
from glob import glob
from itertools import chain
from multiprocessing.dummy import Pool as ThreadPool
from statistics import (
    mean,
    pstdev
)
from typing import (
    Any,
    Dict,
    List,
    Optional,
    Tuple,
)
import argparse
import csv
import gzip
import logging
import math
import multiprocessing
import os
import shlex
import shutil
import subprocess
import sys
import tarfile
import time
import traceback
import urllib.request

# Third-party imports
from Bio import SeqIO
from pysam.utils import SamtoolsError
from scipy.stats import (
    betabinom,
    chi2,
    mannwhitneyu
)
import numpy as np
import psutil
import pysam

# Local imports
from confindr_src.version import __version__
from confindr_src.wrappers import (
    bbtools,
    mash,
)

# Chunking defaults (can be overridden from CLI via confindr.py)
CONTIG_CHUNK_MULTIPLIER = 3  # create threads * multiplier chunks
CONTIG_CHUNK_MAX_BASES = 200_000  # maximum bases per chunk before splitting

# Per-worker FASTQ index state (mutable structure so functions can mutate
# contents without using the 'global' statement). This keeps state local to
# each process (initializer sets indexes in worker processes).
_FASTQ_INDEX_STATE = {'fwd': None, 'rev': None, 'paired': False}

# Minimum base count cutoff for dynamic calculation of expected positions. If
# base_cutoff is set to 0, the expected positions will be calculated based on
# the number of bases at the position, but this minimum cutoff will be applied
# to avoid extremely low cutoffs for low-coverage positions. Setting this to 0
# will disable the safety floor and allow dynamic calculation to produce very
# low cutoffs for low-coverage positions.
MIN_DYNAMIC_CUTOFF = 3


def _format_seconds(
    *,  # Enforce keyword arguments
    s: float
) -> str:
    """
    Format seconds as H:MM:SS

    Args:
        s: Number of seconds (float or int).

    Returns:
        Formatted time string. For example, 3661 seconds becomes '1:01:01'.
        If input cannot be converted to int, returns 'N/A'.
    """
    # Ensure input is numeric and finite; return 'N/A' for invalid inputs
    if not isinstance(s, (int, float)) or not math.isfinite(s):
        return 'N/A'

    try:
        # Round to nearest integer second
        s = int(round(s))
    except (TypeError, ValueError, OverflowError):
        return 'N/A'

    # Calculate hours, minutes, seconds
    h = s // 3600
    m = (s % 3600) // 60
    sec = s % 60

    # Format string based on whether hours are present
    if h:
        return f'{h}:{m:02d}:{sec:02d}'
    return f'{m:02d}:{sec:02d}'


def download_mash_sketch(
    *,  # Enforce keyword arguments
    output_folder: str
) -> None:
    """
    Download the prebuilt RefSeq mash sketch.

    Args:
        output_folder: Directory where the sketch will be saved.

    Returns:
        None
    """
    logging.info('Downloading mash refseq sketch...')
    urllib.request.urlretrieve(
        'https://github.com/OLC-Bioinformatics/ConFindr/raw/master/'
        'refseq_sketch/refseq.msh',
        os.path.join(output_folder, 'refseq.msh')
    )


def download_cgmlst_derived_data(
    *,  # Enforce keyword arguments
    output_folder: str
) -> None:
    """
    Download and extract cgMLST-derived databases.

    The function downloads a tarball containing precomputed databases and
    extracts it into ``output_folder``. The temporary tarball is removed
    after successful extraction.

    Args:
        output_folder: Directory to download and extract files into.

    Returns:
        None
    """
    logging.info(
        'Downloading cgMLST-derived data for Escherichia, Salmonella, '
        'and Listeria...'
    )
    dest = os.path.join(output_folder, 'confindr_db.tar.gz')
    urllib.request.urlretrieve(
        'https://ndownloader.figshare.com/files/14771267',
        dest
    )
    try:
        with tarfile.open(dest) as tar:
            tar.extractall(path=output_folder)
    finally:
        try:
            os.remove(dest)
        except OSError:
            # Best-effort cleanup; ignore failures to remove
            pass
    index(
        output_folder=output_folder,
        genera=['Escherichia', 'Listeria', 'Salmonella'],
        cgderived=True
    )


def index(
    *,  # Enforce keyword arguments
    output_folder: str,
    genera: List[str],
    cgderived: bool = False
) -> None:
    """Ensure genus-specific databases exist and index them.

    Args:
        output_folder: Directory where genus databases live or will be
            created.
        genera: Sequence of genus names to process.
        cgderived: If True, prefer *_db_cgderived.fasta filenames.

    Returns:
        None
    """
    # Iterate through each genus and ensure the database exists and is indexed
    for predominant_genus in genera:
        # Determine the database filename based on cgderived flag
        if cgderived:
            sample_database = os.path.join(
                output_folder,
                f'{predominant_genus}_db_cgderived.fasta'
            )
        else:
            sample_database = os.path.join(
                output_folder,
                f'{predominant_genus}_db.fasta'
            )

        # If the database file does not exist, create it
        if not os.path.isfile(sample_database):

            if os.path.isfile(
                os.path.join(
                    output_folder,
                    'gene_allele.txt'
                )
            ) and os.path.isfile(
                os.path.join(
                    output_folder,
                    'rMLST_combined.fasta'
                )
            ):
                logging.info(
                    'Setting up rMLST genus-specific database for genus %s...',
                    predominant_genus
                )

                # Find the allele list for the predominant genus
                allele_list = find_genus_specific_allele_list(
                    profiles_file=os.path.join(
                        output_folder,
                        'gene_allele.txt'
                    ),
                    target_genus=predominant_genus
                )

                # Create the allele-specific database
                setup_allelespecific_database(
                    fasta_file=sample_database,
                    database_folder=output_folder,
                    allele_list=allele_list
                )

        # Perform the necessary samtools and KMA indexing
        index_databases(sample_database=sample_database)


def run_cmd(
    *,  # Enforce keyword arguments
    cmd: str
) -> Tuple[str, str]:
    """
    Run a shell command and return its stdout and stderr.

    Args:
        cmd: The command to run (as would be typed in the shell).

    Returns:
        A tuple (stdout, stderr) as decoded UTF-8 strings.

    Raises:
        subprocess.CalledProcessError: If the command exits with a non-zero
            return code. `output` and `stderr` attributes will contain the
            captured stdout/stderr.
    """
    p = subprocess.Popen(
        cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    out, err = p.communicate()
    out = out.decode('utf-8')
    err = err.decode('utf-8')
    if p.returncode != 0:
        raise subprocess.CalledProcessError(
            p.returncode,
            cmd,
            output=out,
            stderr=err
        )

    return out, err


def estimate_genome_size(
    *,  # Enforce keyword arguments
    genus: str,
) -> int:
    """
    Return an estimated genome size in bases for a bacterial/archaeal genus.

    These are heuristic genus-level defaults intended for rough parameter
    selection, e.g. expected assembly size, coverage estimation, QC thresholds,
    or downsampling. They are not authoritative strain-level genome sizes.

    If the genus is not recognised, return a broad bacterial default of
    4,000,000 bp.

    Args:
        genus: Genus name, or a taxonomic string beginning with the genus.

    Returns:
        Estimated genome size in bases.
    """
    if not genus or not str(genus).strip():
        return 4_000_000

    g = str(genus).strip().lower().split()[0]

    # Some common historical/renamed genera or spelling variants.
    aliases = {
        "chlamydophila": "chlamydia",
        "propionibacterium": "cutibacterium",
        "ensifer": "sinorhizobium",
        "clostridioides": "clostridioides",
        "lacticaseibacillus": "lacticaseibacillus",
        "lactiplantibacillus": "lactiplantibacillus",
        "lentilactobacillus": "lentilactobacillus",
        "ligilactobacillus": "ligilactobacillus",
        "limosilactobacillus": "limosilactobacillus",
    }
    g = aliases.get(g, g)

    sizes = {
        # ------------------------------------------------------------------
        # Enterobacterales and related common clinical/food genera
        # ------------------------------------------------------------------
        "escherichia": 4_600_000,
        "salmonella": 4_800_000,
        "shigella": 4_600_000,
        "klebsiella": 5_500_000,
        "raoultella": 5_500_000,
        "enterobacter": 4_800_000,
        "cronobacter": 4_500_000,
        "citrobacter": 5_000_000,
        "yersinia": 4_600_000,
        "proteus": 4_000_000,
        "serratia": 5_100_000,
        "morganella": 3_800_000,
        "providencia": 4_300_000,
        "edwardsiella": 3_800_000,
        "hafnia": 4_700_000,
        "kluyvera": 4_800_000,
        "leclercia": 4_800_000,
        "plausibacter": 5_000_000,
        "kosakonia": 5_000_000,
        "pantoea": 4_800_000,
        "erwinia": 4_800_000,
        "pectobacterium": 4_900_000,
        "dickeya": 4_900_000,
        "pragia": 4_500_000,
        "rouxiella": 5_000_000,
        "tatumella": 4_000_000,
        "moellerella": 4_000_000,
        "obesumbacterium": 5_000_000,
        "buttiauxella": 4_700_000,
        "cedecea": 4_800_000,
        "evingella": 5_000_000,
        "grimontia": 5_000_000,
        "sodalis": 4_500_000,
        "xenorhabdus": 4_500_000,
        "photorhabdus": 5_000_000,
        # ------------------------------------------------------------------
        # Pseudomonadota / non-fermenters / environmental opportunists
        # ------------------------------------------------------------------
        "pseudomonas": 6_500_000,
        "acinetobacter": 4_000_000,
        "stenotrophomonas": 4_700_000,
        "burkholderia": 7_500_000,
        "paraburkholderia": 8_000_000,
        "caballeronia": 7_500_000,
        "ralstonia": 5_800_000,
        "cupriavidus": 7_000_000,
        "achromobacter": 6_500_000,
        "alcaligenes": 4_000_000,
        "bordetella": 4_100_000,
        "comamonas": 4_800_000,
        "acidovorax": 5_000_000,
        "delftia": 6_500_000,
        "variovorax": 6_500_000,
        "herbaspirillum": 5_400_000,
        "polaromonas": 5_000_000,
        "janthinobacterium": 6_000_000,
        "collimonas": 5_500_000,
        "massilia": 5_500_000,
        "limnohabitans": 3_500_000,
        "methylobacterium": 6_500_000,
        "methylorubrum": 6_500_000,
        "sphingomonas": 4_200_000,
        "sphingobium": 4_500_000,
        "novosphingobium": 4_200_000,
        "sphingopyxis": 4_200_000,
        "zymomonas": 2_100_000,
        "gluconobacter": 3_300_000,
        "acetobacter": 3_500_000,
        "komagataeibacter": 3_700_000,
        "azospirillum": 7_000_000,
        "magnetospirillum": 5_000_000,
        "rhodospirillum": 4_000_000,
        "caulobacter": 4_000_000,
        "brevundimonas": 3_500_000,
        "phenylobacterium": 4_500_000,
        "maricaulis": 4_000_000,
        "hyphomonas": 3_500_000,
        # ------------------------------------------------------------------
        # Plant-associated Alphaproteobacteria / nitrogen fixers
        # ------------------------------------------------------------------
        "agrobacterium": 5_500_000,
        "rhizobium": 6_700_000,
        "sinorhizobium": 6_700_000,
        "bradyrhizobium": 8_500_000,
        "mesorhizobium": 7_000_000,
        "azorhizobium": 5_300_000,
        "neorhizobium": 6_500_000,
        "allorhizobium": 6_500_000,
        "azotobacter": 6_500_000,
        "beijerinckia": 7_000_000,
        "methylocystis": 4_500_000,
        "methylosinus": 4_500_000,
        # ------------------------------------------------------------------
        # Vibrionales / Aeromonadales / water-associated Gram-negatives
        # ------------------------------------------------------------------
        "vibrio": 4_000_000,
        "photobacterium": 5_000_000,
        "aliivibrio": 4_500_000,
        "aeromonas": 4_700_000,
        "plesiomonas": 3_400_000,
        "tolumonas": 3_500_000,
        "shewanella": 5_000_000,
        "alteromonas": 4_500_000,
        "pseudoalteromonas": 5_000_000,
        "colwellia": 5_000_000,
        "idiomarina": 2_800_000,
        "marinobacter": 4_500_000,
        "thalassomonas": 4_000_000,
        # ------------------------------------------------------------------
        # Xanthomonadales
        # ------------------------------------------------------------------
        "xanthomonas": 5_000_000,
        "xylella": 2_600_000,
        "lysobacter": 4_500_000,
        "dokdonella": 3_800_000,
        "dyella": 4_800_000,
        "luteimonas": 4_500_000,
        # ------------------------------------------------------------------
        # Campylobacterales / related epsilonproteobacteria
        # ------------------------------------------------------------------
        "campylobacter": 1_700_000,
        "helicobacter": 1_700_000,
        "arcobacter": 2_500_000,
        "wolinella": 2_100_000,
        "sulfurospirillum": 3_000_000,
        # ------------------------------------------------------------------
        # Pasteurellaceae and other respiratory/animal-associated genera
        # ------------------------------------------------------------------
        "haemophilus": 1_800_000,
        "pasteurella": 2_300_000,
        "actinobacillus": 2_300_000,
        "mannheimia": 2_600_000,
        "aggregatibacter": 2_200_000,
        "gallibacterium": 2_400_000,
        "histophilus": 2_000_000,
        "avibacterium": 2_400_000,
        # ------------------------------------------------------------------
        # Neisseriaceae and related
        # ------------------------------------------------------------------
        "neisseria": 2_200_000,
        "moraxella": 2_200_000,
        "kingella": 2_000_000,
        "eikenella": 2_200_000,
        "simonsiella": 2_500_000,
        "chromobacterium": 4_700_000,
        # ------------------------------------------------------------------
        # Intracellular / vector-borne / zoonotic Gram-negatives
        # ------------------------------------------------------------------
        "brucella": 3_300_000,
        "bartonella": 1_900_000,
        "rickettsia": 1_300_000,
        "orientia": 2_100_000,
        "ehrlichia": 1_200_000,
        "anaplasma": 1_200_000,
        "wolbachia": 1_300_000,
        "coxiella": 2_000_000,
        "francisella": 1_900_000,
        "legionella": 3_400_000,
        "afipia": 5_000_000,
        # ------------------------------------------------------------------
        # Spirochetes
        # ------------------------------------------------------------------
        "leptospira": 4_600_000,
        "borrelia": 1_500_000,
        "treponema": 1_100_000,
        "brachyspira": 3_200_000,
        "spirochaeta": 3_000_000,
        # ------------------------------------------------------------------
        # Bacillota/Firmicutes: Bacillales and relatives
        # ------------------------------------------------------------------
        "bacillus": 4_200_000,
        "geobacillus": 3_600_000,
        "parageobacillus": 3_600_000,
        "paenibacillus": 6_000_000,
        "brevibacillus": 5_500_000,
        "aneurinibacillus": 5_000_000,
        "virgibacillus": 4_000_000,
        "halobacillus": 4_000_000,
        "oceanobacillus": 4_000_000,
        "lysinibacillus": 4_700_000,
        "solibacillus": 4_500_000,
        "thermobacillus": 3_500_000,
        "staphylococcus": 2_800_000,
        "macrococcus": 2_500_000,
        "mammaliicoccus": 2_500_000,
        "jeotgalicoccus": 2_500_000,
        "salinicoccus": 2_800_000,
        "listeria": 3_000_000,
        "brochothrix": 2_900_000,
        "kurthia": 3_500_000,
        "exiguobacterium": 3_000_000,
        "planococcus": 3_500_000,
        "sporosarcina": 4_000_000,
        # ------------------------------------------------------------------
        # Clostridia / anaerobic Firmicutes
        # ------------------------------------------------------------------
        "clostridium": 4_000_000,
        "clostridioides": 4_300_000,
        "paraclostridium": 4_000_000,
        "perfringens": 3_300_000,
        "desulfotomaculum": 3_500_000,
        "thermoanaerobacter": 2_800_000,
        "thermoanaerobacterium": 3_000_000,
        "caldicellulosiruptor": 2_800_000,
        "acetobacterium": 4_000_000,
        "moorella": 3_000_000,
        "eubacterium": 3_300_000,
        "roseburia": 4_300_000,
        "faecalibacterium": 3_000_000,
        "ruminococcus": 3_500_000,
        "butyrivibrio": 3_500_000,
        "anaerostipes": 3_000_000,
        "coprococcus": 3_000_000,
        "dorea": 3_000_000,
        "blautia": 3_500_000,
        "lachnoclostridium": 3_500_000,
        "oscillibacter": 3_500_000,
        "subdoligranulum": 3_000_000,
        "veillonella": 2_100_000,
        "megasphaera": 2_600_000,
        "megamonas": 2_700_000,
        "selenomonas": 2_500_000,
        "anaerococcus": 2_100_000,
        "peptoniphilus": 2_000_000,
        "peptostreptococcus": 2_000_000,
        "finegoldia": 2_000_000,
        "parvimonas": 1_800_000,
        # ------------------------------------------------------------------
        # Lactic acid bacteria
        # ------------------------------------------------------------------
        "streptococcus": 2_200_000,
        "enterococcus": 3_000_000,
        "lactococcus": 2_500_000,
        "lactobacillus": 2_000_000,
        "lacticaseibacillus": 2_000_000,
        "lactiplantibacillus": 3_000_000,
        "lentilactobacillus": 2_000_000,
        "ligilactobacillus": 2_000_000,
        "limosilactobacillus": 2_000_000,
        "leuconostoc": 2_000_000,
        "pediococcus": 2_000_000,
        "weissella": 2_200_000,
        "carnobacterium": 2_500_000,
        "oenococcus": 1_800_000,
        "tetragenococcus": 2_400_000,
        "aerococcus": 2_000_000,
        "gemella": 1_900_000,
        "vagococcus": 2_200_000,
        # ------------------------------------------------------------------
        # Bacteroidota / gut, oral, environmental
        # ------------------------------------------------------------------
        "bacteroides": 5_200_000,
        "prevotella": 3_500_000,
        "porphyromonas": 2_400_000,
        "parabacteroides": 5_000_000,
        "alistipes": 3_800_000,
        "barnesiella": 3_500_000,
        "odoribacter": 4_000_000,
        "butyricimonas": 4_500_000,
        "paludibacter": 3_500_000,
        "capnocytophaga": 2_800_000,
        "flavobacterium": 3_500_000,
        "chryseobacterium": 4_500_000,
        "elizabethkingia": 4_000_000,
        "weeksella": 2_800_000,
        "pedobacter": 5_000_000,
        "sphingobacterium": 5_000_000,
        "hymenobacter": 5_000_000,
        "maribacter": 4_000_000,
        "zobellia": 5_000_000,
        "tenacibaculum": 4_000_000,
        "formosa": 4_000_000,
        "cytophaga": 4_500_000,
        "runella": 7_000_000,
        # ------------------------------------------------------------------
        # Actinobacteria / Actinomycetota
        # ------------------------------------------------------------------
        "mycobacterium": 4_400_000,
        "mycolicibacterium": 6_000_000,
        "mycolicibacter": 5_500_000,
        "corynebacterium": 2_800_000,
        "cutibacterium": 2_500_000,
        "bifidobacterium": 2_200_000,
        "actinomyces": 3_000_000,
        "nocardia": 6_500_000,
        "rhodococcus": 5_500_000,
        "gordonia": 5_000_000,
        "dietzia": 3_500_000,
        "tsukamurella": 4_500_000,
        "streptomyces": 8_500_000,
        "kitasatospora": 8_500_000,
        "micromonospora": 7_000_000,
        "saccharopolyspora": 7_000_000,
        "amycolatopsis": 9_000_000,
        "nocardiopsis": 5_500_000,
        "frankia": 7_500_000,
        "micrococcus": 2_600_000,
        "arthrobacter": 4_500_000,
        "paenarthrobacter": 4_500_000,
        "pseudarthrobacter": 4_500_000,
        "kocuria": 2_800_000,
        "microbacterium": 3_500_000,
        "leifsonia": 3_500_000,
        "curtobacterium": 3_700_000,
        "plantibacter": 3_500_000,
        "agromyces": 3_500_000,
        "brevibacterium": 4_000_000,
        "cellulomonas": 4_000_000,
        "dermacoccus": 3_000_000,
        "janibacter": 3_500_000,
        "jonesia": 3_000_000,
        "mobiluncus": 2_200_000,
        "gardnerella": 1_700_000,
        "collinsella": 2_000_000,
        "egerthella": 2_000_000,
        "slackia": 2_000_000,
        "olsenella": 2_000_000,
        "atopobium": 1_700_000,
        # ------------------------------------------------------------------
        # Mollicutes / reduced-genome bacteria
        # ------------------------------------------------------------------
        "mycoplasma": 800_000,
        "ureaplasma": 750_000,
        "acholeplasma": 1_500_000,
        "mesoplasma": 900_000,
        "spiroplasma": 1_300_000,
        "phytoplasma": 700_000,
        # ------------------------------------------------------------------
        # Chlamydiae and related intracellular bacteria
        # ------------------------------------------------------------------
        "chlamydia": 1_000_000,
        "parachlamydia": 2_400_000,
        "simkania": 2_500_000,
        "waddlia": 2_100_000,
        # ------------------------------------------------------------------
        # Cyanobacteria
        # ------------------------------------------------------------------
        "synechococcus": 2_700_000,
        "prochlorococcus": 1_800_000,
        "nostoc": 7_500_000,
        "anabaena": 6_500_000,
        "microcystis": 5_000_000,
        "cyanothece": 5_000_000,
        "fischerella": 7_000_000,
        "calothrix": 7_000_000,
        "gleobacter": 4_600_000,
        "arthrospira": 6_500_000,
        "spirulina": 6_000_000,
        "oscillatoria": 6_500_000,
        "planktothrix": 5_000_000,
        "crocosphaera": 6_000_000,
        "trichodesmium": 7_000_000,
        # ------------------------------------------------------------------
        # Deinococcus-Thermus and thermophiles
        # ------------------------------------------------------------------
        "deinococcus": 3_300_000,
        "thermus": 2_200_000,
        "meiothermus": 3_000_000,
        "thermotoga": 1_900_000,
        "petrotoga": 2_000_000,
        "furcifer": 2_000_000,
        "aquifex": 1_600_000,
        "hydrogenobacter": 1_800_000,
        # ------------------------------------------------------------------
        # Planctomycetes / Verrucomicrobia / PVC superphylum
        # ------------------------------------------------------------------
        "planctomyces": 6_000_000,
        "gemmata": 8_000_000,
        "rhodopirellula": 7_000_000,
        "blastopirellula": 7_000_000,
        "pirellula": 7_000_000,
        "akkermansia": 2_700_000,
        "verrucomicrobium": 6_000_000,
        "prosthecobacter": 6_000_000,
        "chthoniobacter": 5_000_000,
        # ------------------------------------------------------------------
        # Acidobacteria and soil-associated groups
        # ------------------------------------------------------------------
        "acidobacterium": 5_000_000,
        "granulicella": 5_000_000,
        "terracidiphilus": 5_000_000,
        "bryobacter": 7_000_000,
        "edaphobacter": 7_000_000,
        "koribacter": 6_000_000,
        "solibacter": 9_000_000,
        # ------------------------------------------------------------------
        # Chloroflexi and related environmental bacteria
        # ------------------------------------------------------------------
        "chloroflexus": 5_000_000,
        "roseiflexus": 5_500_000,
        "herpetosiphon": 6_500_000,
        "anaerolinea": 4_500_000,
        "dehalococcoides": 1_500_000,
        "dehalogenimonas": 1_600_000,
        # ------------------------------------------------------------------
        # Nitrospirae / nitrifiers / sulfur and iron bacteria
        # ------------------------------------------------------------------
        "nitrospira": 4_500_000,
        "nitrobacter": 4_000_000,
        "nitrosomonas": 3_000_000,
        "nitrosospira": 3_500_000,
        "nitrosococcus": 3_500_000,
        "thiobacillus": 3_500_000,
        "acidithiobacillus": 3_000_000,
        "beggiatoa": 5_000_000,
        "thiomicrospira": 2_500_000,
        "allochromatium": 3_500_000,
        "chromatium": 4_000_000,
        # ------------------------------------------------------------------
        # Desulfobacterota / sulfate reducers
        # ------------------------------------------------------------------
        "desulfovibrio": 3_800_000,
        "desulfobacter": 4_000_000,
        "desulfobulbus": 4_000_000,
        "desulfococcus": 4_000_000,
        "desulfotalea": 3_500_000,
        "desulfomicrobium": 3_000_000,
        "desulfosporosinus": 5_000_000,
        "desulfitobacterium": 5_000_000,
        "sulfurovum": 2_500_000,
        # ------------------------------------------------------------------
        # Fusobacteria
        # ------------------------------------------------------------------
        "fusobacterium": 2_400_000,
        "leptotrichia": 2_300_000,
        "streptobacillus": 1_700_000,
        # ------------------------------------------------------------------
        # Archaea, optional
        # ------------------------------------------------------------------
        "methanobrevibacter": 2_000_000,
        "methanococcus": 1_700_000,
        "methanocaldococcus": 1_800_000,
        "methanosarcina": 4_500_000,
        "methanosaeta": 3_000_000,
        "methanothrix": 3_000_000,
        "methanobacterium": 2_800_000,
        "methanoculleus": 2_500_000,
        "halobacterium": 2_600_000,
        "haloferax": 3_900_000,
        "halorubrum": 3_000_000,
        "halococcus": 3_000_000,
        "natronomonas": 3_000_000,
        "sulfolobus": 3_000_000,
        "saccharolobus": 3_000_000,
        "thermococcus": 2_000_000,
        "pyrococcus": 1_900_000,
        "archaeoglobus": 2_200_000,
        "thermoplasma": 1_600_000,
        "picrophilus": 1_600_000,
    }

    return sizes.get(g, 4_000_000)


def estimate_mean_read_length(
    *,  # Enforce keyword arguments
    fastq_path: str,
    sample_reads: int = 1000
) -> int:
    """
    Estimate mean read length by sampling up to `sample_reads` reads from
    the file. Uses gzip if file ends with .gz for streaming.

    Args:
        fastq_path: Path to FASTQ file (can be gzipped).
        sample_reads: Number of reads to sample for estimation.

    Returns:
        Estimated mean read length (integer).
    """
    # Initialize counters
    total = 0
    count = 0

    # Open the FASTQ file, using gzip if necessary
    opener = gzip.open if fastq_path.endswith('.gz') else open
    mode = 'rt' if fastq_path.endswith('.gz') else 'r'

    # Attempt to read the FASTQ file and sample reads
    try:
        with opener(
            fastq_path,
            mode=mode,
            encoding='utf-8',
            errors='ignore'
        ) as fh:
            for _ in fh:
                # only count header + sequence + plus + qual lines as grouped
                # fastq format: sequence lines at every 2nd line if we iterate
                # naive approach: read in blocks of 4 lines
                seq = fh.readline()
                if not seq:
                    break
                # Strip whitespace and count length
                seq = seq.strip()
                total += len(seq)
                count += 1

                # If we've sampled enough reads, stop
                if count >= sample_reads:
                    break
    except (OSError, IOError):
        # If we can't read the file, return a default
        return 150

    return int(max(50, total // max(1, count)))


def count_fastq_reads(
    *,  # Enforce keyword arguments
    fastq_path: str
) -> int:
    """
    Fast count of number of reads in a FASTQ file using wc -l; returns
    number of reads (not lines). Handles gzipped files via zcat.

    Args:
        fastq_path: Path to FASTQ file (can be gzipped).

    Returns:
        Number of reads in the FASTQ file (integer).
    """
    # Construct command based on whether file is gzipped
    if fastq_path.endswith('.gz'):
        cmd = f"zcat {shlex.quote(fastq_path)} | wc -l"
    else:
        cmd = f"wc -l {shlex.quote(fastq_path)} | awk '{'{print $1}'}'"

    # Run the command
    out, _ = run_cmd(cmd=cmd)

    # Parse the output to get number of lines
    try:
        lines = int(out.strip())
    except (ValueError, TypeError):
        # If we can't parse the output, return 0
        return 0

    # Each read consists of 4 lines in FASTQ format
    return max(0, lines // 4)


def downsample_reads(
    *,  # Enforce keyword arguments
    pair: List[str],
    sample_tmp_dir: str,
    sample_name: str,
    target_reads: int,
    log: str,
    xmx: Optional[str] = None,
    seed: Optional[int] = None
) -> List[str]:
    """
    Downsample reads using BBTools reformat.sh (samplereadstarget). Returns
    new pair/list of file paths (same shape as input). Raises
    subprocess.CalledProcessError on failure.

    If `seed` is provided, it will be passed to reformat.sh as `sampleseed`
    to enable deterministic, repeatable downsampling across replicates.

    Args:
        pair: List of one (unpaired) or two (paired) FASTQ file paths.
        sample_tmp_dir: Path to temporary directory to write downsampled
            files to.
        sample_name: Sample name to use for output files.
        target_reads: Target number of reads after downsampling.
        log: Path to logfile to write reformat.sh output to.
        xmx: Optional Java Xmx memory setting (e.g. '4g').
        seed: Optional integer seed for random downsampling.

    Returns:
        A list of one or two file paths to the downsampled FASTQ files.
    """
    # Prepare seed argument if provided
    seed_arg = f" sampleseed={seed}" if seed is not None else ""

    # Paired reads
    if len(pair) == 2:
        # Unpack the pair
        r1, r2 = pair

        # Define output file paths
        out1 = os.path.join(
            sample_tmp_dir,
            f'{sample_name}_downsampled_R1.fastq.gz'
        )
        out2 = os.path.join(
            sample_tmp_dir,
            f'{sample_name}_downsampled_R2.fastq.gz'
        )

        # Construct the reformat.sh command
        cmd = (
            f"reformat.sh in={shlex.quote(r1)} in2={shlex.quote(r2)}"
            f" out={shlex.quote(out1)} out2={shlex.quote(out2)} "
            f"samplereadstarget={target_reads} overwrite=true{seed_arg} "
            f"-Xmx{xmx}"
        )

        logging.debug('Downsampling paired reads with command: %s', cmd)

        # Run the command
        out, err = run_cmd(cmd=cmd)

        # Write output to logfile
        write_to_logfile(logfile=log, out=out, err=err, cmd=cmd)

        return [out1, out2]

    # Unpaired reads
    r = pair[0]

    # Define output file path
    out = os.path.join(
        sample_tmp_dir,
        f'{sample_name}_downsampled.fastq.gz'
    )

    # Construct the reformat.sh command
    cmd = (
        f"reformat.sh in={shlex.quote(r)} out={shlex.quote(out)} "
        f"samplereadstarget={target_reads} overwrite=true{seed_arg} "
        f"-Xmx{xmx}"
    )

    logging.debug('Downsampling unpaired reads with command: %s', cmd)

    # Run the command
    out, err = run_cmd(cmd=cmd)

    # Write output to logfile
    write_to_logfile(logfile=log, out=out, err=err, cmd=cmd)

    return [out]


def write_to_logfile(
    *,  # Enforce keyword arguments
    logfile: str,
    out: str,
    err: str,
    cmd: str
) -> None:
    """
    Append command, stdout and stderr to a logfile.

    Args:
        logfile: Path to file to write output to.
        out: Stdout of program called, as a string.
        err: Stderr of program called, as a string.
        cmd: Command that was used.
    """
    with open(logfile, 'a+', encoding='utf-8') as outfile:
        outfile.write(f'Command used: {cmd}\n\n')
        outfile.write(f'STDOUT: {out}\n\n')
        outfile.write(f'STDERR: {err}\n\n')


def dependency_check(
    *,  # Enforce keyword arguments
    dependency: str
) -> bool:
    """
    Check whether a command-line dependency is available on PATH.

    Args:
        dependency: The executable name to check for (e.g. 'blastn').

    Returns:
        True if the executable is found on PATH, False otherwise.
    """
    return shutil.which(dependency) is not None


def find_paired_reads(
    *,  # Enforce keyword arguments
    fastq_directory: str,
    forward_id: str = '_R1',
    reverse_id: str = '_R2'
) -> List[List[str]]:
    """
    Find paired FASTQ files in a directory.

    Args:
        fastq_directory: Path containing FASTQ files.
        forward_id: Identifier substring for forward reads (default: '_R1').
        reverse_id: Identifier substring for reverse reads (default: '_R2').

    Returns:
        A list of [forward, reverse] pairs (paths as strings).
    """
    # Define list to hold pairs
    pair_list: List[List[str]] = []

    # Find all FASTQ files in the directory
    fastq_files = glob(os.path.join(fastq_directory, '*.f*q*'))
    for name in sorted(fastq_files):
        # Check to see if the file is a forward read and that the corresponding
        # reverse read exists
        if forward_id in name and os.path.isfile(
            name.replace(
                forward_id,
                reverse_id
            )
        ):
            pair_list.append([name, name.replace(forward_id, reverse_id)])

    return pair_list


def find_unpaired_reads(
    *,  # Enforce keyword arguments
    fastq_directory: str,
    forward_id: str = '_R1',
    reverse_id: str = '_R2',
    find_fasta: bool = False
) -> List[List[str]]:
    """
    Find unpaired FASTQ (or FASTA) files in a directory.

    Args:
        fastq_directory: Path to directory containing reads.
        forward_id: Identifier substring for forward reads (default: '_R1').
        reverse_id: Identifier substring for reverse reads (default: '_R2').
        find_fasta: If True, search for FASTA (*.f*a*). Otherwise search for
        FASTQ.

    Returns:
        A list of single-element lists containing paths to unpaired read files.
    """
    # Define list to hold unpaired reads
    read_list = []

    # Find all FASTQ or FASTA files in the directory
    if find_fasta is False:
        sequence_files = glob(os.path.join(fastq_directory, '*.f*q*'))
    else:
        # Find FASTA files
        sequence_files = glob(os.path.join(fastq_directory, '*.f*a*'))

    # Iterate through files, adding them to our list of unpaired reads if:
    for name in sorted(sequence_files):
        # 1) They don't have the forward identifier or the reverse identifier
        # in their name.
        if forward_id not in name and reverse_id not in name:
            read_list.append([name])

        # 2) They have forward but the reverse isn't there.
        elif forward_id in name and not os.path.isfile(
            name.replace(
                forward_id,
                reverse_id
            )
        ):
            read_list.append([name])

        # 3) They have reverse but the forward isn't there.
        elif reverse_id in name and not os.path.isfile(
            name.replace(
                reverse_id,
                forward_id
            )
        ):
            read_list.append([name])

    return read_list


def find_genus_specific_allele_list(
    *,  # Enforce keyword arguments
    profiles_file: str,
    target_genus: str
) -> List[str]:
    """
    Return allele list for a target genus from a profiles file.

    Args:
        profiles_file: Path to profiles file containing genus:allele1,
        allele2,... lines.
        target_genus: Genus for which alleles should be returned.

    Returns:
        A list of allele identifiers for the target genus.
    """
    # Initialize empty allele list
    alleles = []

    # Read through the profiles file to find alleles for the target genus
    with open(profiles_file, encoding='utf-8') as f:
        lines = f.readlines()

    # Parse the lines to find the target genus
    for line in lines:
        line = line.rstrip()
        genus = line.split(':')[0]

        # If this line corresponds to the target genus, extract alleles
        if genus == target_genus:
            alleles = line.split(':')[1].split(',')[:-1]

    return alleles


def setup_allelespecific_database(
    *,  # Enforce keyword arguments
    fasta_file: str,
    database_folder: str,
    allele_list: List[str]
) -> None:
    """
    Create a genus-specific FASTA file containing only the alleles listed.

    Args:
        fasta_file: Path to FASTA file to write allele-specific database to.
        database_folder: Path where rMLST_combined.fasta can be found.
        allele_list: List of allele identifiers to include.
    """
    # Create an index of the rMLST combined FASTA file
    rmlst_index = SeqIO.index(
        os.path.join(
            database_folder,
            'rMLST_combined.fasta'
        ),
        'fasta'
    )

    # Create a list to hold the sequences to write
    seqs = []

    # Iterate through the allele list, adding sequences to the seqs list
    for s in allele_list:
        try:
            seqs.append(rmlst_index[s])
        except KeyError:
            logging.warning(
                'Tried to add %s to allele-specific database, but could not '
                'find it.', s
            )
    try:
        SeqIO.write(seqs, fasta_file, 'fasta')
    except FileNotFoundError:
        pass


def find_cross_contamination(
    *,  # Enforce keyword arguments
    databases: str,
    reads: Any,
    sample_name: str,
    tmpdir: str = 'tmp',
    log: str = 'log.txt',
    threads: int = 1,
    min_matching_hashes: int = 40
) -> str:
    """
    Uses mash to find out whether or not a sample has more than one genus
    present, indicating cross-contamination.

    Args:
        databases: Path to folder containing mash sketch database (refseq.msh).
        reads: Either a string path to single-end reads, or a list of two
        string paths for paired-end reads.
        sample_name: Sample name to use for temporary files.
        tmpdir: Path to temporary directory to use.
        log: Path to logfile to write mash output to.
        threads: Number of threads to use.
        min_matching_hashes: Minimum number of matching hashes to consider a
        genus present.

    Returns:
        A string representing the genera present. If only one genus is found,
        the string is that genus. If no genera are found, the string is 'ND'.
        If more than one genus is found, the string is a list of genera
        present, separated by colons.
    """
    # Initialize empty list to hold genera present
    genera_present = []

    # Only run the MASH analyses if the screen.tab output file does not
    # already exist
    screen_file = os.path.join(
        tmpdir,
        f'{sample_name}_screen.tab'
    )

    # Run mash screen if output file doesn't already exist
    if not os.path.isfile(screen_file):
        # If reads is a string, it's unpaired. If it's a list, it's paired
        if isinstance(reads, str):
            out, err, cmd = mash.screen(
                f'{databases}/refseq.msh',
                reads,
                threads=threads,
                w='',
                i='0.85',
                output_file=screen_file,
                returncmd=True
            )
        else:
            out, err, cmd = mash.screen(
                f'{databases}/refseq.msh',
                reads[0],
                reads[1],
                threads=threads,
                w='',
                i='0.85',
                output_file=screen_file,
                returncmd=True
            )

        # Write mash output to logfile
        write_to_logfile(
            logfile=log,
            out=out,
            err=err,
            cmd=cmd
        )

    # Read mash screen output
    screen_output = mash.read_mash_screen(
        os.path.join(
            tmpdir,
            f'{sample_name}_screen.tab'
        )
    )

    # Parse mash screen output to find genera present
    for item in screen_output:
        # Extract genus from the query ID
        mash_genus = item.query_id.split('/')[-3]

        # Convert Shigella to Escherichia
        if 'Shigella' in mash_genus:
            mash_genus = 'Escherichia'

        # Extract number of matching hashes
        matching_hashes = int(item.shared_hashes.split('/')[0])

        # Only add the genus to the genera_present list of the number of
        # matching hashes exceeds the cutoff
        if matching_hashes >= min_matching_hashes:
            if mash_genus not in genera_present:
                genera_present.append(mash_genus)

    # Format the output string depending on how many genera were found
    if len(genera_present) == 1:
        genera_present = genera_present[0]
    elif len(genera_present) == 0:
        genera_present = 'ND'
    else:
        # Concatenate multiple genera with colons
        tmpstr = ''
        for mash_genus in genera_present:
            tmpstr += mash_genus + ':'
        genera_present = tmpstr[:-1]

    return genera_present


def number_of_bases_above_threshold(
    *,  # Enforce keyword arguments
    high_quality_base_count: Dict[str, int],
    base_count_cutoff: int = 2,
    base_fraction_cutoff: Optional[float] = None
) -> int:
    """
    Return number of bases that meet configured high-quality thresholds.

    Args:
        high_quality_base_count: Mapping of base->count at a position.
        base_count_cutoff: Absolute count cutoff for a base to be considered
        high-quality.
        base_fraction_cutoff: Optional fraction of total coverage required for
        base to count.

    Returns:
        The number of bases that meet the threshold (integer).
    """
    # Dictionary comprehension where values are True or False for each base
    # depending on whether the count meets the threshold.
    # Method differs depending on whether absolute or fraction cutoff is
    # specified
    if base_fraction_cutoff:
        # Calculate total high-quality base count for fraction comparison
        total_hq_base_count = sum(high_quality_base_count.values())

        # Calculate which bases meet both the fraction and count cutoffs
        bases_above_threshold = {
            base:
                float(count) / total_hq_base_count >= base_fraction_cutoff
                and count >= base_count_cutoff for (base, count)
                in high_quality_base_count.items()
        }
    else:
        bases_above_threshold = {
            base:
                count >= base_count_cutoff for (base, count)
                in high_quality_base_count.items()
            }

    # True is equal to 1 so sum of the number of Trues in the
    # bases_above_threshold dict is the number of bases passing threshold
    return sum(bases_above_threshold.values())


def parse_bam(
    *,  # Enforce keyword arguments
    bamfile_name: str,
    contig_name: str,
    pysam_fasta: Any,
    start: Optional[int] = None,
    end: Optional[int] = None
) -> Tuple[Any, Any]:
    """
    Open a BAM file with pysam and produce a pileup iterator for a contig or
    contig region.

    Args:
        bamfile_name: Path to sorted, indexed BAM file.
        contig_name: Contig/sequence name to produce pileup for.
        pysam_fasta: A pysam.FastaFile (or compatible) used for reference base.
        start: Optional 0-based start coordinate (inclusive).
        end: Optional 0-based end coordinate (exclusive).

    Returns:
        A tuple (bamfile, pileup) where bamfile is the opened AlignmentFile and
        pileup is the iterator produced by pysam's pileup() for the contig
        or contig region.
    """
    # Load the sorted BAM-formatted file using pysam
    bamfile = pysam.AlignmentFile(bamfile_name, 'rb')

    # Normalize contig_name to str (some pysam variants or Python versions may
    # expose contig names as bytes which leads to confusing errors like
    # "reference sequence for 'b'...'' not found").
    if isinstance(contig_name, bytes):
        try:
            contig_name = contig_name.decode('utf-8')
        except UnicodeDecodeError:
            # Fall back to a lossy decode rather than catching all Exceptions
            contig_name = contig_name.decode('utf-8', errors='replace')

    # These parameters seem to be fairly undocumented with pysam, but I think
    # that they should make the output that I'm getting to match up with
    # what I'm seeing in Tablet.
    if start is None and end is None:
        pileup = bamfile.pileup(
            contig_name,
            stepper='samtools',
            ignore_orphans=False,
            fastafile=pysam_fasta,
            min_base_quality=0
        )
    else:
        # Provide explicit start/end to pysam.pileup for region-limited pileups
        pileup = bamfile.pileup(
            contig_name,
            start if start is not None else 0,
            end if end is not None else None,
            stepper='samtools',
            ignore_orphans=False,
            fastafile=pysam_fasta,
            min_base_quality=0
        )

    return bamfile, pileup


def characterise_read(
    *,  # Enforce keyword arguments
    column: Any,
    reference_sequence: str,
    fastq_records: Dict[str, Any],
    quality_cutoff: int,
    min_quality: int = 15,
    fasta: bool = False,
    nanopore: bool = False
) -> Dict[str, Any]:
    """
    Extract read-level characteristics from a pileup column.

    Args:
        column: Pysam pileup column object for the position.
        reference_sequence: The full reference sequence string for the contig.
        fastq_records: Mapping of read name -> SeqRecord with quality info.
        quality_cutoff: Minimum base quality used for read trimming and
            initial filtering.
        min_quality: Minimum base quality required to count a base as SNV
            support.

    Keyword Args:
        fasta: If True, operate in FASTA-only mode (no quality checks).
        nanopore: If True, use nanopore-specific heuristics.

    Returns:
        A dict containing read characterisations (bases, qualities, mapping
        qualities, strand information, etc.)."
    """
    # Initialise a dictionary to store the base details
    filtered_read_dict = {
        'congruent_SNV': {},
        'congruent_ref': {},
        'forward_SNV_reverse_SNV1': {},
        'reverse_SNV_forward_SNV1': {},
        'forward_SNV_reverse_ref': {},
        'reverse_SNV_forward_ref': {},
        'forward_SNV_reverse_UM_QF': {},
        'forward_ref_reverse_UM_QF': {},
        'forward_quality_filtered': {},
        'reverse_SNV_forward_UM_QF': {},
        'reverse_ref_forward_UM_QF': {},
        'reverse_quality_filtered': {}
    }

    # Initialise a dictionary to store the details parsed from the pileup
    unfiltered_read_details = {}

    # Initialise a list to store all the phred scores for the bases passing
    # filter in the column
    qualities = []

    # Initialise a dict to store per-base support metadata for statistical
    # tests
    base_support = {}

    # Iterate through every read present in the column of the pileup
    for read in column.pileups:
        # Not entirely sure why this is sometimes None, but it causes bad stuff
        if read.query_position is not None:
            #  Initialise the read name in the dictionary as required
            if read.alignment.qname not in unfiltered_read_details:
                unfiltered_read_details[read.alignment.qname] = {}

            # Extract the sequence of the base in the read
            query_base = read.alignment.query_sequence[read.query_position]

            # Extract the sequence of the base in the reference gene
            ref_base = reference_sequence[column.pos]

            # Create a boolean of whether the query base matches the reference
            # base (is it a SNV?)
            match = query_base == ref_base

            # Read names in the pileup have the direction removed - add it
            # back for future parsing. Not an issue for FASTA files
            if not fasta and not nanopore:
                if read.alignment.is_read1:
                    if read.alignment.qname.split(' ')[0].endswith('/1'):
                        read_name = read.alignment.qname.split(' ')[0]
                    else:
                        read_name = read.alignment.qname.split(' ')[0] + '/1'
                else:
                    if read.alignment.qname.split(' ')[0].endswith('/2'):
                        read_name = read.alignment.qname.split(' ')[0]
                    else:
                        read_name = read.alignment.qname.split(' ')[0] + '/2'
            else:
                read_name = read.alignment.qname

            # Extract the phred quality score from the FASTQ records,
            # using on-demand lookup to avoid keeping the whole FASTQ in
            # memory when running in parallel workers.
            rec = None
            if fastq_records:
                rec = fastq_records.get(read_name)
            else:
                rec = _get_fastq_record(read_name)
            if rec is None:
                # Could not find the read in FASTQ; fallback to quality 0
                quality = 0
            else:
                # Explicit checks for letter_annotations attribute
                quals = None

                # Check for letter_annotations dict and extract phred_quality
                if hasattr(
                    rec,
                    'letter_annotations'
                ) and isinstance(
                    rec.letter_annotations, dict
                ):
                    quals = rec.letter_annotations.get("phred_quality")
                if isinstance(
                    quals,
                    (list, tuple)
                ) and 0 <= read.query_position < len(
                    quals
                ):
                    quality = quals[read.query_position]
                else:
                    quality = 0

            # Initialise a dictionary to store the range
            range_dict = {}

            # Determine whether there are SNVs clustered together - they will
            # be discarded from the analysis. Iterate through a range of the
            # five positions preceding, and five subsequent positions
            for iterator, contig_pos in enumerate(
                chain(
                    range(
                        column.pos - 5, column.pos
                    ),
                    range(
                        column.pos + 1, column.pos + 6
                    )
                )
            ):
                # Ensure that the contig position being examined isn't beyond
                # the length of the gene
                if 0 <= contig_pos < len(reference_sequence) - 1:
                    # Calculate the read position corresponding to the current
                    # column position
                    read_pos = list(
                        chain(
                            range(
                                read.query_position - 5, read.query_position
                            ),
                            range(
                                read.query_position + 1,
                                read.query_position + 6
                            )
                        )
                    )[iterator]
                    # Ensure that read position being examined isn't beyond
                    # the length of the read
                    if 0 <= read_pos < read.alignment.query_alignment_end - 1:
                        # Populate the dictionary with the calculated positions
                        range_dict[contig_pos] = read_pos

            # Initialise a boolean of whether the current base passes filters
            # and should be added to the dictionary
            add_base = True

            # Iterate through the range_dict to extract the sequence of the
            # bases in the range for both the read and the reference sequence
            for contig_pos, read_pos in range_dict.items():
                try:
                    # Extract the sequence of the base
                    reference_base = reference_sequence[contig_pos]
                    base = read.alignment.query_sequence[read_pos]
                    # If any of the downstream or upstream bases do not match,
                    # set the boolean to False
                    if not match and reference_base != base:
                        add_base = False
                except KeyError:
                    pass

            # Populate the dictionary only if there are no other SNVs within
            # five downstream and five upstream bases
            if add_base:
                unfiltered_read_details[
                    read.alignment.qname
                ][read.alignment.is_read1] = {
                    'mate unmapped': read.alignment.mate_is_unmapped,
                    'paired': read.alignment.is_paired,
                    'forward': read.alignment.is_read1,
                    'reverse': read.alignment.is_read2,
                    'match': match,
                    'qbase': query_base,
                    'rbase': ref_base,
                    'qual': quality,
                    'pos': column.pos,
                    'gene': column.reference_name
                }

                # Add the quality of the base only if it meets the SNV
                # support threshold.
                if quality >= min_quality:
                    qualities.append(quality)
                    # record per-base support qualities, mapping qualities and
                    # strand counts for later statistical tests
                    if query_base not in base_support:
                        base_support[query_base] = {
                            'quals': [],
                            'forward_quals': [],
                            'reverse_quals': [],
                            'forward': 0,
                            'reverse': 0,
                            'mapqs': [],
                            'forward_mapqs': [],
                            'reverse_mapqs': [],
                            'positions': [],
                            'forward_positions': [],
                            'reverse_positions': []
                        }
                    base_support[query_base]['quals'].append(quality)

                    # Mapping quality and read position
                    try:
                        mapq = int(read.alignment.mapping_quality)
                    except TypeError:
                        mapq = None

                    # Add mapping quality if it exists
                    if mapq is not None:
                        base_support[query_base]['mapqs'].append(mapq)

                    # Add strand and position information
                    base_support[query_base]['positions'].append(
                        read.query_position
                    )

                    # Strand-specific information
                    if read.alignment.is_read1:
                        base_support[query_base]['forward_quals'].append(
                            quality
                        )

                        # Increment forward strand count
                        base_support[query_base]['forward'] += 1

                        # Add mapping quality if it exists
                        if mapq is not None:
                            base_support[query_base]['forward_mapqs'].append(
                                mapq
                            )

                        # Add read position
                        base_support[query_base]['forward_positions'].append(
                            read.query_position
                        )
                    else:
                        # Reverse strand
                        base_support[query_base]['reverse_quals'].append(
                            quality
                        )

                        # Increment reverse strand count
                        base_support[query_base]['reverse'] += 1

                        # Add mapping quality if it exists
                        if mapq is not None:
                            base_support[query_base]['reverse_mapqs'].append(
                                mapq
                            )

                        # Add read position
                        base_support[query_base]['reverse_positions'].append(
                            read.query_position
                        )

    # Parse the unfiltered reads to characterise the bases
    for _, dir_dict in unfiltered_read_details.items():
        # Check to see if paired reads are present at this position
        if len(dir_dict) > 1:
            # SNV in both forward and reverse reads
            if not dir_dict[True]['match'] and not dir_dict[False]['match']:
                # Same SNV sequence - don't quality filter as the reads
                # agreeing acts as a quality check
                if dir_dict[True]['qbase'] == dir_dict[False]['qbase']:
                    # Forward and reverse SNV - add two matches (for forward
                    # and reverse reads)
                    if (
                        dir_dict[True]['qbase']
                        not in filtered_read_dict['congruent_SNV']
                    ):
                        filtered_read_dict[
                            'congruent_SNV'
                        ][dir_dict[True]['qbase']] = 2
                    else:
                        filtered_read_dict[
                            'congruent_SNV'
                        ][dir_dict[True]['qbase']] += 2

                # Different SNV sequences
                else:
                    # Both SNV sequences pass quality
                    if (
                        dir_dict[True]['qual'] >= min_quality
                        and dir_dict[False]['qual'] >= min_quality
                    ):
                        # Forward SNV1
                        if (
                            dir_dict[True]['qbase']not in
                            filtered_read_dict['forward_SNV_reverse_SNV1']
                        ):
                            filtered_read_dict[
                                'forward_SNV_reverse_SNV1'
                            ][dir_dict[True]['qbase']] = 1
                        else:
                            filtered_read_dict[
                                'forward_SNV_reverse_SNV1'
                            ][dir_dict[True]['qbase']] += 1

                        # Reverse SNV2
                        if (
                            dir_dict[False]['qbase'] not in
                            filtered_read_dict['reverse_SNV_forward_SNV1']
                        ):
                            filtered_read_dict[
                                'reverse_SNV_forward_SNV1'
                            ][dir_dict[False]['qbase']] = 1
                        else:
                            filtered_read_dict[
                                'reverse_SNV_forward_SNV1'
                            ][dir_dict[False]['qbase']] += 1

                    # Only the forward reads pass quality
                    elif dir_dict[True]['qual'] >= min_quality:
                        # Forward SNV reverse QF
                        if (
                            dir_dict[True]['qbase'] not in
                            filtered_read_dict['forward_SNV_reverse_UM_QF']
                        ):
                            filtered_read_dict[
                                'forward_SNV_reverse_UM_QF'
                            ][dir_dict[True]['qbase']] = 1
                        else:
                            filtered_read_dict[
                                'forward_SNV_reverse_UM_QF'
                            ][dir_dict[True]['qbase']] += 1

                        # Reverse QF
                        if (
                            dir_dict[False]['qbase'] not in
                            filtered_read_dict['reverse_quality_filtered']
                        ):
                            filtered_read_dict[
                                'reverse_quality_filtered'
                            ][dir_dict[False]['qbase']] = 1
                        else:
                            filtered_read_dict[
                                'reverse_quality_filtered'
                            ][dir_dict[False]['qbase']] += 1

                    # Only the reverse reads pass quality
                    elif dir_dict[False]['qual'] >= min_quality:
                        # Reverse SNV forward QF
                        if (
                            dir_dict[False]['qbase'] not in
                            filtered_read_dict['reverse_SNV_forward_UM_QF']
                        ):
                            filtered_read_dict[
                                'reverse_SNV_forward_UM_QF'
                            ][dir_dict[False]['qbase']] = 1
                        else:
                            filtered_read_dict[
                                'reverse_SNV_forward_UM_QF'
                            ][dir_dict[False]['qbase']] += 1

                        # Forward QF
                        if (
                            dir_dict[True]['qbase'] not in
                            filtered_read_dict['forward_quality_filtered']
                        ):
                            filtered_read_dict[
                                'forward_quality_filtered'
                            ][dir_dict[True]['qbase']] = 1
                        else:
                            filtered_read_dict[
                                'forward_quality_filtered'
                            ][dir_dict[True]['qbase']] += 1

                    # Neither forward nor reverse reads pass quality
                    else:
                        # Forward QF
                        if (
                            dir_dict[True]['qbase'] not in
                            filtered_read_dict['forward_quality_filtered']
                        ):
                            filtered_read_dict[
                                'forward_quality_filtered'
                            ][dir_dict[True]['qbase']] = 1
                        else:
                            filtered_read_dict[
                                'forward_quality_filtered'
                            ][dir_dict[True]['qbase']] += 1

                        # Reverse QF
                        if (
                            dir_dict[False]['qbase'] not in
                            filtered_read_dict['reverse_quality_filtered']
                        ):
                            filtered_read_dict[
                                'reverse_quality_filtered'
                            ][dir_dict[False]['qbase']] = 1
                        else:
                            filtered_read_dict[
                                'reverse_quality_filtered'
                            ][dir_dict[False]['qbase']] += 1

            # SNV in forward read only
            elif not dir_dict[True]['match'] and dir_dict[False]['match']:
                # Since only the forward read supports the SNV, quality filter
                if dir_dict[True]['qual'] >= min_quality:
                    if (
                        dir_dict[True]['qbase'] not in
                        filtered_read_dict['forward_SNV_reverse_ref']
                    ):
                        filtered_read_dict[
                            'forward_SNV_reverse_ref'
                        ][dir_dict[True]['qbase']] = 1
                    else:
                        filtered_read_dict[
                            'forward_SNV_reverse_ref'
                        ][dir_dict[True]['qbase']] += 1
                    # Since the reverse base matches the reference,
                    # don't quality filter
                    if (
                        dir_dict[False]['qbase'] not in
                        filtered_read_dict['forward_SNV_reverse_ref']
                    ):
                        filtered_read_dict[
                            'forward_SNV_reverse_ref'
                        ][dir_dict[False]['qbase']] = 1
                    else:
                        filtered_read_dict[
                            'forward_SNV_reverse_ref'
                        ][dir_dict[False]['qbase']] += 1

                # Reverse ref forward QF
                else:
                    if (
                        dir_dict[False]['qbase'] not in
                        filtered_read_dict['reverse_ref_forward_UM_QF']
                    ):
                        filtered_read_dict[
                            'reverse_ref_forward_UM_QF'
                        ][dir_dict[False]['qbase']] = 1
                    else:
                        filtered_read_dict[
                            'reverse_ref_forward_UM_QF'
                        ][dir_dict[False]['qbase']] += 1

            # SNV in reverse read only
            elif dir_dict[True]['match'] and not dir_dict[False]['match']:
                # Quality filter
                if dir_dict[False]['qual'] >= min_quality:
                    if (
                        dir_dict[False]['qbase'] not in
                        filtered_read_dict['reverse_SNV_forward_ref']
                    ):
                        filtered_read_dict[
                            'reverse_SNV_forward_ref'
                        ][dir_dict[False]['qbase']] = 1
                    else:
                        filtered_read_dict[
                            'reverse_SNV_forward_ref'
                        ][dir_dict[False]['qbase']] += 1
                    # Since the forward base matches the reference, don't
                    # quality filter
                    if (
                        dir_dict[True]['qbase'] not in
                        filtered_read_dict['reverse_SNV_forward_ref']
                    ):
                        filtered_read_dict[
                            'reverse_SNV_forward_ref'
                        ][dir_dict[True]['qbase']] = 1
                    else:
                        filtered_read_dict[
                            'reverse_SNV_forward_ref'
                        ][dir_dict[True]['qbase']] += 1

                # Forward ref reverse QF
                else:
                    if (
                        dir_dict[True]['qbase'] not in
                        filtered_read_dict['forward_ref_reverse_UM_QF']
                    ):
                        filtered_read_dict[
                            'forward_ref_reverse_UM_QF'
                        ][dir_dict[True]['qbase']] = 1
                    else:
                        filtered_read_dict[
                            'forward_ref_reverse_UM_QF'
                        ][dir_dict[True]['qbase']] += 1

            # Both reads match the reference sequence - don't quality filter,
            # and add two matches (for forward and reverse reads)
            else:
                if (
                    dir_dict[True]['qbase'] not in
                    filtered_read_dict['congruent_ref']
                ):
                    filtered_read_dict[
                        'congruent_ref'
                    ][dir_dict[True]['qbase']] = 2
                else:
                    filtered_read_dict[
                        'congruent_ref'
                    ][dir_dict[True]['qbase']] += 2

        # Either the reads are unpaired, or only a single read aligns to this
        # position on the gene (no overlap)
        else:
            for direction in dir_dict:
                # SNV supported by a single read
                if not dir_dict[direction]['match']:
                    if dir_dict[direction]['qual'] >= min_quality:
                        # Forward
                        if direction:
                            if (
                                dir_dict[direction]['qbase'] not in
                                filtered_read_dict['forward_SNV_reverse_UM_QF']
                            ):
                                filtered_read_dict[
                                    'forward_SNV_reverse_UM_QF'
                                ][dir_dict[direction]['qbase']] = 1
                            else:
                                filtered_read_dict[
                                    'forward_SNV_reverse_UM_QF'
                                ][dir_dict[direction]['qbase']] += 1

                        # Reverse
                        else:
                            if (
                                dir_dict[direction]['qbase'] not in
                                filtered_read_dict['reverse_SNV_forward_UM_QF']
                            ):
                                filtered_read_dict[
                                    'reverse_SNV_forward_UM_QF'
                                ][dir_dict[direction]['qbase']] = 1
                            else:
                                filtered_read_dict[
                                    'reverse_SNV_forward_UM_QF'
                                ][dir_dict[direction]['qbase']] += 1
                    else:
                        # Forward QF
                        if (
                            dir_dict[direction]['qbase'] not in
                            filtered_read_dict['forward_quality_filtered']
                        ):
                            filtered_read_dict[
                                'forward_quality_filtered'
                            ][dir_dict[direction]['qbase']] = 1
                        else:
                            filtered_read_dict[
                                'forward_quality_filtered'
                            ][dir_dict[direction]['qbase']] += 1

                # Match to the reference sequence supported by a single read
                else:
                    if dir_dict[direction]['qual'] >= min_quality:
                        # Forward
                        if direction:
                            if (
                                dir_dict[direction]['qbase'] not in
                                filtered_read_dict['forward_ref_reverse_UM_QF']
                            ):
                                filtered_read_dict[
                                    'forward_ref_reverse_UM_QF'
                                ][dir_dict[direction]['qbase']] = 1
                            else:
                                filtered_read_dict[
                                    'forward_ref_reverse_UM_QF'
                                ][dir_dict[direction]['qbase']] += 1

                        # Reverse
                        else:
                            if (
                                dir_dict[direction]['qbase'] not in
                                filtered_read_dict['reverse_ref_forward_UM_QF']
                            ):
                                filtered_read_dict[
                                    'reverse_ref_forward_UM_QF'
                                ][dir_dict[direction]['qbase']] = 1
                            else:
                                filtered_read_dict[
                                    'reverse_ref_forward_UM_QF'
                                ][dir_dict[direction]['qbase']] += 1
                    else:
                        # Reverse QF
                        if (
                            dir_dict[direction]['qbase'] not in
                            filtered_read_dict['reverse_quality_filtered']
                        ):
                            filtered_read_dict[
                                'reverse_quality_filtered'
                            ][dir_dict[direction]['qbase']] = 1
                        else:
                            filtered_read_dict[
                                'reverse_quality_filtered'
                            ][dir_dict[direction]['qbase']] += 1

    return filtered_read_dict, qualities, base_support


def determine_cutoff(
    *,  # Enforce keyword arguments
    qualities: List[int],
    reference_sequence: str,
    base_cutoff: int,
    error_cutoff: float = 1.0,
    max_expected_positions: float = 0.01
) -> Tuple[int, float, float]:
    """
    Determine the smallest integer k >= base_cutoff such that the probability
    of observing >= k sequencing errors by chance (given per-base qualities)
    is <= alpha, where alpha = (error_cutoff% / 100) / len(reference_sequence).

    Uses exact Poisson-binomial pmf by iterative convolution when depth is
    reasonable; otherwise uses a normal approximation with continuity
    correction.

    Args:
        qualities: Per-base quality scores at the position.
        reference_sequence: Reference sequence string (used for length/GC
        heuristics).
        base_cutoff: Base cutoff parameter (minimum bases required before SNV
        considered).
        error_cutoff: Expected error percentage (as 1.0 for 1%).

    Returns:
        A tuple (k, expected_positions, error_percent_per_site) where:
        - k is the calculated cutoff (int)
        - expected_positions is the expected number of false-positive
          positions across the gene at this cutoff (float)
        - error_percent_per_site is the per-site error probability expressed
          as a percentage (float)
    """
    # Determine the maximum read length for error percentage calculation
    max_len = max(1, len(reference_sequence))

    # Calculate alpha for error percentage calculation
    alpha = (error_cutoff / 100.0) / max_len

    # Handle edge cases: if there are no quality scores available, return
    # the minimum cutoff (at least 1) and a zero error percentage. Do NOT
    # short-circuit when base_cutoff < 1 because base_cutoff==0 indicates
    # the user requested a dynamic calculation based on observed qualities.
    if not qualities:
        return max(1, base_cutoff), 0.0, 0.0

    # Per-observation error probabilities from Phred Q
    p_list = [10 ** (-q / 10.0) for q in qualities]
    n = len(p_list)

    # Instead of modelling errors across the *entire* gene (which leads to
    # overly conservative cutoffs when depth is high), estimate a per-
    # position cutoff using the mean per-position depth. This produces a
    # cutoff k that is appropriate for a single site while still
    # accounting for sequencing quality distribution.
    mean_p = float(sum(p_list)) / n if n > 0 else 0.0
    mean_depth_per_pos = max(1, int(round(n / max_len)))

    # Build a per-position p_list approximation (identical probabilities)
    p_list_pos = [mean_p] * mean_depth_per_pos

    # Use exact poisson-binomial (FFT) when per-position depth is small,
    # otherwise normal approximation with continuity correction.
    def perpos_tail(k):
        m = mean_depth_per_pos
        if m == 0:
            return 0.0
        if m <= 1000:
            # exact via FFT convolution on the identical-p_list_pos
            pmf = _poisson_binomial_pmf_fft(p_list=p_list_pos)
            cumsum = pmf[::-1].cumsum()[::-1]
            return float(cumsum[k]) if k <= m else 0.0
        else:
            mu_pos = m * mean_p
            var_pos = m * mean_p * (1.0 - mean_p)
            sigma_pos = math.sqrt(var_pos) if var_pos > 0 else 0.0
            if sigma_pos == 0.0:
                return 1.0 if mu_pos + 1e-12 >= k else 0.0
            z = (k - 0.5 - mu_pos) / (sigma_pos * math.sqrt(2.0))
            return 0.5 * math.erfc(z)

    # per-position alpha (Bonferroni-like correction across positions)
    alpha_per_pos = alpha

    k = max(1, base_cutoff)

    # reasonable cap: mean depth + 10 sigma (if available) or
    # base_cutoff + 1000
    var_pos = mean_depth_per_pos * mean_p * (1.0 - mean_p)
    sigma_pos = math.sqrt(var_pos) if var_pos > 0 else 0.0
    cap = int(max(mean_depth_per_pos + 10 * sigma_pos, base_cutoff + 1000))

    while k <= cap and perpos_tail(k) > alpha_per_pos:
        k += 1

    # If k exceeds the achievable number of successes at a single site
    # (i.e. mean_depth_per_pos), we cannot satisfy the desired alpha with
    # the observed per-position depth. In that case, fall back to the
    # minimal practical cutoff (1 or user-supplied base_cutoff) so that
    # sites remain discoverable rather than impossible to detect.
    if k > mean_depth_per_pos:
        k = max(1, base_cutoff)

    # Optionally tighten the cutoff until the expected number of false
    # positive positions (per gene) is <= max_expected_positions. Use the
    # per-position tail probability scaled by gene length (expected count).
    expected_positions = perpos_tail(k) * max_len

    while (
        max_expected_positions is not None
        and expected_positions > max_expected_positions
        and k < mean_depth_per_pos
    ):
        k += 1
        expected_positions = perpos_tail(k) * max_len

    # If we've tightened past the per-position depth, fall back to a
    # sensible minimum cutoff
    if k > mean_depth_per_pos:
        k = max(1, base_cutoff)
        expected_positions = perpos_tail(k) * max_len

    # Safety floor for dynamic calculation (make base_cutoff==0 more
    # conservative)
    if base_cutoff == 0:
        k = max(k, MIN_DYNAMIC_CUTOFF)
        expected_positions = perpos_tail(k) * max_len

    # Provide two clear reporting metrics:
    # - expected_positions (float): expected count of false-positive
    #   positions across the gene at this cutoff
    # - error_percent_per_site (float): per-site probability expressed as %
    tail_prob = perpos_tail(k)
    error_percent_per_site = tail_prob * 100.0

    return k, expected_positions, error_percent_per_site


def _poisson_binomial_pmf_fft(
    *,  # Enforce keyword arguments
    p_list: List[float]
) -> 'np.ndarray':
    """
    Compute Poisson-binomial PMF via divide-and-conquer FFT convolution.
    Returns an array of length n+1 where pmf[j] = P(X==j).

    Args:
        p_list: List of success probabilities for independent Bernoulli trials.

    Returns:
        A numpy array containing the PMF for 0..n successes.
    """
    # convert to numpy array
    p_arr = np.asarray(p_list, dtype=float)
    n = p_arr.size
    if n == 0:
        return np.array([1.0], dtype=float)
    if n == 1:
        return np.array([1.0 - p_arr[0], p_arr[0]], dtype=float)

    # recursive divide-and-conquer
    def _rec(pv):
        m = len(pv)
        if m == 0:
            return np.array([1.0], dtype=float)
        if m == 1:
            return np.array([1.0 - pv[0], pv[0]], dtype=float)
        mid = m // 2
        left = _rec(pv[:mid])
        right = _rec(pv[mid:])
        size = left.size + right.size - 1
        nfft = 1 << ((size - 1).bit_length())
        fa = np.fft.rfft(left, nfft)
        fb = np.fft.rfft(right, nfft)
        fr = fa * fb
        conv = np.fft.irfft(fr, nfft)[:size]

        # numerical rounding cleanup
        conv[conv < 0] = 0.0
        return conv
    pmf = _rec(p_arr.tolist())

    # renormalize to sum to 1
    s = pmf.sum()
    if s > 0:
        pmf = pmf / s
    return pmf


def _beta_binomial_tail(
    *,  # Enforce keyword arguments
    k: int,
    n: int,
    a: float,
    b: float
) -> float:
    """
    Compute tail probability P(X>=k) for Beta-Binomial(n, a, b).
    Prefer SciPy's betabinom.sf when available for numerical stability;
    fall back to log-gamma summation

    Args:
        k: Minimum number of successes to include in tail.
        n: Number of trials.
        a: Alpha parameter of Beta prior.
        b: Beta parameter of Beta prior.

    Returns:
        The tail probability as a float.
    """
    # handle edge cases
    if n <= 0:
        return 0.0 if k > 0 else 1.0

    # use SciPy's betabinom.sf to calculate tail probability
    return float(betabinom.sf(k - 1, n, a, b))


def _estimate_beta_params(
    *,  # Enforce keyword arguments
    p_list: List[float]
) -> Optional[Tuple[float, float]]:
    """
    Estimate Beta distribution (alpha, beta) from probabilities list.

    Args:
        p_list: List of observed probabilities.

    Returns:
        (alpha, beta) tuple, or (None, None) if estimation is not possible
        (e.g. zero variance).
    """
    if not p_list:
        return None, None
    m = float(np.mean(p_list))
    v = float(np.var(p_list, ddof=0))
    # ensure positive variance less than m*(1-m)
    denom = m * (1.0 - m)
    if v <= 0 or v >= denom:
        return None, None
    common = (denom / v) - 1.0
    a = m * common
    b = (1.0 - m) * common
    # ensure >0
    if a <= 0 or b <= 0:
        return None, None
    return a, b


def poisson_binomial_tail(
    *,  # Enforce keyword arguments
    k: int,
    p_list: List[float]
) -> float:
    """
    Tail probability for Poisson-Binomial distribution P(X >= k). Uses
    FFT-based convolution for exact pmf when feasible (n <= 1000), otherwise
    normal approximation with continuity correction.

    Args:
        k: Threshold number of successes.
        p_list: List of success probabilities.

    Returns:
        Tail probability (float).
    """
    n = len(p_list)
    if n == 0:
        return 0.0 if k > 0 else 1.0

    # Use exact FFT convolution for moderate sizes for accuracy.
    if n <= 1000:
        pmf = _poisson_binomial_pmf_fft(p_list=p_list)
        cumsum = pmf[::-1].cumsum()[::-1]
        return float(cumsum[k]) if k <= n else 0.0

    # For large n, fall back to normal approximation with continuity
    # correction for numerical stability and performance.
    mu = sum(p_list)
    var = sum(p * (1 - p) for p in p_list)
    sigma = math.sqrt(var) if var > 0.0 else 0.0

    # Handle degenerate variance case
    if sigma == 0.0:
        return 0.0 if k > mu else 1.0

    # Survival function with continuity correction (k - 0.5)
    z = (k - 0.5 - mu) / (sigma * math.sqrt(2.0))
    p = 0.5 * math.erfc(z)

    # Clamp for numerical safety
    return float(min(1.0, max(0.0, p)))


def mann_whitney_u_p(
    *,  # Enforce keyword arguments
    x: List[float],
    y: List[float]
) -> float:
    """
    SciPy's implementation (exact when possible, otherwise normal
    approximation). If SciPy is available we delegate to it to ensure
    consistent handling of ties and exact calculations for small samples.

    Falls back to an internal normal-approximation implementation if SciPy
    is not available.

    Args:
        x, y: Lists of samples.

    Returns:
        Two-sided p-value as float.
    """
    # Initialize sample sizes
    n1 = len(x)
    n2 = len(y)

    # Handle edge cases
    n1 = len(x)
    n2 = len(y)
    if n1 == 0 or n2 == 0:
        return None

    # SciPy's mannwhitneyu supports 'two-sided' alternative; ensure we
    # return a Python float and guard against NaN results from tied ranks.
    pvalue = float(mannwhitneyu(x, y, alternative='two-sided').pvalue)
    return None if math.isnan(pvalue) else pvalue


def fisher_two_sided_p(
    *,  # Enforce keyword arguments
    a: int,
    b: int,
    c: int,
    d: int
) -> float:
    """
    Two-sided Fisher exact test p-value for a 2x2 contingency table.

    Args:
        a, b, c, d: Table counts arranged as [[a, b], [c, d]].

    Returns:
        Two-sided p-value as float.
    """
    # Hypergeometric probabilities for fixed margins
    def hyper_p(x, m1, n1, k1):
        # prob of x successes in sample size k1 from population with m1
        # successes and n1 failures
        return (
            math.comb(m1, x) * math.comb(n1, k1 - x)
        ) / math.comb(m1 + n1, k1)

    m1 = a + c  # successes in population (col1)
    n1 = b + d  # failures in population (col2)
    k1 = a + b  # sample size (row1)
    if m1 < 0 or n1 < 0 or k1 < 0:
        return 1.0
    # compute observed probability
    try:
        p_obs = hyper_p(a, m1, n1, k1)
    except ValueError:
        return 1.0
    # iterate all possible x and sum probabilities <= p_obs
    lo = max(0, k1 - n1)
    hi = min(k1, m1)
    p_sum = 0.0
    for x in range(lo, hi + 1):
        try:
            p_x = hyper_p(x, m1, n1, k1)
        except ValueError:
            p_x = 0.0
        if p_x <= p_obs + 1e-20:
            p_sum += p_x

    return min(1.0, p_sum)


def benjamini_hochberg(
    *,  # Enforce keyword arguments
    pvals: List[float]
) -> List[float]:
    """
    Compute Benjamini-Hochberg adjusted q-values.

    Args:
        pvals: List of p-values.

    Returns:
        List of q-values (same order as input pvals).
    """
    n = len(pvals)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: pvals[i])
    qvals = [0.0] * n
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        unadj = pvals[i]
        adj = min(prev, unadj * n / rank)
        qvals[i] = adj
        prev = adj
    return qvals


def combine_pvalues_fisher(
    *,  # Enforce keyword arguments
    pvals: List[float]
) -> float:
    """
    Combine a list of p-values using Fisher's method. Returns combined p-value

    Args:
        pvals: List of p-values.

    Returns:
        Combined p-value as float.
    """
    if not pvals:
        return None

    # clip p-values to avoid log(0)
    eps = 1e-300
    adj = [max(min(p, 1.0), eps) for p in pvals]
    stat = -2.0 * sum(math.log(p) for p in adj)
    k = len(adj)

    return float(chi2.sf(stat, 2 * k))


def _position_entry_passes_probabilistic_gating(
    *,
    position_stats: Dict[str, Any],
    q_threshold: float = 0.05,
    strand_p_threshold: float = 0.01,
    pos_p_threshold: float = 0.01
) -> bool:
    """
    Return True if the position has at least one alternate base that is
    statistically supported after correction for multiple testing.

    A base is considered supported if it has a BH-adjusted q-value at or
    below ``q_threshold`` and, if strand/position bias tests are available,
    those tests do not indicate significant bias.
    """
    if not position_stats:
        return False

    for stats in position_stats.values():
        q_value = stats.get('q_value')
        if q_value is None or q_value > q_threshold:
            continue

        strand_p = stats.get('strand_p')
        pos_p = stats.get('pos_p')

        if strand_p is not None and strand_p <= strand_p_threshold:
            continue
        if pos_p is not None and pos_p <= pos_p_threshold:
            continue

        return True

    return False


def find_multibase_positions(
    *,  # Enforce keyword arguments
    ref_base: str,
    filtered_read_dict: dict,
    base_cutoff: int,
    base_fraction_cutoff: float,
    base_support=None
) -> Tuple[dict, dict, int]:
    """
    Determine the characterised bases at the current position in the pileup,
    and whether any bases pass the specified cutoffs for SNV calling.

    Args:
        ref_base: String of the sequence of the reference gene at the current
        position
        filtered_read_dict: Dictionary with of base types: read direction:
        count
        base_cutoff: Integer of the number of identical mismatches in a column
        required for a SNV call
        base_fraction_cutoff: Float fraction of bases necessary to support a
        SNV call
        base_support: Optional dictionary of per-base support information

    Returns:
        snv_dict: Dictionary of characterised bases of the current column e.g.
        total, number matching reference
        sequence, number of SNVs in both forward and reverse reads, etc.
        passing_snv_dict: Dictionary summarising counts of categories of bases
        e.g. congruent, forward, reverse,
        paired
        total_coverage: Integer of the total number of bases passing filter in
        the pileup
    """
    # Initialise a dictionary to store the counts of the characterised bases
    snv_dict = {
        'total': 0,
        'total_congruent': 0,
        'total_congruent_SNV': 0,
        'total_forward': 0,
        'total_forward_SNV': 0,
        'total_reverse': 0,
        'total_reverse_SNV': 0,
        'total_SNV': 0
    }
    # Initialise the total depth of the column to zero
    total_coverage = 0

    # Initialise a dictionary to store the count for each individual
    # nucleotide at this position e.g. A:24, G:2
    base_count = {}

    # Iterate through the categories e.g. congruent_ref in the dictionary
    for category, base_dict in filtered_read_dict.items():
        # Iterate through the sequence of each query base, and the
        # corresponding count
        for base, count in base_dict.items():
            # Do not process the base if it has been flagged as 'filtered'
            if 'filtered' not in category:
                # Populate the base counting dictionary with the base sequence
                # and increment the count
                if base not in base_count:
                    base_count[base] = count
                else:
                    base_count[base] += count

                # Update the total coverage
                total_coverage += count
                snv_dict['total'] += count

    # If the reference base is not observed in the pileup, then there is no
    # allele mixture at this position. Do not call a SNV for a pure alternate
    # allele signal.
    if base_count and base_count.get(ref_base, 0) == 0:
        return snv_dict, {}, total_coverage, {}

    # Forward and reverse reads agree on base
    for category, base_dict in filtered_read_dict.items():
        # Iterate through the sequence of each query base, and the
        # corresponding count
        for base, count in base_dict.items():
            if 'filtered' not in category:
                # Forward and reverse reads agree on base
                if 'congruent' in category:
                    # Congruent reference sequence
                    snv_dict['total_congruent'] += count
                    snv_dict['total_forward'] += int(count / 2)
                    snv_dict['total_reverse'] += int(count / 2)

                    # Congruent SNVs
                    if 'SNV' in category and base != ref_base:
                        snv_dict['total_congruent_SNV'] += count
                        snv_dict['total_forward_SNV'] += int(count / 2)
                        snv_dict['total_reverse_SNV'] += int(count / 2)
                        snv_dict['total_SNV'] += count

                # SNV in forward read
                elif category.startswith('forward_SNV'):
                    snv_dict['total_forward'] += count
                    if base != ref_base:
                        snv_dict['total_forward_SNV'] += count
                        snv_dict['total_SNV'] += count

                # Forward read matches reference
                elif category.startswith('forward_ref'):
                    snv_dict['total_forward'] += count

                # SNV in reverse read
                elif category.startswith('reverse_SNV'):
                    snv_dict['total_reverse'] += count
                    if base != ref_base:
                        snv_dict['total_reverse_SNV'] += count
                        snv_dict['total_SNV'] += count

                # Reverse read match reference
                elif category.startswith('reverse_ref'):
                    snv_dict['total_reverse'] += count

    # Initialise a dictionary to store the summary of characterised base types
    passing_snv_dict = {
        'congruent': {'A': 0, 'C': 0, 'G': 0, 'T': 0},
        'forward': {'A': 0, 'C': 0, 'G': 0, 'T': 0},
        'reverse': {'A': 0, 'C': 0, 'G': 0, 'T': 0},
        'paired': {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    }

    # Boolean of whether there are bases passing filter, and the
    # passing_snv_dict should be used
    return_dict = False

    # Iterate through the categories in the dictionary
    for category, base_dict in filtered_read_dict.items():
        # Iterate through each query base and its corresponding count
        for base, count in base_dict.items():
            # Congruent SNVs
            if 'congruent' in category and base != ref_base:
                # If the base_cutoff_fraction has been provided, use it in
                # making SNV calls
                if base_fraction_cutoff:
                    # Ensure that the number of SNVs in the pileup is greater
                    # than the base_cutoff and the fraction of SNVs in the
                    # pileup is greater than the base_cutoff_fraction
                    if (
                        float(
                            base_count[base] / snv_dict['total']
                        ) >= base_fraction_cutoff and base_count[base]
                        >= base_cutoff
                    ):
                        # Update the summary dictionary with the base sequence
                        # and count
                        if base not in passing_snv_dict['forward']:
                            passing_snv_dict['congruent'][base] = count
                        else:
                            passing_snv_dict['congruent'][base] += count
                        return_dict = True

                # If base_cutoff_fraction is not supplied, only the number of
                # SNVs in the pileup must be greater than the base_cutoff
                # value in order for a SNV call
                else:
                    if base_count[base] >= base_cutoff:
                        if base not in passing_snv_dict['congruent']:
                            passing_snv_dict['congruent'][base] = count
                        else:
                            passing_snv_dict['congruent'][base] += count
                    return_dict = True

            if category.startswith('forward_SNV') and base != ref_base:
                if base_fraction_cutoff:
                    if (
                        float(
                            base_count[base] / snv_dict['total']
                        ) >= base_fraction_cutoff
                        and base_count[base] >= base_cutoff
                    ):
                        if base not in passing_snv_dict['forward']:
                            passing_snv_dict['forward'][base] = count
                        else:
                            passing_snv_dict['forward'][base] += count
                        return_dict = True
                else:
                    if base_count[base] >= base_cutoff:
                        if base not in passing_snv_dict['forward']:
                            passing_snv_dict['forward'][base] = count
                        else:
                            passing_snv_dict['forward'][base] += count
                        return_dict = True

            # SNVs in the reverse read
            if category.startswith('reverse_SNV') and base != ref_base:
                if base_fraction_cutoff:
                    if (
                        float(
                            base_count[base] / snv_dict['total']
                        ) >= base_fraction_cutoff
                        and base_count[base] >= base_cutoff
                    ):
                        if base not in passing_snv_dict['reverse']:
                            passing_snv_dict['reverse'][base] = count
                        else:
                            passing_snv_dict['reverse'][base] += count
                        return_dict = True
                else:
                    if base_count[base] >= base_cutoff:
                        if base not in passing_snv_dict['reverse']:
                            passing_snv_dict['reverse'][base] = count
                        else:
                            passing_snv_dict['reverse'][base] += count
                        return_dict = True

            # All unfiltered bases
            if base != ref_base and 'filtered' not in category:
                if base_fraction_cutoff:
                    if (
                        float(
                            base_count[base] / snv_dict['total']
                        ) >= base_fraction_cutoff
                        and base_count[base] >= base_cutoff
                    ):
                        if base not in passing_snv_dict['paired']:
                            passing_snv_dict['paired'][base] = count
                        else:
                            passing_snv_dict['paired'][base] += count
                        return_dict = True
                else:
                    if base_count[base] >= base_cutoff:
                        if base not in passing_snv_dict['paired']:
                            passing_snv_dict['paired'][base] = count
                        else:
                            passing_snv_dict['paired'][base] += count
                        return_dict = True

    # Compute per-base statistical tests if base_support was provided
    position_stats = {}

    if base_support:
        # flatten all qualities to make p_list for poisson-binomial
        # (prob of any error -> specific alt base ~ p/3)
        all_quals = []
        for b in base_support:
            all_quals.extend(base_support[b].get('quals', []))

        # Calculate p-values for each base
        p_list = [
            10 ** (-q / 10.0) / 3.0 for q in all_quals
        ] if all_quals else []

        # Iterate through each base in the base_count dictionary
        for base in base_count:
            if base == ref_base:
                continue
            k = base_count.get(base, 0)
            # p-value from poisson-binomial (exact) and mean-binomial
            # (conservative) - take the more conservative (larger) p
            if p_list:
                p_exact = poisson_binomial_tail(
                    k=k,
                    p_list=p_list
                )

                # Calculate mean probability for mean-binomial
                mean_p = sum(p_list) / len(p_list) if len(p_list) > 0 else 0.0
                p_binom = poisson_binomial_tail(
                    k=k,
                    p_list=[mean_p] * len(p_list)
                )

                # try beta-binomial if overdispersion present
                a, b = _estimate_beta_params(p_list=p_list)

                # Combine p-values
                if a is not None and b is not None:
                    p_bb = _beta_binomial_tail(
                        k=k,
                        n=len(p_list),
                        a=a,
                        b=b
                    )
                    p_value = max(p_exact, p_binom, p_bb)
                else:
                    p_value = max(p_exact, p_binom)
            else:
                p_value = None

            # strand bias
            fcount = base_support.get(base, {}).get('forward', 0)
            rcount = base_support.get(base, {}).get('reverse', 0)
            total_forward = snv_dict.get('total_forward', 0)
            total_reverse = snv_dict.get('total_reverse', 0)

            # Run Fisher's exact test
            strand_p = fisher_two_sided_p(
                a=fcount,
                b=total_forward - fcount,
                c=rcount,
                d=total_reverse - rcount
            )

            # Position bias
            base_pos = base_support.get(base, {}).get('positions', [])
            ref_pos = base_support.get(ref_base, {}).get('positions', [])

            # Mann-Whitney U test
            pos_p = mann_whitney_u_p(
                x=base_pos,
                y=ref_pos
            ) if base_pos and ref_pos else None

            # Mean qualities and mapping qualities
            mean_q_val = mean(
                base_support.get(
                    base, {}
                ).get(
                    'quals', []
                )
            ) if base_support.get(base, {}).get('quals') else None

            # Mean mapping quality
            mapqs = base_support.get(base, {}).get('mapqs', [])
            mean_mapq = mean(mapqs) if mapqs else None

            # Store the statistics for this base
            position_stats[base] = {
                'p_value': p_value,
                'strand_p': strand_p,
                'pos_p': pos_p,
                'mean_qual': mean_q_val,
                'mean_mapq': mean_mapq
            }
    if return_dict:
        return snv_dict, passing_snv_dict, total_coverage, position_stats

    return snv_dict, {}, total_coverage, position_stats


def position_details(
    *,  # Enforce keyword arguments
    actual_position: int,
    passing_snv_dict: Dict[str, int],
    contig_name: str,
    ref_base: str,
    total_coverage: int,
    base_cutoff: int,
    error_perc: float,
    p_value: Optional[float] = None,
    adj_p_value: Optional[float] = None,
    strand_p: Optional[float] = None,
    pos_p: Optional[float] = None,
    mean_q: Optional[float] = None,
    mean_mapq: Optional[float] = None
) -> Dict[str, Any]:
    """
    Format per-position statistics for reporting and return a dict
    representing the row.

    Args:
        actual_position: 1-based position in reference sequence.
        passing_snv_dict: Mapping base->coverage for bases that passed filters.
        contig_name: Contig name.
        ref_base: Reference base at this position.
        total_coverage: Total depth of coverage at the position.
        base_cutoff: Base cutoff used for calling.
        error_perc: Error percentage used in calculations.
        p_value, adj_p_value, strand_p, pos_p: Optional p-values calculated
        elsewhere.
        mean_q, mean_mapq: Optional mean base and mapping qualities.

    Returns:
        A dictionary with keys matching the TSV columns used by the pipeline.
    """
    # List of the base categories present in passing_snv_dict
    read_types = ['congruent', 'paired', 'forward', 'reverse']

    # Update the to_write string with basic information
    to_write = f'{contig_name}\t{actual_position}\t{ref_base}\t'

    # Initialise the coverage value to zero
    snv_coverage = 0

    # Iterate through all the read types
    for read_type in read_types:
        # Boolean of whether a semi-colon needs to be added to the string to
        # separate all the base sequences
        semi_colon = False
        for base, coverage in passing_snv_dict[read_type].items():
            # Ensure that the base isn't empty
            if base:
                if semi_colon:
                    to_write += ';'
                to_write += f'{base}:{coverage}'
                semi_colon = True

                # Increment the coverage only if the read type isn't 'paired'
                if read_type != 'paired':
                    snv_coverage += coverage
        to_write += ','

    # Format the error_perc to be more readable
    error_perc = f'{error_perc:0.2f}' if error_perc else 'ND'

    # Format the statistical fields
    p_str = f'{p_value:0.3e}' if p_value is not None else 'ND'
    q_str = f'{adj_p_value:0.3e}' if adj_p_value is not None else 'ND'
    strand_str = f'{strand_p:0.3e}' if strand_p is not None else 'ND'
    pos_str = f'{pos_p:0.3e}' if pos_p is not None else 'ND'
    mean_q_str = f'{mean_q:0.2f}' if mean_q is not None else 'ND'
    mean_mapq_str = f'{mean_mapq:0.2f}' if mean_mapq is not None else 'ND'

    # Finalise the to_write string with the remaining statistics
    to_write += (
        f'{snv_coverage}\t{total_coverage}\t{base_cutoff}\t{error_perc}\t'
        f'{p_str}\t{q_str}\t{strand_str}\t{pos_str}\t{mean_q_str}\t'
        f'{mean_mapq_str}\n'
    )
    return to_write


def read_contig(
    *,  # Enforce keyword arguments
    contig_name: str,
    bamfile_name: str,
    reference_fasta: str,
    allele_records: Optional[Any] = None,
    fastq_records: Optional[Dict[str, Any]] = None,
    quality_cutoff: int = 20,
    min_quality: int = 15,
    base_cutoff: Optional[int] = None,
    base_fraction_cutoff: Optional[float] = None,
    fasta: bool = False,
    error_cutoff: float = 1.0,
    nanopore: bool = False,
    max_expected_positions: float = 0.001,
    start: Optional[int] = None,
    end: Optional[int] = None
) -> Tuple[Dict[str, Any], str]:

    """
    Analyse a contig in a BAM file and return multibase dict and TSV text.

    Args:
        contig_name: Name of contig to analyse.
        bamfile_name: Path to BAM file.
        reference_fasta: Path to reference FASTA.
        allele_records: SeqIO index or similar mapping of alleles.
        fastq_records: Mapping of read_name -> SeqRecord with quality info.

    Keyword Args:
        quality_cutoff: Minimum base quality used only for read trimming and
            earlier filtering (default: 20).
        min_quality: Minimum base quality required to count a base as SNV
            support (default: 15).
        base_cutoff: Absolute base-count cutoff (optional).
        base_fraction_cutoff: Fractional cutoff of coverage (optional).
        fasta: If True, operate in FASTA mode (no base qualities).
        error_cutoff: Error cutoff used when computing base_cutoff adjustments.
        nanopore: If True, enable nanopore-specific handling.

    Returns:
        A tuple (multibase_dict, tsv_output) where multibase_dict is a dict of
        detected multibase positions and gene summary stats, and tsv_output is
        the position-level TSV text for the contig.
    """
    # Initialize pysam FastaFile object and other variables
    pysam_fasta = pysam.FastaFile(reference_fasta)

    # Ensure contig_name is str (some systems/pysam versions may provide bytes)
    if isinstance(contig_name, bytes):
        try:
            contig_name = contig_name.decode('utf-8')
        except UnicodeDecodeError:
            contig_name = contig_name.decode('utf-8', errors='replace')

    # Quick sanity-check: ensure the contig exists in the reference FASTA. If
    # not present, warn and skip analysis for this contig to avoid faidx/pysam
    # errors during pileup. Use SeqIO as a fallback if pysam doesn't list refs.
    try:
        fasta_refs = set(
            x.decode() if isinstance(
                x,
                bytes
            ) else str(x) for x in pysam_fasta.references
        )

        # Check if contig_name is in the FASTA references
        if contig_name not in fasta_refs:
            # Fallback: parse FASTA with SeqIO to double-check presence (works
            # even if pysam's internal references are not populated).
            try:
                seqio_ids = set(
                    rec.id for rec in SeqIO.parse(reference_fasta, 'fasta')
                )
                if contig_name in seqio_ids:
                    logging.debug(
                        'Contig %s found via SeqIO.parse but not present in '
                        'pysam.FastaFile.references; attempting to build .fai '
                        'index and reopen.',
                        contig_name
                    )
                    # Attempt to create a FASTA index if it doesn't exist
                    try:
                        if not os.path.isfile(reference_fasta + '.fai'):

                            logging.debug(
                                'Creating FASTA index for %s', reference_fasta
                            )

                            # Create the FASTA index
                            pysam.faidx(reference_fasta)

                        # Reopen to refresh references
                        pysam_fasta.close()

                        # Reopen the FastaFile
                        pysam_fasta = pysam.FastaFile(reference_fasta)

                        # Recheck references
                        fasta_refs = set(
                            x.decode() if isinstance(
                                x,
                                bytes
                            ) else str(x) for x in pysam_fasta.references
                        )

                        # Final check for contig presence
                        if contig_name not in fasta_refs:
                            logging.warning(
                                'Contig %s present in FASTA but did not '
                                'appear after index creation; skipping.',
                                contig_name
                            )
                            # Close pysam_fasta before returning
                            pysam_fasta.close()

                            return {}, ''
                    except (SamtoolsError, OSError, ValueError) as exc_idx:
                        logging.debug(
                            'Failed to create/reopen FASTA index: %s',
                            exc_idx
                        )
                        logging.warning(
                            'Contig %s not accessible via pysam for FASTA %s; '
                            'skipping this contig.',
                            contig_name, reference_fasta
                        )

                        # Close pysam_fasta before returning
                        try:
                            pysam_fasta.close()
                        except (OSError, AttributeError):
                            # Best-effort close, ignore inability to close
                            pass
                        return {}, ''
                else:
                    logging.warning(
                        'Contig %s not found in reference FASTA %s; skipping '
                        'this contig.',
                        contig_name, reference_fasta
                    )
                    pysam_fasta.close()
                    return {}, ''
            except (OSError, ValueError) as exc2:
                logging.debug('SeqIO parse check failed: %s', exc2)
                logging.warning(
                    'Contig %s not found in reference FASTA %s; skipping '
                    'this contig.',
                    contig_name, reference_fasta
                )
                pysam_fasta.close()
                return {}, ''
    except (
        SamtoolsError,
        AttributeError,
        TypeError,
        OSError,
        ValueError,
        UnicodeDecodeError
    ) as exc:
        logging.debug(
            'Could not verify contig presence in FASTA via pysam: %s', exc
        )

    # Initialize variables
    multibase_position_dict = {}
    to_write = str()

    # If analysing FASTA files, a single base difference is all that
    # is expected
    if fasta:
        base_cutoff = 1

    # Extract the reference sequence for the contig being analysed.
    # If allele_records dict wasn't passed (to avoid pickling large dicts),
    # load only the requested contig from the FASTA file.
    if allele_records is not None:
        reference_sequence = str(allele_records[contig_name].seq)
    else:
        try:
            idx = SeqIO.index(reference_fasta, 'fasta')
            # Avoid KeyError by checking membership first
            if contig_name in idx:
                reference_sequence = str(idx[contig_name].seq)
            else:
                # Fall back to parsing whole FASTA if contig not in index
                allele_tmp = SeqIO.to_dict(
                    SeqIO.parse(reference_fasta, 'fasta')
                )
                reference_sequence = str(allele_tmp[contig_name].seq)
        except (FileNotFoundError, OSError, ValueError) as exc:
            # Indexing failed due to I/O or parse error; try full-parse
            logging.debug(
                'SeqIO.index failed for %s: %s',
                reference_fasta,
                exc
            )
            allele_tmp = SeqIO.to_dict(
                SeqIO.parse(reference_fasta, 'fasta')
            )
            reference_sequence = str(allele_tmp[contig_name].seq)
        except KeyError:
            # Contig not found even after fallback
            logging.warning(
                'Contig %s not found in reference FASTA %s',
                contig_name, reference_fasta
            )
            return {}, ''

    # Parse the BAM file with pysam to create AlignmentFile, and
    # AlignmentFile.pileup objects (support optional region start/end)
    bamfile, pileup = parse_bam(
        bamfile_name=bamfile_name,
        contig_name=contig_name,
        pysam_fasta=pysam_fasta,
        start=start,
        end=end
    )

    # Define dictionaries and lists to store filtered reads, base support and
    # quality scores
    filtered_read_dict = {}
    base_support_dict = {}
    quality_list = []

    # Containers for report aggregation and multiple-testing correction
    report_entries = []
    all_tests = []

    # Iterate through each column in the pileup
    for i, column in enumerate(pileup):
        filtered_reads, qualities, base_support = characterise_read(
            column=column,
            reference_sequence=reference_sequence,
            fastq_records=fastq_records,
            quality_cutoff=quality_cutoff,
            min_quality=min_quality,
            fasta=fasta,
            nanopore=nanopore
        )

        # Populate the dictionaries and lists with the results from
        # characterise_read
        filtered_read_dict[i] = filtered_reads
        base_support_dict[i] = base_support
        quality_list += qualities

    # Initialise the calculated error percentage to zero
    error_perc = None

    # Determine the maximum read length for error percentage calculation
    max_len = max(1, len(reference_sequence))

    # If the base_cutoff is set to zero, determine the appropriate cutoff value
    if base_cutoff == 0:

        # Set the max expected positions
        max_expected_positions = \
            max_expected_positions if 'max_expected_positions' \
            in locals() else 0.001
        try:
            computed = determine_cutoff(
                qualities=quality_list,
                reference_sequence=reference_sequence,
                base_cutoff=base_cutoff,
                error_cutoff=error_cutoff,
                max_expected_positions=max_expected_positions
            )
            # Determine_cutoff now returns (k, expected_positions,
            # error_percent_per_site)
            if isinstance(computed, tuple) and len(computed) >= 3:
                computed_cutoff = int(computed[0])
                expected_positions = float(computed[1])
                error_percent_per_site = float(computed[2])
                error_perc = error_percent_per_site
            elif isinstance(computed, tuple) and len(computed) == 2:
                computed_cutoff = int(computed[0])
                expected_positions = float(computed[1])
                error_perc = (expected_positions / max_len) * 100.0
            else:
                computed_cutoff = int(computed)
                expected_positions = 0.0
                error_perc = 0.0

            base_cutoff = computed_cutoff
            logging.debug(
                'Contig %s: computed dynamic base_cutoff=%s, '
                'expected_positions=%0.3f, per_site_percent=%0.6f',
                contig_name, base_cutoff, expected_positions, error_perc
            )
        except (ValueError, TypeError):
            logging.debug(
                'Contig %s: error computing dynamic base cutoff: %s',
                contig_name, traceback.format_exc()
            )
    bamfile.close()

    # It seems that the pileup (generator?) is used up above, so it must be
    # recreated (support optional region start/end)
    bamfile, pileup = parse_bam(
        bamfile_name=bamfile_name,
        contig_name=contig_name,
        pysam_fasta=pysam_fasta,
        start=start,
        end=end
    )

    # Iterate through each column in the pileup
    for i, column in enumerate(pileup):
        # Extract the sequence of the reference gene at the current position
        ref_base = reference_sequence[column.pos]

        # Summarise the pileup (now collecting statistics per position)
        _, passing_snv_dict, total_coverage, position_stats = \
            find_multibase_positions(
                ref_base=ref_base,
                filtered_read_dict=filtered_read_dict[i],
                base_cutoff=base_cutoff,
                base_fraction_cutoff=base_fraction_cutoff,
                base_support=base_support_dict.get(i, {})
            )

        # If there are any SNVs called for the gene, aggregate stats and
        # record the position for downstream filtering after multiple-test
        # correction.
        if passing_snv_dict:
            # Pysam starts counting at 0, whereas we actually want to start
            # counting at 1.
            actual_position = column.pos + 1

            # Store entry for per-gene reporting and FDR correction
            report_entries.append({
                'gene': column.reference_name,
                'position': actual_position,
                'ref_base': ref_base,
                'passing_snv_dict': passing_snv_dict,
                'total_coverage': total_coverage,
                'base_cutoff': base_cutoff,
                'error_perc': error_perc,
                'position_stats': position_stats
            })

            # Collect p-values for BH adjustment (per-base tests)
            for base, stats in (position_stats or {}).items():
                pval = stats.get('p_value')
                if pval is not None:
                    all_tests.append(
                        {
                            'entry_idx': len(report_entries) - 1,
                            'base': base, 'p': pval
                        }
                    )

    # Multiple-testing correction across bases in this gene
    if all_tests:
        pvals = [t['p'] for t in all_tests]
        qvals = benjamini_hochberg(pvals=pvals)
        for idx, t in enumerate(all_tests):
            t['q'] = qvals[idx]
            entry = report_entries[t['entry_idx']]
            base = t['base']

            # Attach test summary to entry
            if 'tests' not in entry:
                entry['tests'] = []
            entry['tests'].append({'base': base, 'p': t['p'], 'q': t['q']})

        # Annotate q-values back into position_stats
        for entry in report_entries:
            pos_stats = entry.get('position_stats') or {}
            for test in entry.get('tests', []):
                b = test['base']
                if b in pos_stats:
                    pos_stats[b]['q_value'] = test['q']

    # Filter positions using probabilistic evidence. If per-base statistics are
    # available, require at least one alternate base to be statistically
    # supported before the position is retained.
    valid_report_entries = []
    for entry in report_entries:
        position_stats = entry.get('position_stats') or {}
        if not position_stats or _position_entry_passes_probabilistic_gating(
            position_stats=position_stats
        ):
            valid_report_entries.append(entry)

    # Rebuild the multibase dict only from valid positions.
    multibase_position_dict = {}
    for entry in valid_report_entries:
        gene = entry['gene']
        pos = entry['position']
        if gene not in multibase_position_dict:
            multibase_position_dict[gene] = {}
        multibase_position_dict[gene][pos] = entry['passing_snv_dict']

    # Compute per-gene combined p-value (Fisher) and a gene-level score
    # (sum -log10(q)) p-values used for Fisher should be the unadjusted
    # position p-values; q-values are used for scoring.
    pvals_for_gene = [
        test['p']
        for entry in valid_report_entries
        for test in entry.get('tests', [])
    ] if valid_report_entries else []
    qvals_for_gene = [
        test.get('q')
        for entry in valid_report_entries
        for test in entry.get('tests', [])
        if 'q' in test
    ] if valid_report_entries else []

    # Calculate combined p-value and gene score
    combined_p = combine_pvalues_fisher(
        pvals=pvals_for_gene
    ) if pvals_for_gene else None

    # Initialize gene_score and num_sig_positions
    gene_score = None

    # Calculate gene_score and num_sig_positions if q-values are available
    if qvals_for_gene:
        gene_score = sum([-math.log10(max(q, 1e-300)) for q in qvals_for_gene])

    # Calculate number of significant positions based on q-value threshold
    num_sig_positions = sum(
        1 for q in qvals_for_gene if q is not None and q <= 0.05
    ) if qvals_for_gene else 0

    # Store gene-level stats in the multibase_position_dict under a
    # special key for this gene
    if contig_name not in multibase_position_dict:
        multibase_position_dict[contig_name] = {}
    multibase_position_dict[contig_name]['_gene_stats'] = {
        'combined_p': combined_p,
        'gene_score': gene_score,
        'num_sig_positions': num_sig_positions
    }

    # Build a per-gene summary TSV file for this sample. Aggregate across all
    # per-gene dicts returned later. (Actual file will be written once in
    # find_contamination after collecting all genes.) build output lines
    to_write = ''
    for entry in valid_report_entries:
        pos_stats = entry.get('position_stats') or {}
        if pos_stats:
            # Choose the base with smallest p-value for summary metrics
            min_base = min(
                pos_stats.items(), key=lambda kv: kv[1].get('p_value', 1.0)
            )[0]
            stats = pos_stats[min_base]
            pval = stats.get('p_value')
            qval = stats.get('q_value')
            strand_p = stats.get('strand_p')
            pos_p = stats.get('pos_p')
            mean_q = stats.get('mean_qual')
            mean_mapq = stats.get('mean_mapq')
        else:
            pval = qval = strand_p = pos_p = mean_q = mean_mapq = None
        to_write += position_details(
            actual_position=entry['position'],
            passing_snv_dict=entry['passing_snv_dict'],
            contig_name=entry['gene'],
            ref_base=entry['ref_base'],
            total_coverage=entry['total_coverage'],
            base_cutoff=entry['base_cutoff'],
            error_perc=entry['error_perc'],
            p_value=pval,
            adj_p_value=qval,
            strand_p=strand_p,
            pos_p=pos_p,
            mean_q=mean_q,
            mean_mapq=mean_mapq
        )
    bamfile.close()

    return multibase_position_dict, to_write


def _build_contig_chunks(
    contig_lengths: Dict[str, int],
    contigs: List[str],
    threads: int,
    multiplier: int = CONTIG_CHUNK_MULTIPLIER,
    max_chunk_bases: int = CONTIG_CHUNK_MAX_BASES
) -> List[List[Dict[str, Optional[int]]]]:
    """
    Build balanced chunks of contigs (first-fit decreasing) splitting large
    contigs into subranges of size <= max_chunk_bases.

    Returns a list of chunks where each chunk is a list of items of the
    form {'contig': name, 'start': start_or_None, 'end': end_or_None,
    'length': length}.
    """
    # Build items (split large contigs)
    items = []
    for c in contigs:
        length = contig_lengths.get(c)
        if length is None:
            # Fallback length 1 to avoid crashes; calling code should log
            # and recover when contig is missing.
            length = 1

        # If the contig is longer than max_chunk_bases, split it
        if length > max_chunk_bases:
            # Split into subranges
            for s in range(0, length, max_chunk_bases):
                # Define end position as the minimum of the chunk size or
                # the contig length
                e = min(s + max_chunk_bases, length)

                # Append the subrange item
                items.append(
                    {'contig': c, 'start': s, 'end': e, 'length': e - s}
                )
        # Otherwise, add the whole contig as a single item
        else:
            items.append(
                {'contig': c, 'start': None, 'end': None, 'length': length}
            )

    # Sort descending by length
    items.sort(key=lambda x: x['length'], reverse=True)

    # Initialize chunks as the greatest of (threads * multiplier) or 1
    n_chunks = max(1, threads * max(1, multiplier))

    # Initialize empty chunks
    chunks = [{'total': 0, 'items': []} for _ in range(n_chunks)]

    # First-fit decreasing: put each item to the chunk with smallest total
    for it in items:
        # Find chunk with minimum total
        idx = min(range(len(chunks)), key=lambda i: chunks[i]['total'])
        chunks[idx]['items'].append(it)
        chunks[idx]['total'] += it['length']

    # Return only the list of items per chunk
    return [c['items'] for c in chunks if c['items']]


def _read_contig_chunk_dispatch(
    kwargs: Dict[str, Any]
) -> Tuple[Dict[str, Any], str]:
    """
    Worker helper that processes a chunk of contigs/ranges.

    Args:
        kwargs: Dictionary of keyword arguments including 'contig_chunk', which
        is a list of dicts with keys 'contig', 'start', 'end'.

    Returns:
        Tuple of (combined_multibase_dict, combined_report_text).
    """
    # Extract contig_chunk from kwargs
    chunk = kwargs.get('contig_chunk') or []

    # Get the PID for logging
    pid = os.getpid()

    logging.debug('Worker %s: starting chunk with %s items', pid, len(chunk))

    # Start timing
    t0 = time.time()

    # Initialize combined results
    combined = {}
    combined_reports = []
    try:
        # Process each item in the chunk
        for item in chunk:
            # Extract contig, start, end
            contig = item.get('contig')
            start = item.get('start')
            end = item.get('end')

            # Build per-contig kwargs for read_contig
            per_kwargs = kwargs.copy()
            per_kwargs.pop('contig_chunk', None)
            per_kwargs['contig_name'] = contig

            # Set start/end if present
            if start is not None:
                per_kwargs['start'] = start

            # Set end if present
            if end is not None:
                per_kwargs['end'] = end

            # Call read_contig
            multibase_dict, report_write = read_contig(**per_kwargs)

            # Merge multibase_dicts
            for k, v in multibase_dict.items():
                combined.setdefault(k, {}).update(v)

            # Append report text if present
            if report_write:
                combined_reports.append(report_write)

        # Calculate elapsed time
        elapsed = time.time() - t0

        logging.debug('Worker %s: finished chunk (%.2fs)', pid, elapsed)

        return combined, '\n'.join(combined_reports)
    except (
        SamtoolsError,
        OSError,
        ValueError,
        KeyError,
        IndexError,
        TypeError,
        RuntimeError
    ) as exc:
        # Log expected worker-level errors and re-raise
        logging.exception('Worker %s: error processing chunk: %s', pid, exc)
        raise


def _read_contig_dispatch(
    *,  # Enforce keyword arguments
    kwargs: Dict[str, Any]
) -> Tuple[Dict[str, Any], str]:
    """
    Helper for multiprocessing that calls `read_contig` with kwargs.

    multiprocessing.Pool.map passes a single argument to the worker
    function. Since ``read_contig`` enforces keyword-only arguments, we
    build a dict of kwargs and dispatch via this helper.

    Args:
        kwargs: Dictionary of keyword arguments for `read_contig`.

    Returns:
        Tuple of (multibase_dict, report_text) returned by `read_contig`.
    """
    # Extract contig name for logging
    contig = kwargs.get('contig_name')

    # Get the PID for logging
    pid = os.getpid()

    logging.debug('Worker %s: starting contig %s', pid, contig)

    # Start timing
    t0 = time.time()
    try:
        # Defensive checks: ensure kwargs is a dict and has the required key(s)
        if not isinstance(kwargs, dict):
            raise TypeError('Expected kwargs to be a dict')
        if 'contig_name' not in kwargs:
            raise TypeError("Missing required keyword 'contig_name'")

        # Call read_contig with the provided kwargs
        result = read_contig(**kwargs)
        elapsed = time.time() - t0

        logging.debug(
            'Worker %s: finished contig %s (%.2fs)', pid, contig, elapsed
        )

        return result
    except TypeError as exc:
        logging.exception(
            "Worker %s: invalid arguments for read_contig %s: %s",
            pid, contig, exc
        )
        raise
    except (
        SamtoolsError,
        OSError,
        ValueError,
        KeyError,
        IndexError,
        RuntimeError
    ) as exc:
        logging.exception(
            "Worker %s: error processing contig %s: %s", pid, contig, exc
        )
        raise


def count_multibase_positions(
    *,  # Enforce keyword arguments
    multibase_dict_list: List[Dict[str, Any]]
) -> int:
    """
    Count multibase positions across the list of per-gene multibase dicts.

    The per-gene dicts have the shape {gene: {position: {...}, '_gene_stats':
    {...}}}.

    This function excludes meta keys that start with '_' (e.g. '_gene_stats').

    Args:
        multibase_dict_list: List of dicts returned by `read_contig` for each
        gene.

    Returns:
        Total number of positions (int).
    """
    # Initialize total to zero
    total = 0

    # Iterate through each gene's multibase dict
    for multibase_position_dict in multibase_dict_list:
        # Iterate through each gene and its SNP positions
        for _, snp_positions in multibase_position_dict.items():
            # Increment total by the number of positions excluding meta keys
            total += sum(
                1 for k in snp_positions.keys() if not str(k).startswith('_')
            )

    return total


def find_rmlst_type(
    *,  # Enforce keyword arguments
    kma_report: str, rmlst_report: str
) -> List[str]:
    """
    Uses a report generated by KMA to determine what allele is present for
    each rMLST gene.

    Args:
        kma_report: Path to KMA report file.
        rmlst_report: Path to output rMLST report file.

    Returns:
        List of strings representing the gene and allele called e.g.
        ['abcZ_1', 'adk_4', ...]
    """
    # Initialize a dictionary to store the best scoring allele for each gene
    genes_to_use = {}

    # Initialize a dictionary to store the highest score for each gene
    score_dict = {}

    # Initialize a list to store the final gene_allele strings
    gene_alleles = []

    # Read through the KMA report file
    with open(kma_report, encoding='utf-8') as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')

        # Iterate through each row in the KMA report
        for row in reader:
            # Extract the gene_allele and score from the row
            gene_allele = row['#Template']

            # Extract the score from the row
            score = int(row['Score'])

            # Split the gene_allele into gene and allele components
            gene = gene_allele.split('_')[0]
            allele = gene_allele.split('_')[1]

            # Determine if this allele has the highest score for the gene
            if gene not in score_dict:
                score_dict[gene] = score
                genes_to_use[gene] = allele
            else:
                if score > score_dict[gene]:
                    score_dict[gene] = score
                    genes_to_use[gene] = allele

    # Create the gene_allele strings and write to the rMLST report file
    for gene, allele in genes_to_use.items():
        gene_alleles.append(gene + '_' + allele.replace(' ', ''))

    # Sort the gene_alleles list before writing to file
    gene_alleles = sorted(gene_alleles)

    # Write the rMLST report file
    with open(rmlst_report, 'w', encoding='utf-8') as f:
        # rMLST report is a TSV: Gene\tAllele
        f.write('Gene\tAllele\n')

        # Write each gene and allele to the report file
        for gene_allele in gene_alleles:
            gene = gene_allele.split('_')[0]
            allele = gene_allele.split('_')[1]
            f.write(f'{gene}\t{allele}\n')
    return gene_alleles


def base_dict_to_string(
    *,  # Enforce keyword arguments
    base_dict: Dict[str, int]
) -> str:
    """
    Converts a dictionary to a string. {'C': 12, 'A':4} gets converted to
    C:12;A:4

    Args:
        base_dict: Dictionary of bases and counts created by find_if_multibase

    Returns:
        A string representation of the base_dict
    """
    # Initialize the output string
    outstr = ''

    # First, sort base_dict so that major allele always comes first -
    # makes output report nicer to look at.
    base_list = sorted(base_dict.items(), key=lambda kv: kv[1], reverse=True)
    for base in base_list:
        outstr += f'{base[0]}:{base[1]};'

    return outstr[:-1]


def find_total_sequence_length(
    *,  # Enforce keyword arguments
    fasta_file: str
) -> int:
    """
    Totals up number of bases in a fasta file.

    Args:
        fasta_file: Path to FASTA file.

    Returns:
        Total number of bases in the FASTA file.
    """
    # Initialize total_length to zero
    total_length = 0

    # Iterate through each sequence in the FASTA file and sum lengths
    for sequence in SeqIO.parse(fasta_file, 'fasta'):
        total_length += len(sequence.seq)

    return total_length


def load_fastq_records(
    *,  # Enforce keyword arguments
    gz: Any,
    paired: bool,
    forward: bool
) -> Dict[str, Any]:
    """
    Use SeqIO to load FASTQ records from file

    Args:
        gz: Path to FASTQ file (can be gzipped) or open file handle.
        paired: Boolean of whether reads are paired.
        forward: Boolean of whether reads are forward reads.

    Returns:
        Dictionary of SeqIO records keyed by read ID.
    """
    # Initialise a dictionary to store the FASTQ records
    records = {}

    close_handle = False
    if isinstance(gz, str) and gz.endswith('.gz'):
        handle = gzip.open(gz, 'rt', encoding='utf-8', errors='ignore')
        close_handle = True
    else:
        handle = gz

    try:
        for record in SeqIO.parse(handle, 'fastq'):
            # Only update the naming scheme for paired reads
            if paired:
                if forward:
                    # Don't worry if the record.id already has a /1
                    if not record.id.endswith('/1'):
                        record.id = record.id + '/1'
                # Process reverse reads in a similar fashion to forward reads
                else:
                    if not record.id.endswith('/2'):
                        record.id = record.id + '/2'
            records.update(SeqIO.to_dict([record]))
    finally:
        if close_handle:
            handle.close()
    return records


def _init_fastq_index(
    *,  # Enforce keyword arguments
    fwd_path: Optional[str],
    rev_path: Optional[str],
    paired: bool,
    tmpdir: str
) -> None:
    """
    Pool initializer to create per-worker on-disk FASTQ indexes using
    Bio.SeqIO.index_db. Index objects are stored in module-level globals
    so worker processes can access reads by ID without large in-memory
    dicts being pickled to each worker.

    Args:
        fwd_path: Path to forward (R1) FASTQ file (may be gzipped).
        rev_path: Path to reverse (R2) FASTQ file (may be gzipped) or None.
        paired: Whether reads are paired.
        tmpdir: Directory where index DB files will be created.

    Returns:
        None
    """
    # Use the mutable state dict instead of module-level globals
    _FASTQ_INDEX_STATE['paired'] = bool(paired)

    try:
        # Check if forward FASTQ file exists
        if fwd_path and os.path.exists(fwd_path):
            # Define path for forward index DB
            fwd_db = os.path.join(
                tmpdir,
                os.path.basename(fwd_path) + '.fwd.idx'
            )

            # Build or open the index DB for the forward reads
            _FASTQ_INDEX_STATE['fwd'] = SeqIO.index_db(
                fwd_db,
                fwd_path,
                'fastq'
            )
        else:
            _FASTQ_INDEX_STATE['fwd'] = None
    except (OSError, ValueError, RuntimeError) as exc:
        # More specific exceptions for I/O and SeqIO failures
        logging.warning(
            'Failed to build forward FASTQ index (%s): %s', fwd_path, exc
        )
        _FASTQ_INDEX_STATE['fwd'] = None

    # Build reverse index if paired
    if paired and rev_path:
        try:
            # Check if reverse FASTQ file exists
            if rev_path and os.path.exists(rev_path):
                # Define path for reverse index DB
                rev_db = os.path.join(
                    tmpdir,
                    os.path.basename(rev_path) + '.rev.idx'
                )
                # Set up the index DB for the reverse reads
                _FASTQ_INDEX_STATE['rev'] = SeqIO.index_db(
                    rev_db,
                    rev_path,
                    'fastq'
                )
            else:
                _FASTQ_INDEX_STATE['rev'] = None
        except (OSError, ValueError, RuntimeError) as exc:
            logging.warning(
                'Failed to build reverse FASTQ index (%s): %s', rev_path, exc
            )
            _FASTQ_INDEX_STATE['rev'] = None
    else:
        _FASTQ_INDEX_STATE['rev'] = None


def _get_fastq_record(
    read_name: str
):
    """
    Retrieve a SeqRecord for `read_name` using per-worker indexes if
    available. Returns None if not found.
    """
    # Read indexes from the shared mutable dict (no 'global' needed)
    idx_fwd = _FASTQ_INDEX_STATE.get('fwd')
    idx_rev = _FASTQ_INDEX_STATE.get('rev')
    paired = _FASTQ_INDEX_STATE.get('paired', False)

    if idx_fwd is None and idx_rev is None:
        return None

    # Normalize read_name (strip any trailing description fields)
    base_name = read_name.split(' ')[0]

    # Try exact lookup in forward index first
    try:
        if idx_fwd is not None and base_name in idx_fwd:
            return idx_fwd[base_name]
    except (KeyError, TypeError, AttributeError, IndexError):
        # Handle lookup/index implementations that raise different errors
        pass

    # If paired, try the appropriate index based on /1 or /2 suffix
    if paired:
        if base_name.endswith('/1'):
            try:
                if idx_fwd is not None and base_name in idx_fwd:
                    return idx_fwd[base_name]
            except (KeyError, TypeError, AttributeError, IndexError):
                pass
            # Try without suffix
            try:
                key2 = base_name[:-2]
                if idx_fwd is not None and key2 in idx_fwd:
                    return idx_fwd[key2]
            except (KeyError, TypeError, AttributeError, IndexError):
                pass
        if base_name.endswith('/2'):
            try:
                if idx_rev is not None and base_name in idx_rev:
                    return idx_rev[base_name]
            except (KeyError, TypeError, AttributeError, IndexError):
                pass
            try:
                key2 = base_name[:-2]
                if idx_rev is not None and key2 in idx_rev:
                    return idx_rev[key2]
            except (KeyError, TypeError, AttributeError, IndexError):
                pass

    # As a last resort, check both indexes for the base name without suffix
    try:
        if idx_fwd is not None and base_name in idx_fwd:
            return idx_fwd[base_name]
    except (KeyError, TypeError, AttributeError, IndexError):
        pass
    try:
        if idx_rev is not None and base_name in idx_rev:
            return idx_rev[base_name]
    except (KeyError, TypeError, AttributeError, IndexError):
        pass

    return None


def _ensure_fastq_index(
    fwd_path: Optional[str],
    rev_path: Optional[str],
    paired: bool,
    tmpdir: str
) -> None:
    """
    Ensure FASTQ index DB files exist in tmpdir. This is intended to be
    called once by the master process before spawning worker processes so
    workers can quickly open existing indexes instead of racing to create
    them concurrently.

    Args:
        fwd_path: Path to forward FASTQ file (may be gzipped).
        rev_path: Path to reverse FASTQ file (may be gzipped) or None.
        paired: Whether reads are paired.
        tmpdir: Directory where index DB files will be created.

    Returns:
        None
    """
    try:
        # Check if forward FASTQ file exists
        if fwd_path and os.path.exists(fwd_path):
            # Define path for forward index DB
            fwd_db = os.path.join(
                tmpdir, os.path.basename(fwd_path) + '.fwd.idx'
            )

            # Build the index DB for the forward reads if it doesn't exist
            if not os.path.exists(fwd_db):
                logging.debug('Building forward FASTQ index at %s', fwd_db)
                SeqIO.index_db(fwd_db, fwd_path, 'fastq')

    except (OSError, ValueError, RuntimeError) as exc:
        logging.warning(
            'Failed to build forward FASTQ index (%s): %s', fwd_path, exc
        )

    # Build reverse index  if paired
    if paired and rev_path:
        try:
            # Check if reverse FASTQ file exists
            if rev_path and os.path.exists(rev_path):
                # Define path for reverse index DB
                rev_db = os.path.join(
                    tmpdir, os.path.basename(rev_path) + '.rev.idx'
                )

                # Build the index DB for the reverse reads if it doesn't exist
                if not os.path.exists(rev_db):
                    logging.debug(
                        'Building reverse FASTQ index at %s', rev_db
                    )
                    SeqIO.index_db(rev_db, rev_path, 'fastq')
        except (OSError, ValueError, RuntimeError) as exc:
            logging.warning(
                'Failed to build reverse FASTQ index (%s): %s', rev_path, exc
            )


def index_databases(
    *,  # Enforce keyword arguments
    sample_database: str
) -> str:
    """
    Index the database file with pysam and kma

    Args:
        sample_database: Path to FASTA database file.

    Returns:
        Path to KMA-indexed database prefix.
    """
    # Don't bother re-indexing, this only needs to happen once.
    if not os.path.isfile(sample_database + '.fai'):
        try:
            pysam.faidx(sample_database)
        except pysam.utils.SamtoolsError:
            pass

    # Set the KMA database name
    kma_database = sample_database.replace('.fasta', '') + '_kma'

    # The .name is one of the files KMA creates when making a database.
    if not os.path.isfile(kma_database + '.name'):
        logging.info(
            'Since this is the first time you are using this database, it '
            'needs to be indexed by KMA. This might take a while'
        )

        # Set the KMA indexing command
        cmd = f'kma index -i {sample_database} -o {kma_database}'

        # Set output variable
        out = str()

        # Run the KMA indexing command
        try:
            out, err = run_cmd(cmd=cmd)
        except subprocess.CalledProcessError as exc:
            err = exc

        # Write to logfile
        log = sample_database + '_log.txt'
        write_to_logfile(
            logfile=log,
            out=out,
            err=err,
            cmd=cmd
        )

    return kma_database


def find_contamination(
    *,  # Enforce keyword arguments
    pair: List[str],
    output_folder: str,
    databases_folder: str,
    base_cutoff: int,
    xmx: str,
    forward_id: str = '_R1',
    threads: int = 1,
    keep_files: bool = False,
    quality_cutoff: int = 20,
    min_quality: int = 15,
    base_fraction_cutoff: float = 0.05,
    cgmlst_db: Optional[str] = None,
    tmpdir: Optional[str] = None,
    data_type: str = 'Illumina',
    use_rmlst: bool = False,
    min_matching_hashes: int = 40,
    fasta: bool = False,
    error_cutoff: float = 1.0,
    use_prob_scoring: bool = False,
    score_threshold: float = 2.0,
    max_expected_positions: float = 0.001,
    downsample_depth: Optional[int] = None,
    subreplicates: int = 1,
    subreplicate_seed: Optional[int] = None,
    subreplicate_consensus: float = 0.5
) -> Optional[Tuple[str, bool]]:
    """
    Run contamination detection for a sample (paired or single reads).

    Args:
        pair: Pair or single read list for the sample (e.g. [R1, R2] or
        [single]).
        output_folder: Output folder for results.
        databases_folder: Path to ConFindr databases.
        base_cutoff: Number of bases required to support a variant call.
        forward_id: Identifier for forward reads.
        threads: Number of threads to use.
        keep_files: Keep temporary files if True.
        quality_cutoff: Base quality cutoff (default 20).
        base_fraction_cutoff: Fractional cutoff for variant support.
        cgmlst_db: Optional cgMLST DB name.
        xmx: Optional memory string for BBMap tools.
        tmpdir: Temporary directory to use.
        data_type: 'Illumina' or 'Nanopore'.
        use_rmlst: Prefer rMLST DB when True.
        min_matching_hashes: Mash matching threshold.
        fasta: If True, operate in FASTA-only mode.
        error_cutoff: Error cutoff percentage (default 1.0).
        debug: Enable debug logging.
        use_prob_scoring: Use probabilistic scoring by requiring supported
            positions with statistical evidence.
        score_threshold: Minimum number of statistically supported positions
            required for a sample to be called contaminated.
        subreplicates: Number of independent downsample replicates to run. If
            >1, reads will be downsampled `subreplicates` times and the
            results aggregated.
        subreplicate_seed: Optional integer seed to make downsampling
            reproducible. When set, each replicate will use seed + replicate
            index to vary sampling deterministically.
        subreplicate_consensus: Fraction (0-1) of replicates that must agree
            on a multibase position for it to be reported in the final
            aggregated output.

    Returns:
        Optional tuple (sample_report_tsv, contamination_boolean), or None on
        failure.
    """
    # Set the name and path for the database download date file
    download_date_file = os.path.join(databases_folder, 'download_date.txt')
    if os.path.isfile(download_date_file):
        with open(download_date_file, 'r', encoding='utf-8') as f:
            database_download_date = f.readline().rstrip()
    else:
        database_download_date = 'ND'

    # Define the log file for this sample
    log = os.path.join(output_folder, 'confindr_log.txt')
    if len(pair) == 2:
        sample_name = os.path.split(pair[0])[-1].split(forward_id)[0]
        paired = True
        logging.debug('Sample is paired. Sample name is %s', sample_name)
    else:
        sample_name = os.path.split(pair[0])[-1].split('.')[0]
        paired = False
        logging.debug('Sample is unpaired. Sample name is %s', sample_name)

    # Create directory for this sample
    sample_tmp_dir = os.path.join(output_folder, sample_name)
    os.makedirs(sample_tmp_dir, exist_ok=True)

    logging.info('Checking for cross-species contamination...')

    # Check if samples are paired or single-end and call appropriate function
    if paired:
        genus = find_cross_contamination(
            databases=databases_folder,
            reads=pair,
            sample_name=sample_name,
            tmpdir=sample_tmp_dir,
            log=log,
            threads=threads,
            min_matching_hashes=min_matching_hashes
        )
    else:
        genus = find_cross_contamination(
            databases=databases_folder,
            reads=pair[0],
            sample_name=sample_name,
            tmpdir=sample_tmp_dir,
            log=log,
            threads=threads,
            min_matching_hashes=min_matching_hashes
        )

    # Setup genus-specific databases, if necessary.
    if cgmlst_db is not None:
        # Sanity check that the DB specified is actually a file, otherwise,
        # quit with appropriate error message.
        if not os.path.isfile(cgmlst_db):
            logging.error(
                'ERROR: Specified cgMLST file (%s) does not exist. Please '
                'check the path and try again.', cgmlst_db
            )
            sys.exit(1)
        sample_database = cgmlst_db
    else:
        db_folder = databases_folder
        if genus != 'ND':
            # Logic here is as follows: users can either have both rMLST
            # databases, which cover all of bacteria, cgmlst-derived databases
            # which cover only Escherichia, Salmonella, and Listeria
            # (may add more at some point), or they can have both. They can
            # also set priority to either always use rMLST, or to use my
            # core-genome derived stuff and fall back on rMLST if they're
            # trying to look at a genus I haven't created a scheme for.
            if len(genus.split(':')) > 1:
                predominant_genus = genus.split(':')[0]
            else:
                predominant_genus = genus
            # In the event rmlst databases have priority, always use them.
            if use_rmlst is True:
                sample_database = os.path.join(
                    db_folder,
                    f'{predominant_genus}_db.fasta'
                )

                # Create genus specific database if it doesn't already exist
                # and we have the necessary rMLST files.
                if not os.path.isfile(sample_database):
                    if os.path.isfile(
                        os.path.join(
                            db_folder,
                            'gene_allele.txt'
                            )
                    ) and os.path.isfile(
                        os.path.join(
                            db_folder,
                            'rMLST_combined.fasta'
                        )
                    ):
                        logging.info(
                            'Setting up rMLST genus-specific database for '
                            'genus %s...', predominant_genus
                        )

                        # Calculate the list of alleles
                        allele_list = find_genus_specific_allele_list(
                            profiles_file=os.path.join(
                                db_folder,
                                'gene_allele.txt'
                            ),
                            target_genus=predominant_genus
                        )

                        # Create the database in a temporary directory if
                        # specified
                        if tmpdir:
                            logging.info(
                                'Using temporary directory for database '
                                'creation: %s', tmpdir
                            )

                            # Create the temporary directory if necessary
                            os.makedirs(tmpdir, exist_ok=True)

                            # Set the sample database path
                            sample_database = os.path.join(
                                tmpdir,
                                f'{predominant_genus}_db.fasta'
                            )

                            # Create the allele-specific database
                            setup_allelespecific_database(
                                fasta_file=sample_database,
                                database_folder=db_folder,
                                allele_list=allele_list
                            )
                        else:
                            setup_allelespecific_database(
                                fasta_file=sample_database,
                                database_folder=db_folder,
                                allele_list=allele_list
                            )
            else:
                # Check if a cgderived database is available. If not, try to
                # use rMLST database.
                sample_database = os.path.join(
                    db_folder,
                    f'{predominant_genus}_db_cgderived.fasta'
                )

                # If cgderived database doesn't exist, fall back on rMLST db.
                if not os.path.isfile(sample_database):
                    sample_database = os.path.join(
                        db_folder,
                        f'{predominant_genus}_db.fasta'
                    )

                    # Create genus specific database if it doesn't already
                    # exist and we have the necessary rMLST files.
                    if os.path.isfile(
                        os.path.join(
                            db_folder,
                            'rMLST_combined.fasta'
                        )
                        and os.path.isfile(
                            os.path.join(
                                db_folder,
                                'gene_allele.txt'
                            )
                        ) and not os.path.isfile(
                            sample_database
                        )
                    ):
                        logging.info(
                            'Setting up core genome genus-specific database '
                            'for genus %s...', predominant_genus
                        )

                        # Calculate the list of alleles
                        allele_list = find_genus_specific_allele_list(
                            profiles_file=os.path.join(
                                db_folder,
                                'gene_allele.txt'
                            ),
                            target_genus=predominant_genus
                        )

                        # Create the database in a temporary directory if
                        # specified
                        if tmpdir:
                            logging.info(
                                'Using temporary directory for database '
                                'creation: %s', tmpdir
                            )
                            # Create the temporary directory if necessary
                            os.makedirs(tmpdir, exist_ok=True)

                            # Set the sample database path
                            sample_database = os.path.join(
                                tmpdir,
                                f'{predominant_genus}_db.fasta'
                            )

                            # Create the allele-specific database
                            setup_allelespecific_database(
                                fasta_file=sample_database,
                                database_folder=db_folder,
                                allele_list=allele_list
                            )
                        else:
                            setup_allelespecific_database(
                                fasta_file=sample_database,
                                database_folder=db_folder,
                                allele_list=allele_list
                            )

        else:
            sample_database = os.path.join(db_folder, 'rMLST_combined.fasta')

    # If a user has gotten to this point and they don't have any database
    # available to do analysis because they don't have rMLST downloaded and
    # we don't have a cg-derived database available, boot them with a helpful
    # message.
    if not os.path.isfile(sample_database):
        write_output(
            output_report=os.path.join(output_folder, 'confindr_report.tsv'),
            sample_name=sample_name,
            multi_positions=0,
            genus=genus,
            total_gene_length=0,
            database_download_date=database_download_date
        )
        logging.info(
            'Did not find databases for genus %s. You can download the rMLST '
            'database to get access to all genera (see https://'
            'olc-bioinformatics.github.io/ConFindr/install/). Alternatively, '
            'if you have a high-quality core-genome derived database for your '
            'genome of interest, we would be happy to add it - open an issue '
            'at https://github.com/OLC-Bioinformatics/ConFindr/issues with '
            'the title "Add genus-specific database: %s"\n', genus, genus
            )
        if keep_files is False:
            shutil.rmtree(sample_tmp_dir)
        return

    # Build a list of (pair, tmp_dir, sample_name) to process. This will
    # contain either a single tuple (the original pair) or multiple tuples
    # corresponding to independent downsampled replicates.
    pairs_to_process: List[Tuple[List[str], str, str]] = [
        (pair, sample_tmp_dir, sample_name)
    ]

    if downsample_depth is not None:
        # Only attempt if genus is a single genus and not 'ND'
        if genus == 'ND' or ':' in genus:
            logging.info(
                'Skipping downsampling: genus ambiguous or undetermined (%s).',
                genus
            )
        else:
            # Use first genus token if mash reported multiple with separators
            predominant_genus = genus.split(':')[0]
            genome_size = estimate_genome_size(genus=predominant_genus)
            # Estimate mean read length from the forward read file
            mean_len = estimate_mean_read_length(fastq_path=pair[0])
            if mean_len <= 0:
                logging.warning(
                    'Could not estimate mean read length, skipping '
                    'downsampling.'
                )
            else:
                target_reads = int(
                    (downsample_depth * genome_size) / float(mean_len)
                )
                if target_reads <= 0:
                    logging.warning(
                        'Computed non-positive target reads (%s), skipping '
                        'downsampling.',
                        target_reads
                    )
                else:
                    # Only downsample if actual reads exceed target
                    actual_reads = count_fastq_reads(fastq_path=pair[0])

                    logging.debug(
                        'Downsampling check: actual reads=%s.',
                        actual_reads
                    )
                    if actual_reads and actual_reads > target_reads:
                        logging.info(
                            'Downsampling reads to approximate depth %s '
                            '(target reads: %s)...',
                            downsample_depth, target_reads
                        )

                        # If subreplicates are requested, create multiple
                        # independent downsampled pairs. Otherwise, just
                        # create a single downsampled pair in the sample tmpdir
                        if subreplicates is None or subreplicates <= 1:
                            try:
                                ds_pair = downsample_reads(
                                    pair=pair,
                                    sample_tmp_dir=sample_tmp_dir,
                                    sample_name=sample_name,
                                    target_reads=target_reads,
                                    log=log,
                                    xmx=xmx
                                )
                                pairs_to_process = [
                                    (ds_pair, sample_tmp_dir, sample_name)
                                ]
                            except subprocess.CalledProcessError as exc:
                                logging.warning(
                                    'Downsampling failed: %s. Continuing with '
                                    'original reads.',
                                    str(exc)
                                )
                        else:
                            # Create multiple replicate tmpdirs and run
                            # downsampling for each replicate. If a seed is
                            # provided, use deterministic sampling per
                            # replicate by adding the replicate index.
                            pairs_to_process = []
                            for i in range(int(subreplicates)):
                                rep_idx = i + 1
                                rep_name = f'{sample_name}_rep{rep_idx}'
                                rep_tmp = os.path.join(
                                    sample_tmp_dir, rep_name
                                )
                                os.makedirs(rep_tmp, exist_ok=True)
                                seed = (
                                    subreplicate_seed + i
                                ) if subreplicate_seed is not None else None
                                try:
                                    ds_pair = downsample_reads(
                                        pair=pair,
                                        sample_tmp_dir=rep_tmp,
                                        sample_name=rep_name,
                                        target_reads=target_reads,
                                        log=os.path.join(
                                            rep_tmp, 'downsample.log'
                                        ),
                                        xmx=xmx,
                                        seed=seed
                                    )
                                except subprocess.CalledProcessError as exc:
                                    logging.warning(
                                        'Downsampling failed for replicate %s:'
                                        ' %s. Falling back to original reads.',
                                        rep_name, str(exc)
                                    )
                                    ds_pair = pair
                                pairs_to_process.append(
                                    (ds_pair, rep_tmp, rep_name)
                                )
                    else:
                        logging.debug(
                            'No downsampling required: actual reads (%s) <= '
                            'target (%s).',
                            actual_reads, target_reads
                        )

    # Extract rMLST reads and quality trim.
    logging.info('Extracting conserved core genes...')
    # If multiple downsample replicates were requested, run each replicate
    # through a simplified version of the pipeline (bait -> trim -> map ->
    # read_contig) and aggregate results. For efficiency and isolation each
    # replicate writes into its own tmpdir. After aggregation we write the
    # final contamination TSV and gene summary and return early.
    if len(pairs_to_process) > 1:
        per_run_results = []
        for pair_run, run_tmpdir, run_name in pairs_to_process:
            os.makedirs(run_tmpdir, exist_ok=True)
            local_log = os.path.join(run_tmpdir, os.path.basename(log))

            logging.debug(
                'Baiting paired reads %s with %s',
                pair_run, sample_database
            )

            # Bait
            if paired:
                forward_bait = os.path.join(
                    run_tmpdir,
                    f'{run_name}_baited_R1.fastq.gz'
                )

                reverse_bait = forward_bait.replace('_R1', '_R2')
                out, err, cmd = bbtools.bbduk_bait(
                    reference=sample_database,
                    forward_in=pair_run[0],
                    reverse_in=pair_run[1],
                    forward_out=forward_bait,
                    reverse_out=reverse_bait,
                    threads=threads,
                    Xmx=xmx,
                    returncmd=True
                )
            else:
                if data_type == 'Nanopore' or fasta:
                    unpaired_bait = os.path.join(
                        run_tmpdir,
                        f'{run_name}_baited_trimmed.fastq.gz'
                    )
                else:
                    unpaired_bait = os.path.join(
                        run_tmpdir,
                        f'{run_name}_baited.fastq.gz'
                    )
                out, err, cmd = bbtools.bbduk_bait(
                    reference=sample_database,
                    forward_in=pair_run[0],
                    forward_out=unpaired_bait,
                    returncmd=True,
                    threads=threads,
                    Xmx=xmx
                )
            write_to_logfile(logfile=local_log, out=out, err=err, cmd=cmd)

            # Trim

            logging.debug(
                'Quality trimming baited reads for sample %s...', run_name
            )

            if paired:
                forward_trimmed = os.path.join(
                    run_tmpdir,
                    f'{run_name}_trimmed_R1.fastq.gz'
                )
                reverse_trimmed = os.path.join(
                    run_tmpdir,
                    f'{run_name}_trimmed_R2.fastq.gz'
                )
                out, err, cmd = bbtools.bbduk_trim(
                    forward_in=forward_bait,
                    reverse_in=reverse_bait,
                    forward_out=forward_trimmed,
                    reverse_out=reverse_trimmed,
                    qtrim=quality_cutoff,
                    threads=threads,
                    Xmx=xmx,
                    returncmd=True
                )
            else:
                unpaired_trimmed = os.path.join(
                    run_tmpdir,
                    f'{run_name}_trimmed.fastq.gz'
                )
                out, err, cmd = bbtools.bbduk_trim(
                    forward_in=unpaired_bait,
                    forward_out=unpaired_trimmed,
                    qtrim=quality_cutoff,
                    threads=threads,
                    Xmx=xmx,
                    returncmd=True
                )
            write_to_logfile(logfile=local_log, out=out, err=err, cmd=cmd)

            # Index database
            kma_database = index_databases(sample_database=sample_database)

            # KMA mapping and allele selection
            run_kma_report = os.path.join(run_tmpdir, f'{run_name}_kma')
            if paired:
                cmd = (
                    f'kma -ipe {forward_trimmed} {reverse_trimmed} '
                    f'-o {run_kma_report} -t_db {kma_database} -t {threads}'
                )
            else:
                cmd = (
                    f'kma -i {unpaired_trimmed} -o {run_kma_report} '
                    f'-t_db {kma_database} -t {threads}'
                )

            logging.debug(
                'Running KMA with command: %s', cmd
            )

            # Run KMA command
            out, err = run_cmd(cmd=cmd)

            # Write to logfile
            write_to_logfile(logfile=local_log, out=out, err=err, cmd=cmd)

            # Determine alleles present
            gene_alleles_run = find_rmlst_type(
                kma_report=run_kma_report + '.res',
                rmlst_report=os.path.join(
                    run_tmpdir,
                    f'{run_name}_alleles.tsv')
            )

            logging.debug(
                'Gene alleles for sample %s: %s',
                run_name,
                gene_alleles_run
            )

            # Diagnostic: check that KMA-reported alleles actually exist
            # in the source database. If some are missing, log a warning and
            # filter them out so they don't cause downstream errors (faidx
            # failures while parsing BAM/pileups).
            try:
                # Load database contig IDs
                db_contigs = set(
                    c.id for c in SeqIO.parse(sample_database, 'fasta')
                )

                # Find any missing alleles
                missing_in_db = [
                    g for g in gene_alleles_run if g not in db_contigs
                ]

                # Log any missing alleles and filter them out
                if missing_in_db:
                    logging.warning(
                        'KMA reported %s alleles that are not present in %s '
                        '(first 20): %s',
                        len(missing_in_db),
                        sample_database,
                        missing_in_db[:20]
                    )

                    # Filter out missing alleles to avoid downstream errors
                    gene_alleles_run = [
                        g for g in gene_alleles_run if g in db_contigs
                    ]

                    logging.info(
                        'Proceeding with %s alleles present in the database '
                        '(filtered out %s).',
                        len(gene_alleles_run), len(missing_in_db)
                    )
            except (TypeError, ValueError, AttributeError) as exc:
                logging.debug(
                    'Could not validate KMA allele presence in DB: %s', exc
                )

            # If no alleles remain after filtering, abort this replicate/sample
            if not gene_alleles_run:
                logging.warning(
                    'No valid alleles found for run %s after filtering; '
                    'skipping this replicate.',
                    run_name
                )

                # Record the run name and empty results so we can track
                # per-run counts
                per_run_results.append((run_name, [], [], 0))
                continue

            # Set the path to the allele-specific FASTA
            run_rmlst_fasta = os.path.join(
                run_tmpdir,
                f'{run_name}_alleles.fasta'
            )

            # Create allele-specific FASTA if it doesn't already exist
            if not os.path.isfile(run_rmlst_fasta):
                with open(run_rmlst_fasta, 'w', encoding='utf-8') as f:
                    for contig in SeqIO.parse(sample_database, 'fasta'):
                        if contig.id in gene_alleles_run:
                            SeqIO.write(contig, f, 'fasta')

                # Ensure the FASTA is faidx-indexed for pysam operations
                try:
                    if not os.path.isfile(run_rmlst_fasta + '.fai'):
                        pysam.faidx(run_rmlst_fasta)
                        logging.debug(
                            'Created FASTA index for %s', run_rmlst_fasta
                        )
                except (OSError, pysam.utils.SamtoolsError) as exc:
                    logging.debug(
                        'Could not create FASTA index for %s: %s',
                        run_rmlst_fasta, exc
                    )

            # Calculate total length of alleles in this run
            rmlst_gene_length_run = find_total_sequence_length(
                fasta_file=run_rmlst_fasta
            )

            logging.debug(
                'Total rMLST gene length for sample %s: %s',
                run_name,
                rmlst_gene_length_run
            )

            # Map back and parse with read_contig
            try:
                # Set the path to the sorted BAM file
                sorted_bam_run = os.path.join(
                    run_tmpdir,
                    f'{run_name}_contamination_sorted.bam'
                )

                # Set the path to the output BAM file
                outbam_run = os.path.join(
                    run_tmpdir,
                    f'{run_name}_contamination.bam'
                )

                # Perform mapping
                if paired:
                    cmd = (
                        f'bbmap.sh ref={run_rmlst_fasta} in={forward_trimmed} '
                        f'in2={reverse_trimmed} out={outbam_run} '
                        f'threads={threads} -Xmx{xmx}  mdtag nodisk'
                    )
                else:
                    if data_type == 'Illumina' and not fasta:
                        cmd = (
                            f'bbmap.sh ref={run_rmlst_fasta} '
                            f'in={unpaired_trimmed} out={outbam_run} '
                            f'threads={threads} -Xmx{xmx}  mdtag nodisk'
                        )
                    else:
                        ax = 'asm5' if fasta else 'map-ont'
                        cmd = (
                            f'minimap2 --MD -t {threads} -ax {ax} '
                            f'{run_rmlst_fasta} {unpaired_bait} '
                            f' | samtools view -@ {threads} -h '
                            f'-bT {run_rmlst_fasta} - | '
                            f'samtools sort - -@ {threads} -o {sorted_bam_run}'
                        )

                logging.debug(
                    'Mapping reads for sample %s with command: %s',
                    run_name,
                    cmd
                )

                # Run mapping command
                out, err = run_cmd(cmd=cmd)

                # Write to logfile
                write_to_logfile(
                    logfile=local_log,
                    out=out,
                    err=err,
                    cmd=cmd
                )

                # Sort and index BAM if necessary
                if (
                    not os.path.isfile(sorted_bam_run)
                    and os.path.isfile(outbam_run)
                ):
                    pysam.sort('-o', sorted_bam_run, outbam_run)

                # Index BAM if necessary
                if not os.path.isfile(sorted_bam_run + '.bai'):
                    pysam.index(sorted_bam_run)

                # Load allele records
                allele_records_run = SeqIO.to_dict(
                    SeqIO.parse(
                        run_rmlst_fasta,
                        'fasta'
                    )
                )

                logging.debug(
                    'Loaded %s allele records for sample %s',
                    len(allele_records_run),
                    run_name
                )

                # Load FASTQ records into memory (reverted behaviour) so worker
                # tasks can access read qualities quickly and avoid per-worker
                # index contention. This follows the old, fast behaviour (may
                # increase memory usage on low-RAM systems).
                if paired:
                    with gzip.open(pair_run[0], 'rt') as gz:
                        fastq_records_run = load_fastq_records(
                            gz=gz,
                            paired=True,
                            forward=True
                        )
                    with gzip.open(pair_run[1], 'rt') as gz:
                        fastq_records_run.update(
                            load_fastq_records(
                                gz=gz,
                                paired=True,
                                forward=False
                            )
                        )
                else:
                    with gzip.open(pair_run[0], 'rt') as gz:
                        fastq_records_run = load_fastq_records(
                            gz=gz,
                            paired=False,
                            forward=True
                        )

                logging.debug(
                    'Loaded %s FASTQ records for sample %s (in-memory)',
                    len(fastq_records_run),
                    run_name
                )

                # Prepare multiprocessing arguments. First ensure the alleles
                # we're about to analyse are present in the BAM header to avoid
                # pysam 'reference sequence not found' errors.
                try:
                    bam_for_check = pysam.AlignmentFile(sorted_bam_run, 'rb')

                    # Normalize to strings in case pysam returns bytes
                    bam_refs = set(
                        x.decode() if isinstance(
                            x,
                            bytes
                        ) else str(x) for x in bam_for_check.references
                    )

                    # Close the BAM file after reading references
                    bam_for_check.close()
                except (
                    OSError,
                    ValueError,
                    SamtoolsError,
                    AttributeError,
                    UnicodeDecodeError
                ) as exc:
                    logging.debug(
                        'Could not read BAM references for %s: %s',
                        sorted_bam_run, exc
                    )
                    bam_refs = set()

                # Ensure gene IDs are strings and filter by BAM header refs
                gene_alleles_run = [
                    g.decode() if isinstance(
                        g, bytes
                    ) else str(g) for g in gene_alleles_run
                ]
                present_alleles = [
                    g for g in gene_alleles_run if g in bam_refs
                ]
                missing_alleles = [
                    g for g in gene_alleles_run if g not in bam_refs
                ]

                # Log any missing alleles
                if missing_alleles:
                    logging.warning(
                        'The following alleles reported by KMA are not '
                        'present in BAM header and will be skipped for run '
                        '%s (first 20): %s',
                        run_name, missing_alleles[:20]
                    )
                if not present_alleles:
                    logging.warning(
                        'No alleles remain after filtering by BAM header '
                        'for run %s; skipping.',
                        run_name
                    )

                    # Preserve run name even when no alleles remained so we
                    # keep column alignment
                    per_run_results.append((run_name, [], [], 0))
                    continue

                logging.debug(
                    'Preparing to process %s alleles for sample %s',
                    len(gene_alleles_run),
                    run_name
                )

                # Build contig-length map using FASTA index where possible
                contig_lengths_run = {}
                try:
                    # Set the path to the FASTA index
                    fai_path = run_rmlst_fasta + '.fai'

                    # Load lengths from .fai if present, otherwise parse FASTA
                    if os.path.isfile(fai_path):
                        with open(fai_path, 'r', encoding='utf-8') as fh:
                            # .fai format: contig_name \t length \t ...
                            for line in fh:
                                # Parse first two columns
                                cols = line.strip().split('\t')

                                # Store contig length
                                if len(cols) >= 2:
                                    contig_lengths_run[cols[0]] = int(cols[1])
                    # Fallback: parse FASTA directly
                    else:
                        for rec in SeqIO.parse(run_rmlst_fasta, 'fasta'):
                            contig_lengths_run[rec.id] = len(rec.seq)
                except (OSError, ValueError, AttributeError) as exc:
                    logging.debug(
                        'Could not obtain contig lengths from '
                        'FASTA/.fai: %s',
                        exc
                    )

                # Partition contigs into balanced chunks
                chunks = _build_contig_chunks(
                    contig_lengths=contig_lengths_run,
                    contigs=present_alleles,
                    threads=threads,
                    multiplier=CONTIG_CHUNK_MULTIPLIER,
                    max_chunk_bases=CONTIG_CHUNK_MAX_BASES
                )
                logging.info(
                    'Partitioned %s alleles into %s chunks for run %s '
                    '(multiplier=%s, max_chunk_bases=%s)',
                    len(present_alleles), len(chunks), run_name,
                    CONTIG_CHUNK_MULTIPLIER, CONTIG_CHUNK_MAX_BASES
                )

                # Build kwargs list (one chunk per task)
                kwargs_list = []
                for chunk in chunks:
                    kwargs_list.append({
                        'contig_chunk': chunk,
                        'bamfile_name': sorted_bam_run,
                        'reference_fasta': run_rmlst_fasta,
                        'fastq_records': fastq_records_run,
                        'quality_cutoff': quality_cutoff,
                        'base_cutoff': base_cutoff,
                        'base_fraction_cutoff': base_fraction_cutoff,
                        'fasta': fasta,
                        'error_cutoff': error_cutoff,
                        'nanopore': True if data_type == 'Nanopore' else False,
                        'max_expected_positions': max_expected_positions
                    })

                # If we have the FASTQ records loaded in memory for this run
                # (debug single-thread mode or user-forced), prefer ThreadPool
                # so we can share the in-memory dict without pickling.
                use_threadpool = fastq_records_run is not None

                if use_threadpool:
                    logging.info(
                        'Using in-memory ThreadPool for run %s (fastq '
                        'records loaded in master).',
                        run_name
                    )
                    p = ThreadPool(processes=threads)
                else:
                    # Ensure index DBs exist (build once in master to avoid
                    # worker contention), then spawn pool with per-worker index
                    _ensure_fastq_index(
                        forward_trimmed if paired else unpaired_trimmed,
                        reverse_trimmed if paired else None,
                        paired,
                        run_tmpdir
                    )
                    p = multiprocessing.Pool(
                        processes=threads,
                        initializer=_init_fastq_index,
                        initargs=(
                            forward_trimmed if paired else unpaired_trimmed,
                            reverse_trimmed if paired else None, paired,
                            run_tmpdir
                        )
                    )

                # Initialize results lists
                multibase_dicts_run = []
                report_lines_run = []

                # Process the results with robust error handling
                try:
                    # Set the total contigs as the number of present alleles
                    total_contigs = len(present_alleles)

                    # Initialize processed contig counter and timing
                    processed_contigs = 0
                    tstart = time.time()

                    # Determine logging frequency
                    log_every = max(1, total_contigs // 20)  # log ~20 updates

                    # Process chunks in parallel
                    for multibase_dict, report_write in p.map(
                        _read_contig_chunk_dispatch,
                        kwargs_list,
                        chunksize=1
                    ):
                        # Iterate over returned dict and append per-contig
                        # results
                        for k, v in multibase_dict.items():
                            multibase_dicts_run.append({k: v})

                            # Update processed contig count
                            processed_contigs += 1

                        # Append report lines for this chunk if specified
                        if report_write:
                            # Split lines and filter by contig IDs in this
                            # chunk
                            lines = report_write.splitlines(True)

                            # Append lines for each contig in the chunk
                            for contig in multibase_dict.keys():
                                # Filter lines for this contig
                                contig_lines = [
                                    ln for ln in lines if ln.startswith(
                                        contig + '\t'
                                    )
                                ]

                                # Append contig lines to report
                                report_lines_run.append(contig_lines)

                        # Log progress periodically
                        if (
                            processed_contigs % log_every == 0
                            or processed_contigs == total_contigs
                        ):
                            # Compute elapsed time and ETA
                            elapsed = time.time() - tstart
                            avg = elapsed / max(1, processed_contigs)
                            eta = avg * (total_contigs - processed_contigs)

                            logging.info(
                                'Run %s progress: %s/%s (%.1f%%). '
                                'Elapsed: %s, ETA: %s',
                                run_name,
                                processed_contigs,
                                total_contigs,
                                processed_contigs / total_contigs * 100.0,
                                _format_seconds(s=elapsed),
                                _format_seconds(s=eta)
                            )
                except (
                    SamtoolsError,
                    OSError,
                    ValueError,
                    KeyError,
                    IndexError,
                    RuntimeError
                ) as exc:
                    logging.exception(
                        'Pool processing failed for run %s; terminating '
                        'workers: %s',
                        run_name, exc
                    )
                    try:
                        p.terminate()
                    except (OSError, RuntimeError) as exc_term:
                        logging.debug('Failed to terminate pool: %s', exc_term)
                    try:
                        p.join()
                    except (OSError, RuntimeError) as exc_join:
                        logging.debug('Failed to join pool: %s', exc_join)
                    raise
                finally:
                    try:
                        p.close()
                        p.join()
                    except (OSError, RuntimeError) as exc:
                        logging.debug('Pool close/join failed: %s', exc)
            except SamtoolsError:
                multibase_dicts_run = []
                report_lines_run = []

            # Store per-run results
            per_run_results.append(
                (
                    run_name,
                    multibase_dicts_run,
                    report_lines_run,
                    rmlst_gene_length_run
                )
            )

        # Aggregate per-run results
        pos_counter = {}
        genes_seen = set()

        # Process each run's results and track per-run counts
        for run_idx, run_entry in enumerate(per_run_results):
            # Unpack run tuple; supports new format (run_name, multibase,
            # report_lines, length)
            try:
                run_name, _multibase, report_lines, rmlst_len = run_entry
            except ValueError:
                # Backwards compatibility: fall back to old tuple format
                _, report_lines, rmlst_len = run_entry

            # Process each line of the report
            for item in report_lines:
                # Handle both list of lines and single string
                lines = item if isinstance(
                    item,
                    list
                ) else item.splitlines(True)

                # Process each line
                for contamination_info in lines:
                    # Parse the TSV line
                    fields = contamination_info.rstrip('\n').split('\t')

                    # Sanity check line has enough fields
                    if len(fields) < 2:
                        continue

                    # Extract gene and position
                    gene = fields[0]
                    pos = fields[1]

                    # Create unique key
                    key = (gene, pos)

                    # Update position counter (and ensure per-run dict exists)
                    pos_entry = pos_counter.setdefault(
                        key,
                        {
                            'count': 0,
                            'line': contamination_info,
                            'pvals': [],
                            'qvals': [],
                            'per_run_counts': {},
                            'per_run_metrics': {}
                        }
                    )
                    pos_entry['count'] += 1

                    # Increment run-specific count for this position
                    pos_entry['per_run_counts'][run_idx] = \
                        pos_entry['per_run_counts'].get(run_idx, 0) + 1

                    # Robustly extract numeric fields from the end of the line
                    # (see other aggregation branch)
                    n_fields = len(fields)
                    mean_mapq_field = fields[-1] if n_fields >= 1 else 'ND'
                    mean_q_field = fields[-2] if n_fields >= 2 else 'ND'
                    adjp_field = fields[-5] if n_fields >= 5 else 'ND'
                    p_field = fields[-6] if n_fields >= 6 else 'ND'
                    total_cov_field = fields[-9] if n_fields >= 9 else 'ND'
                    snv_cov_field = fields[-10] if n_fields >= 10 else 'ND'

                    def safe_int(x):
                        try:
                            return int(x)
                        except (ValueError, TypeError):
                            return None

                    def safe_float(x):
                        try:
                            return float(x)
                        except (ValueError, TypeError):
                            return None

                    # Parse numeric fields
                    snv_coverage = safe_int(snv_cov_field)
                    total_coverage = safe_int(total_cov_field)
                    mean_q = safe_float(mean_q_field)
                    mean_mapq = safe_float(mean_mapq_field)

                    # Extract forward/reverse SNV counts from the read-type
                    # block(s)
                    num_start = max(4, n_fields - 10)

                    # Combine read-type fields in case of extra tabs
                    read_types_combined = '\t'.join(
                        fields[3:num_start]
                    ) if num_start > 3 else fields[3] if n_fields > 3 else ''

                    # Split by commas and filter empty parts
                    parts = [
                        p for p in read_types_combined.split(',') if p != ''
                    ]

                    def sum_base_counts(part):
                        if not part:
                            return 0
                        total = 0
                        for token in part.split(';'):
                            token = token.strip()
                            if not token or ':' not in token:
                                continue
                            try:
                                total += int(token.split(':', 1)[1])
                            except (ValueError, TypeError):
                                # Skip tokens with non-integer or unexpected
                                # values
                                continue
                        return total

                    # Sum forward and reverse SNV counts
                    forward_total = sum_base_counts(
                        parts[2]
                    ) if len(parts) > 2 else 0

                    reverse_total = sum_base_counts(
                        parts[3]
                    ) if len(parts) > 3 else 0

                    # Store per-run metrics
                    pos_entry['per_run_metrics'][run_idx] = {
                        'snv_coverage': snv_coverage,
                        'total_coverage': total_coverage,
                        'mean_qual': mean_q,
                        'mean_mapq': mean_mapq,
                        'forward_snvs': forward_total,
                        'reverse_snvs': reverse_total
                    }

                    # Extract p-value and q-value. Use reverse indexing as well
                    try:
                        pval = safe_float(p_field) if p_field != 'ND' else None
                    except (ValueError, IndexError):
                        pval = None
                    try:
                        qval = safe_float(
                            adjp_field
                        ) if adjp_field != 'ND' else None
                    except (ValueError, IndexError):
                        qval = None

                    # Append p-value and q-value if available
                    if pval is not None:
                        pos_entry['pvals'].append(pval)
                    if qval is not None:
                        pos_entry['qvals'].append(qval)
                    genes_seen.add(gene)

        # Determine threshold for reporting positions
        threshold = math.ceil(
            subreplicate_consensus * max(1, len(pairs_to_process))
        )

        # Initialize final outputs
        final_report_lines = []
        final_multibase = {}
        final_rmlst_len = 0

        # Compile final results based on threshold
        for (gene, pos), entry in pos_counter.items():
            if entry['count'] >= threshold:
                # Append per-subreplicate detailed metrics to the
                # representative line if subreplicates were run
                if len(per_run_results) > 1:
                    cols = []
                    for i in range(len(per_run_results)):
                        # Count of runs observing this SNV (may be 0)
                        run_count = entry.get('per_run_counts', {}).get(i, 0)
                        metrics = entry.get('per_run_metrics', {}).get(i, {})

                        snv_cov = metrics.get('snv_coverage')
                        total_cov = metrics.get('total_coverage')
                        mean_q = metrics.get('mean_qual')
                        mean_mapq = metrics.get('mean_mapq')
                        fwd = metrics.get('forward_snvs')
                        rev = metrics.get('reverse_snvs')

                        # Format mean_q and mean_mapq
                        mean_q_str = f'{mean_q:0.2f}' if isinstance(
                            mean_q, (int, float)
                        ) else 'ND'
                        mean_mapq_str = f'{mean_mapq:0.2f}' if isinstance(
                            mean_mapq, (int, float)
                        ) else 'ND'

                        # Append metrics for this subreplicate
                        cols.extend([
                            str(run_count),
                            str(snv_cov) if snv_cov is not None else 'ND',
                            str(total_cov) if total_cov is not None else 'ND',
                            mean_q_str,
                            mean_mapq_str,
                            str(fwd) if fwd is not None else 'ND',
                            str(rev) if rev is not None else 'ND'
                        ])

                    # Build augmented line
                    augmented_line = entry['line'].rstrip(
                        '\n'
                    ) + '\t' + '\t'.join(cols) + '\n'
                    final_report_lines.append(augmented_line)
                else:
                    # No subreplicates; append line as-is
                    final_report_lines.append(entry['line'])

                # Initialize position in final multibase dict
                gdict = final_multibase.setdefault(gene, {})
                gdict[pos] = {}

        # Determine final allele length (support both old and new per-run
        # tuple formats)
        for entry in per_run_results:
            # entry may be (run_name, multibase, report_lines, rmlst_len) or
            # (multibase, report_lines, rmlst_len)
            if isinstance(entry, tuple) and len(entry) == 4:
                rmlst_len = entry[3]
            elif isinstance(entry, tuple) and len(entry) == 3:
                rmlst_len = entry[2]
            else:
                continue
            if rmlst_len:
                final_rmlst_len = rmlst_len
                break

        # Compute per-gene statistics
        for gene in genes_seen:
            # Initialize lists for p-values and q-values
            pvals = []
            qvals = []

            # Collect p-values and q-values for this gene
            for key, entry in pos_counter.items():
                if key[0] != gene:
                    continue
                pvals.extend(entry['pvals'])
                qvals.extend(entry['qvals'])

            # Compute combined statistics
            combined_p = combine_pvalues_fisher(pvals=pvals) if pvals else None

            # Compute gene score and number of significant positions
            gene_score = sum(
                [-math.log10(q + 1e-300) for q in qvals]
            ) if qvals else None

            # Count significant positions (q <= 0.05)
            num_sig_positions = sum(
                1 for q in qvals if q is not None and q <= 0.05
            )

            # Store in final multibase dictionary
            if gene not in final_multibase:
                final_multibase[gene] = {}

            final_multibase[gene]['_gene_stats'] = {
                'combined_p': combined_p,
                'gene_score': gene_score,
                'num_sig_positions': num_sig_positions
            }

        # Write aggregated contamination TSV
        report_file = os.path.join(
            output_folder,
            sample_name + '_contamination.tsv'
        )
        with open(report_file, 'w', encoding='utf-8') as r:
            base_header = (
                'Gene\tPosition\tRefBase\tCongruentSNVs\tTotalSNVs\t'
                'ForwardSNVs\tReverseSNVs\tSNVCoverage\tTotalCoverage\t'
                'BaseCutoff\tErrorPercent\tPValue\tAdjPValue\tStrandP\t'
                'PosP\tMeanQual\tMeanMapQ'
            )
            if len(per_run_results) > 1:
                # For each subreplicate add multiple per-run metric columns
                sub_header_parts = []
                for i in range(len(per_run_results)):
                    idx = i + 1
                    sub_header_parts.extend([
                        f'Subreplicate_{idx}_Count',
                        f'Subreplicate_{idx}_SNVcov',
                        f'Subreplicate_{idx}_TotalCov',
                        f'Subreplicate_{idx}_MeanQual',
                        f'Subreplicate_{idx}_MeanMapQ',
                        f'Subreplicate_{idx}_ForwardSNVs',
                        f'Subreplicate_{idx}_ReverseSNVs'
                    ])
                sub_headers = '\t' + '\t'.join(sub_header_parts)
            else:
                sub_headers = ''
            r.write(base_header + sub_headers + '\n')
            for line in final_report_lines:
                r.write(line)

        # Write aggregated gene summary
        gene_summary_file = os.path.join(
            output_folder, sample_name + '_gene_summary.tsv'
        )

        # Write gene summary with computed statistics
        with open(gene_summary_file, 'w', encoding='utf-8') as gf:
            gf.write(
                'Gene\tCombinedP\tGeneScore\tNumSigPositions\tNumPositions\n'
            )

            # Iterate over genes and write stats
            for gene, content in final_multibase.items():
                stats = content.get('_gene_stats', {})
                combined_p = stats.get('combined_p')
                gene_score_val = stats.get('gene_score')
                num_sig_positions = stats.get('num_sig_positions', 0)
                num_positions = sum(
                    1 for k in content.keys() if not str(k).startswith('_')
                )
                combined_p_str = f'{combined_p:0.3e}' if combined_p \
                    is not None else 'ND'
                gene_score_str = f'{gene_score_val:0.3f}' if gene_score_val \
                    is not None else 'ND'
                gf.write(
                    f'{gene}\t{combined_p_str}\t{gene_score_str}\t'
                    f'{num_sig_positions}\t{num_positions}\n'
                )

        # Compute final summary metrics
        multi_positions = sum(
            1 for content in final_multibase.values()
            for k in content.keys() if not str(k).startswith('_')
        )
        rmlst_gene_length = final_rmlst_len
        pysam_pass = True

        # Use the number of supported positions as the probabilistic score.
        sample_score = float(multi_positions)

        # Final write to summary report
        if keep_files is False:
            shutil.rmtree(sample_tmp_dir)

        # Compute per-run SNV counts for summary report (one count
        # per subreplicate)
        per_run_snv_counts = []
        for run_entry in per_run_results:
            try:
                run_name, _multibase, report_lines, rmlst_len = run_entry
            except ValueError:
                # Backwards compatibility: older entry format
                _multibase, report_lines, rmlst_len = run_entry

            # Initialize position set
            positions = set()

            # Process each line of the report
            for item in report_lines:
                # Handle both list of lines and single string
                lines = item if isinstance(
                    item,
                    list
                ) else item.splitlines(True)

                # Process each line
                for contamination_info in lines:
                    # Parse the TSV line
                    fields = contamination_info.rstrip('\n').split('\t')

                    # Sanity check line has enough fields
                    if len(fields) < 2:
                        continue

                    # Extract gene and position and add to set
                    positions.add((fields[0], fields[1]))

            # Append count of unique positions for this run
            per_run_snv_counts.append(len(positions))

        write_output(
            output_report=os.path.join(output_folder, 'confindr_report.tsv'),
            sample_name=sample_name,
            multi_positions=multi_positions,
            genus=genus,
            total_gene_length=rmlst_gene_length,
            snp_cutoff=math.ceil(
                rmlst_gene_length / 10000
            ) + 1 if cgmlst_db is None else (1 if fasta else 10),
            database_download_date=database_download_date,
            pysam_pass=pysam_pass,
            sample_score=sample_score,
            use_probabilistic=use_prob_scoring,
            score_threshold=score_threshold,
            subreplicate_counts=per_run_snv_counts
        )

        logging.info(
            'Done (replicate aggregated run). Number of contaminating '
            'SNVs found: %s\n',
            multi_positions
        )

        return

    # Non-replicate processing - initialize variables
    out, err, cmd = '', '', ''
    forward_bait = str()
    forward_trimmed = str()
    reverse_bait = str()
    reverse_trimmed = str()
    unpaired_bait = str()
    unpaired_trimmed = str()

    # Bait reads matching the database
    if paired:
        forward_bait = os.path.join(
            sample_tmp_dir,
            f'{sample_name}_baited_R1.fastq.gz'
        )
        reverse_bait = forward_bait.replace('_R1', '_R2')

        # Only run if the baited files don't already exist
        if not os.path.isfile(forward_bait):
            if xmx is None:
                out, err, cmd = bbtools.bbduk_bait(
                    reference=sample_database,
                    forward_in=pair[0],
                    reverse_in=pair[1],
                    forward_out=forward_bait,
                    reverse_out=reverse_bait,
                    threads=threads,
                    Xmx=xmx,
                    returncmd=True
                )
            else:
                out, err, cmd = bbtools.bbduk_bait(
                    reference=sample_database,
                    forward_in=pair[0],
                    reverse_in=pair[1],
                    forward_out=forward_bait,
                    reverse_out=reverse_bait,
                    threads=threads,
                    Xmx=xmx,
                    returncmd=True
                )
    else:
        # Still name the file '_baited_trimmed' even if the file won't be
        # trimmed
        if data_type == 'Nanopore' or fasta:
            unpaired_bait = os.path.join(
                sample_tmp_dir,
                f'{sample_name}_baited_trimmed.fastq.gz'
            )
        else:
            unpaired_bait = os.path.join(
                sample_tmp_dir,
                f'{sample_name}_baited.fastq.gz'
            )

        # Only run if the baited file doesn't already exist
        if not os.path.isfile(unpaired_bait):
            if xmx is None:
                out, err, cmd = bbtools.bbduk_bait(
                    reference=sample_database,
                    forward_in=pair[0],
                    forward_out=unpaired_bait,
                    Xmx=xmx,
                    returncmd=True, threads=threads
                )
            else:
                out, err, cmd = bbtools.bbduk_bait(
                    reference=sample_database,
                    forward_in=pair[0],
                    forward_out=unpaired_bait,
                    Xmx=xmx,
                    returncmd=True,
                    threads=threads
                )
    if out:
        write_to_logfile(
            logfile=log,
            out=out,
            err=err,
            cmd=cmd
        )

    # Verify baiting output exists before trimming. If the expected bait
    # files are missing, write a failure line to the summary report and
    # abort this sample early.
    if paired:
        if not (os.path.isfile(forward_bait) and os.path.isfile(reverse_bait)):
            logging.error(
                'Baiting failed to create files for %s: %s, %s',
                sample_name, forward_bait, reverse_bait
            )
            write_output(
                output_report=os.path.join(
                    output_folder, 'confindr_report.tsv'
                ),
                sample_name=sample_name,
                multi_positions=0,
                genus='Error processing sample',
                total_gene_length=0,
                database_download_date=database_download_date
            )
            if keep_files is False:
                try:
                    shutil.rmtree(sample_tmp_dir)
                except OSError:
                    pass
            return
    else:
        if not os.path.isfile(unpaired_bait):
            logging.error(
                'Baiting failed to create file for %s: %s',
                sample_name, unpaired_bait
            )
            write_output(
                output_report=os.path.join(
                    output_folder, 'confindr_report.tsv'
                ),
                sample_name=sample_name,
                multi_positions=0,
                genus='Error processing sample',
                total_gene_length=0,
                database_download_date=database_download_date
            )
            if keep_files is False:
                try:
                    shutil.rmtree(sample_tmp_dir)
                except OSError:
                    pass
            return

    # Run quality trimming on the baited reads
    logging.info('Quality trimming...')
    out, err, cmd = '', '', ''

    # Handle Illumina and Nanopore data differently
    if data_type == 'Illumina':
        if paired:
            forward_trimmed = os.path.join(
                sample_tmp_dir,
                f'{sample_name}_baited_trimmed_R1.fastq.gz'
            )

            # Set the reverse trimmed filename
            reverse_trimmed = forward_trimmed.replace('_R1', '_R2')

            # Only run trimming if the trimmed files don't already exist
            if not os.path.isfile(forward_trimmed):
                # Use the appropriate bbduk trimming command based on whether
                # xmx is set
                if xmx is None:
                    out, err, cmd = bbtools.bbduk_trim(
                        forward_in=forward_bait,
                        reverse_in=reverse_bait,
                        forward_out=forward_trimmed,
                        reverse_out=reverse_trimmed,
                        threads=str(threads),
                        Xmx=xmx,
                        returncmd=True
                    )
                else:
                    out, err, cmd = bbtools.bbduk_trim(
                        forward_in=forward_bait,
                        reverse_in=reverse_bait,
                        forward_out=forward_trimmed,
                        reverse_out=reverse_trimmed,
                        Xmx=xmx,
                        threads=str(threads),
                        returncmd=True
                    )

            # Load the trimmed FASTQ records into a dictionary only in
            # debug mode. For parallel execution, use per-worker FASTQ
            # indexes to avoid loading all reads into the master process.
            # Revert to previous behaviour: load FASTQ records into memory to
            # reduce per-worker overhead. Note: this may increase memory usage.
            if paired:
                with gzip.open(forward_trimmed, 'rt') as gz:
                    fastq_records = load_fastq_records(
                        gz=gz,
                        paired=True,
                        forward=True
                    )
                with gzip.open(reverse_trimmed, 'rt') as gz:
                    fastq_records.update(
                        load_fastq_records(
                            gz=gz,
                            paired=True,
                            forward=False
                        )
                    )
            else:
                with gzip.open(unpaired_trimmed, 'rt') as gz:
                    fastq_records = load_fastq_records(
                        gz=gz,
                        paired=False,
                        forward=True
                    )
            logging.debug(
                'Loaded %s FASTQ records for sample %s (in-memory)',
                len(fastq_records),
                sample_name
            )
        # Handle unpaired reads
        else:
            unpaired_trimmed = os.path.join(
                sample_tmp_dir,
                f'{sample_name}_baited_trimmed.fastq.gz'
            )

            # Process FASTA mode differently - no trimming
            if not fasta:
                # Only run trimming if the trimmed file doesn't already exist
                if not os.path.isfile(unpaired_trimmed):
                    # Use the appropriate bbduk trimming command based on
                    # whether xmx is set
                    if xmx is None:
                        out, err, cmd = bbtools.bbduk_trim(
                            forward_in=unpaired_bait,
                            forward_out=unpaired_trimmed,
                            returncmd=True,
                            Xmx=xmx,
                            threads=threads
                        )
                    else:
                        out, err, cmd = bbtools.bbduk_trim(
                            forward_in=unpaired_bait,
                            forward_out=unpaired_trimmed,
                            returncmd=True,
                            threads=threads,
                            Xmx=xmx
                        )

                # Load the trimmed FASTQ records into a dictionary only in
                # debug single-threaded mode; otherwise workers will use index
                # DBs on disk.
                if threads == 1:
                    with gzip.open(unpaired_trimmed, 'rt') as gz:
                        # Load the FASTQ records
                        fastq_records = load_fastq_records(
                            gz=gz,
                            paired=False,
                            forward=True
                        )
                else:
                    fastq_records = None
                    logging.debug(
                        'Not loading unpaired FASTQ records into memory; '
                        'workers will use index_db on disk'
                    )
            else:
                # Unpaired_bait
                if threads == 1:
                    with gzip.open(unpaired_bait, 'rt') as gz:
                        # Load the FASTQ records
                        fastq_records = load_fastq_records(
                            gz=gz,
                            paired=False,
                            forward=True
                        )
                else:
                    fastq_records = None
                    logging.debug(
                        'Not loading unpaired FASTQ records into memory; '
                        'workers will use index_db on disk'
                    )

        # Write to logfile
        write_to_logfile(
            logfile=log,
            out=out,
            err=err,
            cmd=cmd
        )

    # If Nanopore data, no trimming - just load the baited reads
    else:
        if paired:
            if threads == 1:
                with gzip.open(forward_bait, 'rt') as gz:
                    # Load the FASTQ records
                    fastq_records = load_fastq_records(
                        gz=gz,
                        paired=True,
                        forward=True
                    )
                with gzip.open(reverse_bait, 'rt') as gz:
                    # Load the FASTQ records
                    fastq_records.update(load_fastq_records(
                        gz=gz,
                        paired=True,
                        forward=False
                    ))
            else:
                fastq_records = None
                logging.debug(
                    'Not loading Nanopore FASTQ records into memory; '
                    'workers will use index_db on disk'
                )
        else:
            if threads == 1:
                with gzip.open(unpaired_bait, 'rt') as gz:
                    # Load the FASTQ records
                    fastq_records = load_fastq_records(
                        gz=gz,
                        paired=False,
                        forward=True
                    )
            else:
                fastq_records = None
                logging.debug(
                    'Not loading Nanopore FASTQ records into memory; '
                    'workers will use index_db on disk'
                )

    logging.info('Detecting contamination...')

    # Now do mapping in two steps - first, map reads back to database with
    # ambiguous reads matching all - this will be used to get a count of
    # number of reads aligned to each gene/allele so we can create a custom
    # rMLST file with only the most likely allele for each gene.
    kma_report = os.path.join(
        sample_tmp_dir,
        f'{sample_name}_kma'
    )

    # Set the KMA database name
    kma_database = sample_database.replace('.fasta', '') + '_kma'

    # Index the database if necessary
    index_databases(sample_database=sample_database)

    # Run KMA.
    if paired:
        if not os.path.isfile(kma_report + '.res'):
            cmd = (
                f'kma -ipe {forward_trimmed} {reverse_trimmed} '
                f'-t_db {kma_database} -o {kma_report} -t {threads}'
            )
            out, err = run_cmd(cmd=cmd)

            # Write to logfile
            write_to_logfile(
                logfile=log,
                out=out,
                err=err,
                cmd=cmd
            )

    # Unpaired reads
    else:
        if not os.path.isfile(kma_report + '.res'):
            if data_type == 'Illumina':
                # Use the FASTA file (rather than the reads) as the input
                if fasta:
                    cmd = (
                        f'kma -i {pair[0]} -t_db {kma_database} -mem_mode '
                        f'-ID 100 -ConClave 2 -ex_mode -o {kma_report} '
                        f'-t {threads}'
                    )
                else:
                    cmd = (
                        f'kma -i {unpaired_trimmed} -t_db {kma_database} '
                        f'-o {kma_report} -t {threads}'
                    )
            else:
                # Recommended Nanopore settings from KMA repo:
                # https://bitbucket.org/genomicepidemiology/kma
                cmd = (
                    f'kma -i {unpaired_bait} -t_db {kma_database} '
                    f'-o {kma_report} -mem_mode -mp 20 -mrs 0.0 -bcNano '
                    f'-t {threads}'
                )
            out, err = run_cmd(cmd=cmd)
            write_to_logfile(
                logfile=log,
                out=out,
                err=err,
                cmd=cmd
            )

    # Set the rMLST report path
    rmlst_report = os.path.join(output_folder, sample_name + '_alleles.tsv')

    # Parse the KMA report to find the best allele for each rMLST gene
    gene_alleles = find_rmlst_type(
        kma_report=kma_report + '.res',
        rmlst_report=rmlst_report
    )

    # Set the rMLST FASTA path
    rmlst_fasta = os.path.join(
        sample_tmp_dir,
        f'{sample_name}_alleles.fasta'
    )

    # Check if the rMLST FASTA file already exists
    if not os.path.isfile(rmlst_fasta):
        # Create a FASTA file that has only one allele per rMLST gene
        with open(rmlst_fasta, 'w', encoding='utf-8') as f:
            # Use SeqIO to parse the database and write out only the
            # relevant alleles
            for contig in SeqIO.parse(sample_database, 'fasta'):
                if contig.id in gene_alleles:
                    SeqIO.write(contig, f, 'fasta')

    # Get total gene length for later reporting
    rmlst_gene_length = find_total_sequence_length(
        fasta_file=rmlst_fasta
    )

    logging.debug('Total gene length is %s', rmlst_gene_length)

    # Initialize contamination boolean
    pysam_pass = True
    # Second step of mapping - Do a mapping of our baited reads against a
    # fasta file that has only one allele per rMLST gene.
    try:
        # Index the rMLST FASTA file if necessary
        if not os.path.isfile(rmlst_fasta + '.fai'):
            pysam.faidx(rmlst_fasta)

        # Set the path to the contamination BAM, sorted BAM, and SAM files
        outbam = os.path.join(
            sample_tmp_dir,
            f'{sample_name}_contamination.bam'
        )
        sorted_bam = os.path.join(
            sample_tmp_dir,
            f'{sample_name}_contamination_sorted.bam'
        )

        # Run the mapping if the sorted BAM doesn't already exist
        if not os.path.isfile(sorted_bam):
            # Perform mapping differently for paired and unpaired reads
            if paired:
                cmd = (
                    f'bbmap.sh ref={rmlst_fasta} in={forward_trimmed} '
                    f'in2={reverse_trimmed} out={outbam} threads={threads} '
                    'mdtag nodisk'
                )
                # Add subfilter for cgMLST databases
                if cgmlst_db is not None:
                    # Lots of core genes seem to have relatives within a
                    # genome that are at ~70% identity. This means that reads
                    # that shouldn't map do, and cause false positives. Adding
                    # in this sub-filter means that reads can only have one
                    # mismatch, so they have to be from the right gene for
                    # this to work.
                    cmd += ' subfilter=1'
                # Add in memory string if specified
                if xmx:
                    cmd += f' -Xmx{xmx}'

                # Run the command
                out, err = run_cmd(cmd=cmd)
                write_to_logfile(
                    logfile=log,
                    out=out,
                    err=err,
                    cmd=cmd
                )
            else:
                # If Illumina FASTQ:
                if data_type == 'Illumina' and not fasta:
                    cmd = (
                        f'bbmap.sh ref={rmlst_fasta} in={unpaired_trimmed} '
                        f'out={outbam} threads={threads} mdtag nodisk'
                    )

                    # Add subfilter for cgMLST databases
                    if cgmlst_db is not None:
                        # Core genes can have relatives within a genome that
                        # are at ~70 percent identity. This means that reads
                        # that shouldn't map do, and cause false positives.
                        # Adding in this sub-filter means reads can only have
                        # one mismatch, so they have to be from the right gene
                        # for this to work.
                        cmd += ' subfilter=1'
                    # Add in memory string if specified
                    if xmx:
                        cmd += f' -Xmx{xmx}'

                    # Run the command
                    out, err = run_cmd(cmd=cmd)
                    write_to_logfile(
                        logfile=log,
                        out=out,
                        err=err,
                        cmd=cmd
                    )
                else:
                    if fasta:
                        ax = 'asm5'
                    # If Nanopore FASTQ:
                    else:
                        ax = 'map-ont'
                    # Use minimap2 for Nanopore data or FASTA data
                    cmd = (
                        f'minimap2 --MD -t {threads} -ax {ax} {rmlst_fasta} '
                        f'{unpaired_bait}'
                    )
                    # Convert SAM to sorted BAM in one step to save time/disk
                    # space
                    cmd += (
                        f' | samtools view -@ {threads} -h -bT {rmlst_fasta} -'
                        f' | samtools sort - -@ {threads} -o {sorted_bam}'
                    )

                    # Run the command
                    out, err = run_cmd(cmd=cmd)
                    write_to_logfile(
                        logfile=log,
                        out=out,
                        err=err,
                        cmd=cmd
                    )

        # Ensure the BAM is sorted and indexed
        if not os.path.isfile(sorted_bam):
            # Use pysam to sort and index the BAM
            pysam.sort('-o', sorted_bam, outbam)
        if not os.path.isfile(sorted_bam + '.bai'):
            pysam.index(sorted_bam)

        # Now find number of multi-positions for each rMLST gene/allele
        # combination
        multi_positions = 0

        # Prefer an in-memory ThreadPool when FASTQ records are available in
        # the core process to avoid pickling large dicts. Otherwise ensure
        # index DBs exist (build once in core to avoid worker contention),
        # then spawn pool with per-worker index initializers.
        use_threadpool = fastq_records is not None
        if use_threadpool:
            logging.info(
                'Using in-memory ThreadPool for sample %s (fastq records '
                'loaded in core process).',
                sample_name
            )
            p = ThreadPool(processes=threads)
        else:
            _ensure_fastq_index(
                forward_trimmed if paired else unpaired_trimmed,
                reverse_trimmed if paired else None,
                paired,
                sample_tmp_dir
            )
            p = multiprocessing.Pool(
                processes=threads,
                initializer=_init_fastq_index,
                initargs=(
                    forward_trimmed if paired else unpaired_trimmed,
                    reverse_trimmed if paired else None, paired, sample_tmp_dir
                )
            )

        # Build argument lists for parallel processing
        nanopore = True if data_type == 'Nanopore' else False
        nanopore_list = [nanopore] * len(gene_alleles)
        bamfile_list = [sorted_bam] * len(gene_alleles)
        reference_fasta_list = [rmlst_fasta] * len(gene_alleles)
        fasta_list = [fasta] * len(gene_alleles)
        quality_cutoff_list = [quality_cutoff] * len(gene_alleles)
        min_quality_list = [min_quality] * len(gene_alleles)
        base_cutoff_list = [base_cutoff] * len(gene_alleles)
        base_fraction_list = [base_fraction_cutoff] * len(gene_alleles)
        fastq_records_list = [fastq_records] * len(gene_alleles)
        error_cutoff_list = [error_cutoff] * len(gene_alleles)

        logging.debug(
            'Processing %s alleles for sample %s (parallelized)',
            len(gene_alleles), sample_name
        )

        # Initialize result lists
        multibase_dict_list = []
        report_write_list = []
        if threads == 1:
            for i, gene in enumerate(gene_alleles):
                multibase_dict, report_write = read_contig(
                    contig_name=gene,
                    bamfile_name=bamfile_list[i],
                    reference_fasta=reference_fasta_list[i],
                    fastq_records=fastq_records_list[i],
                    quality_cutoff=quality_cutoff_list[i],
                    min_quality=min_quality_list[i],
                    base_cutoff=base_cutoff_list[i],
                    base_fraction_cutoff=base_fraction_list[i],
                    fasta=fasta_list[i],
                    nanopore=nanopore_list[i],
                    error_cutoff=error_cutoff_list[i]
                )
                multibase_dict_list.append(multibase_dict)

                # Keep representation consistent with parallel branch: a
                # list of lines per gene
                report_write_list.append(report_write.splitlines(True))
        else:
            # Filter gene list to those present in the BAM header to avoid
            # pysam errors during pileup
            try:
                with pysam.AlignmentFile(sorted_bam, 'rb') as bam_for_check:
                    bam_refs = set(
                        x.decode() if isinstance(x, bytes) else str(x)
                        for x in bam_for_check.references
                    )
            except (
                OSError,
                FileNotFoundError,
                ValueError,
                AttributeError,
                UnicodeDecodeError,
                SamtoolsError
            ) as exc:
                logging.debug(
                    'Could not read BAM references for %s: %s',
                    sorted_bam, exc
                )
                bam_refs = set()

            # Normalize gene IDs and filter by BAM header refs
            gene_alleles = [
                g.decode() if isinstance(
                    g, bytes
                ) else str(g) for g in gene_alleles
            ]
            present_genes = [g for g in gene_alleles if g in bam_refs]
            missing_genes = [g for g in gene_alleles if g not in bam_refs]

            # Log any missing genes
            if missing_genes:
                logging.warning(
                    'Skipping %s alleles not present in BAM header for '
                    'sample %s (first 20): %s',
                    len(missing_genes), sample_name, missing_genes[:20]
                )
            if not present_genes:
                logging.warning(
                    'No genes present in BAM header for sample %s; '
                    'aborting analysis for this sample.',
                    sample_name
                )

                # Clean up and return early for this sample
                if keep_files is False:
                    try:
                        shutil.rmtree(sample_tmp_dir)
                    except FileNotFoundError:
                        # Directory already removed — ignore
                        pass
                    except OSError as exc:
                        # Log failure to remove dir (permission/IO errors)
                        logging.debug(
                            'Failed to remove temp dir %s: %s',
                            sample_tmp_dir, exc
                        )

                # Write empty report
                write_output(
                    output_report=os.path.join(
                        output_folder, 'confindr_report.tsv'
                    ),
                    sample_name=sample_name,
                    multi_positions=0,
                    genus=genus,
                    total_gene_length=0,
                    database_download_date=database_download_date
                )
                return

            # Build contig-length map using FASTA index (.fai) where possible
            contig_lengths = {}
            try:
                fai_path = reference_fasta_list[0] + '.fai'
                if os.path.isfile(fai_path):
                    with open(fai_path, 'r', encoding='utf-8') as fh:
                        for line in fh:
                            cols = line.strip().split('\t')
                            if len(cols) >= 2:
                                contig_lengths[cols[0]] = int(cols[1])
                else:
                    for rec in SeqIO.parse(reference_fasta_list[0], 'fasta'):
                        contig_lengths[rec.id] = len(rec.seq)
            except (
                OSError,
                ValueError,
                AttributeError,
                UnicodeDecodeError
            ) as exc:
                logging.debug(
                    'Could not obtain contig lengths from FASTA/.fai: %s',
                    exc
                )

            # Partition contigs into balanced chunks (first-fit decreasing)
            chunks = _build_contig_chunks(
                contig_lengths=contig_lengths,
                contigs=present_genes,
                threads=threads,
                multiplier=CONTIG_CHUNK_MULTIPLIER,
                max_chunk_bases=CONTIG_CHUNK_MAX_BASES
            )

            logging.info(
                'Partitioned %s contigs into %s chunks (multiplier=%s, '
                'max_chunk_bases=%s)',
                len(present_genes),
                len(chunks),
                CONTIG_CHUNK_MULTIPLIER,
                CONTIG_CHUNK_MAX_BASES
            )

            # Build kwargs list (one chunk per task)
            kwargs_list = []
            for chunk in chunks:
                kwargs_list.append({
                    'contig_chunk': chunk,
                    'bamfile_name': sorted_bam,
                    'reference_fasta': rmlst_fasta,
                    'fastq_records': fastq_records,
                    'quality_cutoff': quality_cutoff,
                    'min_quality': min_quality,
                    'base_cutoff': base_cutoff,
                    'base_fraction_cutoff': base_fraction_cutoff,
                    'fasta': fasta,
                    'error_cutoff': error_cutoff,
                    'nanopore': True if data_type == 'Nanopore' else False,
                    'max_expected_positions': max_expected_positions
                })

            # Use map with the dispatch helper which unpacks kwargs
            try:
                total_contigs = len(present_genes)
                processed_contigs = 0
                tstart = time.time()
                log_every = max(1, total_contigs // 20)
                for multibase_dict, report_write in p.map(
                    _read_contig_chunk_dispatch, kwargs_list, chunksize=1
                ):
                    # multibase_dict contains possibly many contigs
                    for k, v in multibase_dict.items():
                        multibase_dict_list.append({k: v})
                        processed_contigs += 1
                    # Split combined report by contig and append per-contig
                    if report_write:
                        lines = report_write.splitlines(True)
                        for contig in multibase_dict.keys():
                            contig_lines = [
                                ln for ln in lines if ln.startswith(
                                    contig + '\t'
                                )
                            ]
                            report_write_list.append(contig_lines)

                    # Log progress
                    if (
                        processed_contigs % log_every == 0
                        or processed_contigs == total_contigs
                    ):
                        # Compute elapsed time, average per contig, and ETA
                        elapsed = time.time() - tstart
                        avg = elapsed / max(1, processed_contigs)
                        eta = avg * (total_contigs - processed_contigs)
                        logging.info(
                            'Sample %s progress: %s/%s (%.1f%%). Elapsed: %s, '
                            'ETA: %s',
                            sample_name,
                            processed_contigs,
                            total_contigs,
                            processed_contigs / total_contigs * 100.0,
                            _format_seconds(s=elapsed),
                            _format_seconds(s=eta)
                        )
            except (
                SamtoolsError,
                OSError,
                ValueError,
                KeyError,
                IndexError,
                RuntimeError
            ) as exc:
                logging.exception(
                    'Pool processing failed for sample %s; terminating '
                    'workers: %s',
                    sample_name, exc
                )
                try:
                    p.terminate()
                except (OSError, RuntimeError) as exc_term:
                    logging.debug('Failed to terminate pool: %s', exc_term)
                try:
                    p.join()
                except (OSError, RuntimeError) as exc_join:
                    logging.debug('Failed to join pool: %s', exc_join)
                raise
            finally:
                try:
                    p.close()
                    p.join()
                except (OSError, RuntimeError) as exc:
                    logging.debug('Pool close/join failed: %s', exc)

    except SamtoolsError:
        # Handle Samtools errors gracefully
        pysam_pass = False
        multi_positions = 0
        multibase_dict_list = []
        report_write_list = []

    # Write out report info.
    report_file = os.path.join(
        output_folder,
        sample_name + '_contamination.tsv'
    )

    # Write the contamination report TSV
    with open(report_file, 'w', encoding='utf-8') as r:
        # Contamination report TSV header (tab-separated columns):
        # Gene, Position, RefBase, CongruentSNVs, TotalSNVs, ForwardSNVs,
        # ReverseSNVs, SNVCoverage, TotalCoverage, BaseCutoff, ErrorPercent,
        # PValue, AdjPValue, StrandP, PosP, MeanQual, MeanMapQ
        # Note: the per-position output from `position_details` concatenates
        # the four SNV "read-type" fields (Congruent/Paired/Forward/Reverse)
        # into a single field separated by commas — this means consumers that
        # index columns by numeric position should be careful (prefer matching
        # headers by name when possible).
        r.write(
            'Gene\tPosition\tRefBase\tCongruentSNVs\tTotalSNVs\tForwardSNVs\t'
            'ReverseSNVs\tSNVCoverage\tTotalCoverage\tBaseCutoff\tErrorPercent'
            '\tPValue\tAdjPValue\tStrandP\tPosP\tMeanQual\tMeanMapQ\n'
        )

        # Iterate through the report_write_list and write each line. Items may
        # be lists (parallel branch) or strings (debug branch); normalize to
        # a list of lines before writing.
        for item in report_write_list:
            lines = item if isinstance(item, list) else item.splitlines(True)
            for contamination_info in lines:
                r.write(contamination_info)

    # Total up the number of multibase positions using helper (excludes meta
    # keys)
    multi_positions = count_multibase_positions(
        multibase_dict_list=multibase_dict_list
    )

    # Determine SNP cutoff based on database type
    if cgmlst_db is None:
        snp_cutoff = math.ceil(rmlst_gene_length / 10000) + 1
    elif fasta:
        snp_cutoff = 1
    else:
        snp_cutoff = 10

    # Compute a simple per-sample score: the number of supported
    # multibase positions. This aligns the probabilistic model with the
    # sample-level decision rule.
    sample_score = float(multi_positions)

    # Write gene-summary TSV file
    gene_summary_file = os.path.join(
        output_folder,
        sample_name + '_gene_summary.tsv'
    )
    with open(gene_summary_file, 'w', encoding='utf-8') as gf:
        gf.write('Gene\tCombinedP\tGeneScore\tNumSigPositions\tNumPositions\n')
        for multibase_dict in multibase_dict_list:
            for gene, content in multibase_dict.items():
                # Skip meta keys
                if gene.startswith('_'):
                    continue

                # Get stats
                stats = content.get('_gene_stats', {})

                # Skip if no stats
                if not stats:
                    continue

                # Extract relevant stats
                combined_p = stats.get('combined_p')
                gene_score_val = stats.get('gene_score')
                num_sig_positions = stats.get('num_sig_positions', 0)
                num_positions = sum(
                    1 for k in content.keys() if not str(k).startswith('_')
                )
                combined_p_str = (
                    f'{combined_p:0.3e}' if combined_p is not None else 'ND'
                )
                gene_score_str = (
                    f'{gene_score_val:0.3f}' if gene_score_val is not None
                    else 'ND'
                )
                gf.write(
                    f'{gene}\t{combined_p_str}\t{gene_score_str}\t'
                    f'{num_sig_positions}\t{num_positions}\n'
                )

    logging.info(
        'Done! Number of contaminating SNVs found: %s\n', multi_positions
    )

    # Write summary report
    write_output(
        output_report=os.path.join(output_folder, 'confindr_report.tsv'),
        sample_name=sample_name,
        multi_positions=multi_positions,
        genus=genus,
        total_gene_length=rmlst_gene_length,
        snp_cutoff=snp_cutoff,
        database_download_date=database_download_date,
        pysam_pass=pysam_pass,
        sample_score=sample_score,
        use_probabilistic=use_prob_scoring,
        score_threshold=score_threshold
    )

    # Clean up temporary files unless keep_files is set
    if keep_files is False:
        shutil.rmtree(sample_tmp_dir)


def write_output(
    *,  # Enforce keyword-only arguments
    output_report: str,
    sample_name: str,
    multi_positions: int,
    genus: str,
    total_gene_length: int,
    database_download_date: str,
    snp_cutoff: int = 3,
    pysam_pass: bool = True,
    sample_score: Optional[float] = None,
    use_probabilistic: bool = False,
    score_threshold: Optional[float] = None,
    subreplicate_counts: Optional[List[int]] = None
) -> None:
    """
    Write ConFindr summary report

    Args:
        output_report: Path to CSV/TSV output report file. If it does not
        exist, a header will be written automatically.
        sample_name: Name of the sample.
        multi_positions: Number of multibase positions found for the sample.
        genus: Genus string (may include ':' when multiple genera are present).
        total_gene_length: Total bases examined across genes.
        database_download_date: Date string indicating database download date.
        *  # Enforce keyword-only arguments for the following params
        snp_cutoff: Number of cSNVs required to call a sample contaminated.
        pysam_pass: Whether pysam mapping completed successfully.
        sample_score: Optional per-sample score (diagnostic), now equal
            to the number of statistically supported positions.
        use_probabilistic: If True and score_threshold supplied, decide
            contamination based on sample_score >= score_threshold.
        score_threshold: Minimum number of supported positions required
            when use_probabilistic is True.

    Returns:
        None
    """
    # If the report file hasn't been created, make it, with appropriate header.
    if not os.path.isfile(output_report):
        with open(os.path.join(output_report), 'w', encoding='utf-8') as f:
            # Summary report is a TSV: Sample\tGenus\tNumContamSNVs\t
            # ContamStatus\tBasesExamined\tDatabaseDownloadDate\tScore
            # Place subreplicate summary columns immediately after ContamStatus
            # Primary summary columns: mean, median, stddev of
            # per-subreplicate SNV counts
            # Write header depending on whether subreplicates were run
            if subreplicate_counts and len(subreplicate_counts) > 0:
                # Primary summary columns: mean, median, stddev of
                # per-subreplicate SNV counts
                base_header = (
                    'Sample\tGenus\tMeanContamSNVs\tMedianContamSNVs\t'
                    'StdDevContamSNVs\tContamStatus'
                )
                subparts = ['NumSubreplicates'] + [
                    f'Subreplicate_{i+1}' for i in range(
                        len(
                            subreplicate_counts
                        )
                    )
                ]
                f.write(
                    base_header + '\t' + '\t'.join(
                        subparts
                    ) + '\tBasesExamined\tDatabaseDownloadDate\tScore\n'
                )
            else:
                # Default header for single-run analyses
                base_header = (
                    'Sample\tGenus\tNumContamSNVs\tContamStatus\t'
                    'BasesExamined\tDatabaseDownloadDate\tScore'
                )
                f.write(base_header + '\n')

    # Determine contamination status
    if pysam_pass:
        # Check contamination based on probabilistic or deterministic method
        if use_probabilistic and score_threshold is not None:
            # Probabilistic decision based on supported-position count
            contaminated = bool(
                sample_score is not None and sample_score >= score_threshold
            )
        else:
            contaminated = bool(
                multi_positions >= snp_cutoff or len(genus.split(':')) > 1
            )
    else:
        contaminated = 'Pysam SamtoolsError'
        multi_positions = 'ND'

    # Write out the summary line
    score_str = f'{sample_score:0.3f}' if sample_score is not None else 'ND'
    with open(output_report, 'a+', encoding='utf-8') as f:
        # If subreplicates were run, show mean/median/stddev as primary fields
        if subreplicate_counts and len(subreplicate_counts) > 0:
            # Compute mean, median, stddev
            n_sub = len(subreplicate_counts)
            mean_contam = sum(subreplicate_counts) / float(n_sub)
            sorted_vals = sorted(subreplicate_counts)

            # Median calculation
            if n_sub % 2 == 1:
                median_contam = float(sorted_vals[n_sub // 2])
            else:
                median_contam = (
                    sorted_vals[n_sub // 2 - 1] + sorted_vals[n_sub // 2]
                ) / 2.0

            # Standard deviation calculation
            sd_contam = pstdev(subreplicate_counts)

            # Base line (mean, median, stddev)
            base_line = (
                f'{sample_name}\t{genus}\t{mean_contam:0.2f}\t'
                f'{median_contam:0.2f}\t{sd_contam:0.2f}\t{contaminated}'
            )

            # Subreplicate counts
            subparts = [str(n_sub)] + [
                str(int(c)) for c in subreplicate_counts
            ]

            # Rest of the line
            rest = (
                f'\t{total_gene_length}\t{database_download_date}\t{score_str}'
            )
            f.write(base_line + '\t' + '\t'.join(subparts) + rest + '\n')
        else:
            # Single-run (no subreplicates): write legacy header/line with
            # NumContamSNVs
            base_line = (
                f'{sample_name}\t{genus}\t{multi_positions}\t{contaminated}\t'
                f'{total_gene_length}\t{database_download_date}\t{score_str}'
            )
            f.write(base_line + '\n')


def check_for_databases_and_download(
    *,  # Enforce keyword-only arguments
    database_location: str
) -> None:
    """
    Check for necessary ConFindr databases, download if not present.

    Args:
        database_location: Path to ConFindr database folder.

    Returns:
        None
    """
    # Check for the files necessary - should have rMLST_combined.fasta,
    # gene_allele.txt, profiles.txt, and refseq.msh
    necessary_files = [
        'Escherichia_db_cgderived.fasta',
        'Listeria_db_cgderived.fasta',
        'Salmonella_db_cgderived.fasta',
        'refseq.msh'
    ]

    # Also check for optional rMLST files
    optional_files = [
        'rMLST_combined.fasta',
        'gene_allele.txt',
        'profiles.txt'
    ]

    # Initialize flag
    all_files_present = True

    # Iterate through necessary files
    for necessary_file in necessary_files:
        if not os.path.isfile(os.path.join(database_location, necessary_file)):
            logging.warning('Could not find %s', necessary_file)
            all_files_present = False

    # Download necessary files if not present
    if not all_files_present:
        logging.warning(
            'Databases not present - downloading basic databases now...'
        )

        # Create database directory if it doesn't exist
        os.makedirs(database_location, exist_ok=True)

        # Download necessary database files
        download_mash_sketch(output_folder=database_location)
        download_cgmlst_derived_data(output_folder=database_location)

    # Check for optional files
    optional_files_present = True
    for optional_file in optional_files:
        if not os.path.isfile(os.path.join(database_location, optional_file)):
            optional_files_present = False
    if not optional_files_present:
        logging.warning(
            'Did not find rMLST databases, if you want to use ConFindr on '
            'genera other than Listeria, Salmonella, and Escherichia, '
            'you\'ll need to download them. Instructions are available at '
            'https://olc-bioinformatics.github.io/ConFindr/install/'
            '#downloading-confindr-databases\n'
        )


def check_valid_base_fraction(
    *,  # Enforce keyword-only arguments
    base_fraction: Optional[float]
) -> bool:
    """
    Validate that a base fraction is None or between 0 and 1 inclusive.

    Args:
        base_fraction: The fraction to validate (or None).

    Returns:
        True if valid, False otherwise.
    """
    if base_fraction is None:
        return True
    if 0 <= base_fraction <= 1:
        return True

    return False


def check_acceptable_xmx(
    *,  # Enforce keyword-only arguments
    xmx_string: str
) -> bool:
    """
    Validate a BBTools -Xmx memory specification string.

    Args:
        xmx_string: Memory string such as '20g', '800m', '1024K'.

    Returns:
        True if acceptable, False otherwise.
    """
    # Initialize flag
    acceptable_xmx = True

    # Set of acceptable suffixes
    acceptable_suffixes = ['K', 'M', 'G']

    # Check suffix
    if xmx_string[-1].upper() not in acceptable_suffixes:
        acceptable_xmx = False
        logging.error(
            'ERROR: Memory must be specified as K (kilobytes), M (megabytes), '
            'or G (gigabytes). Your specified suffix was %s.', xmx_string[-1]
        )

    # Check that the rest is an integer
    if '.' in xmx_string:
        acceptable_xmx = False
        logging.error(
            'ERROR: Xmx strings must be integers, floating point numbers '
            'are not accepted.'
        )

    # Check that the rest is digits
    if not str.isdigit(xmx_string[:-1]):
        acceptable_xmx = False
        logging.error(
            'ERROR: The amount of memory requested was not an integer.'
        )

    return acceptable_xmx


def get_version() -> str:
    """
    Get ConFindr version string.

    Returns:
        Version string.
    """
    return f'ConFindr {__version__}'


def _valid_downsample_depth(value: str) -> int:
    """
    argparse type-checker for --downsample_depth.
    Accepts integers 10..100 (inclusive)
    """
    try:
        v = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            'Downsample depth must be an integer'
        ) from exc

    # Check range
    if v < 10 or v > 100:
        raise argparse.ArgumentTypeError(
            'Downsample depth must be between 10 and 100'
        )

    return v


def recommend_xmx(
    *,
    fraction: float = 0.8,
    round_gb: bool = True
) -> str:
    """
    Calculate recommended -Xmx value based on available system memory.

    Args:
        fraction: Fraction of available memory to recommend (default 0.8).
        round_gb: If True, round down to nearest integer GB (default True).

    Returns:
        Recommended -Xmx string (e.g. '16g').
    """
    # Use psutil to calculate available memory
    mem_bytes = psutil.virtual_memory().available

    # Compute recommended bytes and convert to GB
    recommended_bytes = int(mem_bytes * fraction)
    recommended_gb = recommended_bytes / (1024**3)

    # Round down to nearest integer GB if specified
    if round_gb:
        recommended_gb = int(recommended_gb)

    return f'{int(recommended_gb)}g'
