#!/usr/bin/env python3
"""
database_setup.py

Download and prepare ConFindr databases (rMLST and cgMLST-derived).
"""

# Standard library imports
from glob import glob
from typing import (
    Set
)
import argparse
import csv
import datetime
import json
import logging
import os
import re
import shutil
import ssl
import sys

# Third-party imports
from Bio import SeqIO
from Bio.Seq import Seq
try:
    from rauth import OAuth1Session
except ImportError:  # pragma: no cover
    OAuth1Session = None

# Local imports
from confindr_src.methods import (
    download_cgmlst_derived_data,
    download_mash_sketch,
    index,
)


class RmlstRest:
    """
    Class to interact with the rMLST REST API for downloading loci and
    scheme data.
    """

    def get_session_token(self):
        """
        Get a session token from the rMLST REST API.
        """
        # Create an OAuth1 session
        session_request = OAuth1Session(
            self.consumer_key,
            self.consumer_secret,
            access_token=self.access_token,
            access_token_secret=self.access_secret
        )

        # Set up the URL for getting a session token
        url = self.test_rest_url + '/oauth/get_session_token'

        # If unverified SSL, set verify=False
        if self.unverified:
            # Perform a GET request with verify=False
            r = session_request.get(url, verify=False)
        else:
            # Perform a GET request
            r = session_request.get(url)

        # If the status code is '200' (OK), proceed
        if r.status_code == 200:
            self.session_token = r.json()['oauth_token']
            self.session_secret = r.json()['oauth_token_secret']
        # If we couldn't get a session token, exit with error.
        else:
            logging.error(
                'ERROR: Couldn\'t get a session token for rMLST database '
                'download. Check that your consumer secret and access token '
                'files have valid credentials and try again.'
            )
            sys.exit(1)

    def get_loci_and_scheme_url(self):
        """
        Get the URLs for loci and scheme data from the rMLST REST API.
        """
        # Create an OAuth1 session
        session = OAuth1Session(
            self.consumer_key,
            self.consumer_secret,
            access_token=self.session_token,
            access_token_secret=self.session_secret
        )

        # Make the GET request to the test REST URL
        if self.unverified:
            r = session.get(self.test_rest_url, verify=False)
        else:
            r = session.get(self.test_rest_url)

        # Check if the request was successful
        if r.status_code in [200, 201]:
            # Decode the response based on content type
            if re.search('json', r.headers['content-type'], flags=0):
                decoded = r.json()
            else:
                try:
                    decoded = json.loads(r.text)
                except json.JSONDecodeError:
                    decoded = r.text

            # Extract the URLs from the returned data
            self.loci = decoded['loci']
            self.profile = decoded['schemes']
        # If we couldn't get the loci and scheme URLs, exit with error.
        else:
            logging.error(
                'ERROR: Could not find URLs for rMLST download, they may have '
                'moved. Please open an issue at '
                'https://github.com/OLC-Bioinformatics/ConFindr/issues '
                'and we\'ll get things sorted out.'
            )
            sys.exit(1)

    def download_loci(self):
        """
        Download all rMLST loci from the REST API.
        """
        # Create an OAuth1 session
        session = OAuth1Session(
            self.consumer_key,
            self.consumer_secret,
            access_token=self.session_token,
            access_token_secret=self.session_secret
        )

        # Make the GET request to the loci URL
        if self.unverified:
            r = session.get(self.loci, verify=False)
        else:
            r = session.get(self.loci)

        # Check if the request was successful
        if r.status_code in [200, 201]:
            # Decode the response based on content type
            if re.search('json', r.headers['content-type'], flags=0):
                decoded = r.json()
            else:
                decoded = r.text

            # Extract all the URLs in the decoded dictionary under the key loci
            for locus_url in decoded['loci']:
                # Set up output file path
                output_file = os.path.join(
                    self.output_folder,
                    f'{os.path.split(locus_url)[1]}.tfa'
                )

                logging.info('Downloading %s...', os.path.split(locus_url)[1])

                # Make the GET request to download the locus FASTA
                if self.unverified:
                    download = session.get(
                        locus_url + '/alleles_fasta',
                        verify=False
                    )
                else:
                    download = session.get(locus_url + '/alleles_fasta')

                # Check if the download was successful
                if download.status_code in [200, 201]:
                    # Decode based on content type
                    if re.search(
                        'json',
                        download.headers['content-type'],
                        flags=0
                    ):
                        decoded = download.json()
                    else:
                        decoded = download.text

                    # Write the locus FASTA to disk
                    with open(
                        output_file,
                        'w',
                        encoding='utf-8'
                    ) as locus_fasta:
                        locus_fasta.write(decoded)
        # If we couldn't get the loci, exit with error
        else:
            logging.error(
                'ERROR: Could not find URLs for rMLST download, they may have '
                'moved. Please open an issue at https://github.com/'
                'OLC-Bioinformatics/ConFindr/issues and we\'ll get things '
                'sorted out.'
            )
            sys.exit(1)

    def download_profile(self):
        """
        Download the rMLST profiles from the REST API.
        """
        # Set up output profile file path
        profile_file = os.path.join(self.output_folder, 'profiles.txt')

        # Create an OAuth1 session
        session = OAuth1Session(
            self.consumer_key,
            self.consumer_secret,
            access_token=self.session_token,
            access_token_secret=self.session_secret
        )

        # Make the GET request to download the profiles CSV
        if self.unverified:
            r = session.get(self.profile + '/1/profiles_csv', verify=False)
        else:
            r = session.get(self.profile + '/1/profiles_csv')

        logging.info('Downloading rMLST profiles...')

        # Check if the download was successful
        if r.status_code in [200, 201]:
            # Decode based on content type
            if re.search('json', r.headers['content-type'], flags=0):
                decoded = r.json()
            else:
                decoded = r.text

            # Write the profile file to disk
            with open(profile_file, 'w', encoding='utf-8') as profile:
                profile.write(decoded)

    def get_request_token(self):
        """
        Get a request token from the rMLST REST API.
        """
        # Create an OAuth1 session
        session = OAuth1Session(
            consumer_key=self.consumer_key,
            consumer_secret=self.consumer_secret
        )

        # Use the test URL in the GET request
        r = session.request(
            method='GET',
            url=self.request_token_url,
            params={'oauth_callback': 'oob'}
        )

        # If the status code is '200' (OK), proceed
        if r.status_code == 200:
            self.request_token = r.json()['oauth_token']
            self.request_secret = r.json()['oauth_token_secret']

    def get_access_token(self):
        """
        Get an access token from the rMLST REST API.
        """
        # Create the authorization URL
        authorize_url = (
            self.test_web_url + '&page=authorizeClient&oauth_token='
            + self.request_token
        )

        print('Visit this URL in your browser: ' + authorize_url)

        # Accept the oauth_verifier from the user
        verifier = input('Enter oauth_verifier from browser: ')
        session_request = OAuth1Session(
            consumer_key=self.consumer_key,
            consumer_secret=self.consumer_secret,
            access_token=self.request_token,
            access_token_secret=self.request_secret
        )

        # Perform a GET request with the appropriate keys and tokens
        if self.unverified:
            r = session_request.get(
                self.access_token_url,
                verify=False,
                params={
                    'oauth_verifier': verifier
                }
            )
        else:
            r = session_request.get(
                self.access_token_url,
                params={
                    'oauth_verifier': verifier
                }
            )

        # If the status code is '200' (OK), proceed
        if r.status_code == 200:
            # Save the JSON-decoded token secret and token
            self.access_token = r.json()['oauth_token']
            self.access_secret = r.json()['oauth_token_secret']

    def __init__(self, consumer_secret_file, output_folder, unverified=False):
        self.test_rest_url = 'https://rest.pubmlst.org/db/pubmlst_rmlst_seqdef'
        self.test_web_url = (
            'https://pubmlst.org/cgi-bin/bigsdb/bigsdb.pl?'
            'db=pubmlst_rmlst_seqdef'
        )
        self.request_token_url = (
            self.test_rest_url + '/oauth/get_request_token'
        )
        self.access_token_url = self.test_rest_url + '/oauth/get_access_token'
        self.authorize_url = self.test_web_url + '&page=authorizeClient'
        self.output_folder = output_folder
        self.unverified = unverified

        # Get the consumer secret set up.
        if not os.path.isfile(consumer_secret_file):
            logging.error(
                'ERROR: Could not find consumer secret file. Please make sure '
                'the file you specified (%s) exists and try again.',
                consumer_secret_file
            )
            sys.exit(1)

        # Read in the consumer key and secret from the file
        with open(consumer_secret_file, encoding='utf-8') as f:
            lines = f.readlines()
        try:
            self.consumer_key = lines[0].rstrip()
            self.consumer_secret = lines[1].rstrip()
        # Only two lines should be in the file
        except IndexError:
            logging.error(
                'ERROR: Could not parse your consumer secret file. File '
                'should have supplied consumer key on first line, and '
                'consumer secret on the second line.'
            )
            sys.exit(1)

        # Initialize other variables
        self.session_secret = str()
        self.session_token = str()
        self.loci = str()
        self.profile = str()
        self.request_token = str()
        self.request_secret = str()
        self.access_token = str()
        self.access_secret = str()


def create_gene_allele_file(
    *,  # Enforce keyword arguments
    profiles_file: str,
    gene_allele_file: str
) -> Set[str]:
    """
    Create a mapping of genera to alleles from an rMLST profiles file.

    Args:
        profiles_file: Path to TSV profiles file from rMLST.
        gene_allele_file: Output path to write genus:allele lists.

    Returns:
        Set of genera observed in the profiles file.
    """
    # Initialize dictionary to hold genus to allele mappings and set for genera
    genus_allele_info = {}
    genera = set()

    # Read in the profiles file
    with open(profiles_file, encoding='utf-8') as tsvfile:
        # Use DictReader to parse the TSV file
        reader = csv.DictReader(tsvfile, delimiter='\t')

        # Iterate through each row in the TSV
        for row in reader:
            genus = row['genus']
            # If the genus is uncertain e.g. Escherichia/Shigella, split on
            # the /, and use Escherichia as the genus
            if '/' in genus:
                genus = genus.split('/')[0]

            # Add genus to set of genera
            genera.add(genus)

            # If genus not already in dictionary, add it with empty list
            if genus not in genus_allele_info:
                genus_allele_info[genus] = []

            # Iterate through each of the 65 rMLST loci
            for i in range(1, 66):
                if i < 10:
                    gene = 'BACT00000' + str(i)
                else:
                    gene = 'BACT0000' + str(i)

                # If the gene is in the row, get the allele number
                if gene in row:
                    allele_number = row[gene]
                    gene_allele = f'{gene}_{allele_number}'

                    # If allele number is not 'N' and allele not already in
                    # list, add it to the list for that genus
                    if (
                        allele_number != 'N'
                        and gene_allele not in genus_allele_info[genus]
                    ):
                        genus_allele_info[genus].append(gene_allele)

    # Write the genus to allele mapping to the output file
    with open(gene_allele_file, 'w', encoding='utf-8') as f:
        for genus, allele_list in genus_allele_info.items():
            # Write genus:allele1,allele2,...
            f.write(str(genus) + ':')
            for allele in allele_list:
                f.write(str(allele) + ',')
            f.write('\n')

    return genera


def setup_confindr_database(
    *,  # Enforce keyword arguments
    output_folder: str,
    consumer_secret: str,
    index_databases: bool = False,
    unverified: bool = False
) -> None:
    """
    Set up ConFindr databases by downloading and preparing rMLST data.

    Args:
        output_folder: Path to download databases to.
        consumer_secret: Path to the consumer secret credentials file
        index_databases: If True, index genus-specific databases (slow).
        unverified: If True, disable HTTPS certificate verification

    Returns:
        None
    """
    # Go through the REST API in order to get profiles downloaded.
    rmlst_rest = RmlstRest(
        consumer_secret_file=consumer_secret,
        output_folder=output_folder,
        unverified=unverified
    )

    # Get request token, access token, session token, loci URL, and scheme URL
    rmlst_rest.get_request_token()
    rmlst_rest.get_access_token()
    rmlst_rest.get_session_token()
    rmlst_rest.get_loci_and_scheme_url()
    rmlst_rest.download_loci()
    rmlst_rest.download_profile()

    # With the sequences downloaded, make a file of all rMLST sequences
    logging.info('Combining rMLST files...')

    # Set the name of the combined rMLST FASTA file
    combined_fasta_file = os.path.join(output_folder, 'rMLST_combined.fasta')
    with open(combined_fasta_file, 'w', encoding='utf-8') as f:

        # Get a sorted list of all locus files
        locus_files = sorted(glob(os.path.join(output_folder, 'BACT*.tfa')))

        # Iterate through each locus file
        for locus_file in locus_files:
            for record in SeqIO.parse(locus_file, 'fasta'):
                record.id = record.id.replace('-', '_')

                # Normalise sequence without touching protected members.
                # Handles: gaps, Ns (both upper/lower), and byte-like reprs
                # (e.g. "b'ACGT...'") that some inputs may contain.
                seq_text = str(record.seq)

                # If the sequence looks like a byte-repr, strip the b'...'
                if seq_text.startswith("b'") and seq_text.endswith("'"):
                    seq_text = seq_text[2:-1]

                # Remove gaps and Ns (handle lowercase 'n' too), keep as str
                seq_text = seq_text.replace('-', '').replace('N', '').replace(
                    'n', ''
                )

                # Assign back using a Seq object (no protected access)
                record.seq = Seq(seq_text)

                # Clear out name and description fields
                record.name = ''
                record.description = ''

                # Write the cleaned record to the combined FASTA file
                SeqIO.write(record, f, 'fasta')

            # Clean up individual file.
            try:
                os.remove(locus_file)
            except OSError:
                logging.warning(
                    'WARNING: Could not delete %s. This won\'t affect '
                    'ConFindr performance, but  you may want to delete it to '
                    'save on disk space.', locus_file
                )

    logging.info('Assigning alleles to genera...')

    # Parse profiles so that we know what alleles are found with each genus.
    genera = create_gene_allele_file(
        profiles_file=os.path.join(output_folder, 'profiles.txt'),
        gene_allele_file=os.path.join(output_folder, 'gene_allele.txt')
    )

    # Index the rMLST database if specified
    if index_databases:
        index(
            output_folder=output_folder,
            genera=sorted(list(genera)),
            cgderived=False
        )


def main():
    """
    Main function for setting up ConFindr databases.
    """
    # Set up logging
    logging.basicConfig(
        format='\033[92m \033[1m %(asctime)s \033[0m %(message)s ',
        level=logging.INFO,
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # Set up argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-o', '--output_folder',
        default=os.environ.get(
            'CONFINDR_DB',
            os.path.expanduser(
                '~/.confindr_db'
            )
        ),
        help=(
            'Path to download databases to - if folder does not exist, will '
            'be created. If folder does exist, will be deleted and updated '
            'sequences downloaded. Defaults to ~/.confindr_db, or the '
            'CONFINDR_DB environmental variable.'
        )
    )
    parser.add_argument(
        '-s', '--secret_file',
        type=str,
        help='Path to consumer secret file for rMLST database.'
    )
    parser.add_argument(
        '-i', '--index_databases',
        action='store_true',
        help=(
            'Enable this option if you are installing the databases to a '
            'drive that will be read-only after the installation. The script '
            'will create and index all the necessary genus-specific database '
            'files. Note that this is very slow for the rMLST database.'
        )
    )
    parser.add_argument(
        '-u', '--unverified',
        action='store_true',
        help=(
            'Enable this option if you plan on running ConFindr behind a '
            'firewall and/or have a self-signed certificate. '
            'Adds \'verify=False\' during session requests.'
        )
    )

    # Parse the arguments
    args = parser.parse_args()

    # Create output folder, removing old one if it exists
    if os.path.isdir(args.output_folder):
        logging.info('Removing old databases...')
        shutil.rmtree(args.output_folder)
    os.makedirs(args.output_folder)

    # Set the context for unverified SSL if specified
    if args.unverified:
        ssl._create_default_https_context = ssl._create_unverified_context

    # Download the cgMLST-derived data
    download_cgmlst_derived_data(output_folder=args.output_folder)

    # Download and set up databases based on presence of secret file
    if args.secret_file is None:
        logging.warning(
            'WARNING: Without an rMLST secret file, data will only be '
            'downloaded for Escherichia, Salmonella, and Listeria. See '
            'https://olc-bioinformatics.github.io/ConFindr/install/#'
            'downloading-confindr-databases for instructions on how to get '
            'access to rMLST databases so ConFindr can be used for other '
            'species as well'
        )
    else:
        setup_confindr_database(
            output_folder=args.output_folder,
            consumer_secret=args.secret_file,
            index_databases=args.index_databases,
            unverified=args.unverified
        )
    download_mash_sketch(output_folder=args.output_folder)

    # Determine current date for download date file
    current_year = datetime.datetime.utcnow().year
    current_month = datetime.datetime.utcnow().month
    current_day = datetime.datetime.utcnow().day

    # Set the name of the download date file
    download_date_file = os.path.join(args.output_folder, 'download_date.txt')

    # Write the download date to a file
    with open(download_date_file, 'w', encoding='utf-8') as f:
        f.write(f'{current_year}-{current_month}-{current_day}')
    logging.info('Done downloading ConFindr databases!')


if __name__ == '__main__':
    main()
