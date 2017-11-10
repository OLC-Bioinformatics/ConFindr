# ConFindr

This program is designed to find bacterial intra-species contamination in raw NGS data. It does this
 by looking for multiple copies of rMLST genes, which are known to be universal across the bacterial kingdom
 and present only in single copies. 

### Program Requirements
- bbmap (>= 37.23) installed and present on your $PATH
- jellyfish (>= 2.2.6) installed and on your $PATH
- Python 3.5
- NCBI BLAST+ (>=2.2.31) 
- mash (>=2.0 - https://github.com/marbl/Mash) installed and on your $PATH


### Installation

- Clone this repository - you'll need the `databases` folder (Git LFS is required for this download.)
- Install the confindr script with pip: `pip install confindr`. With this done, you'll be able to call the script by 
simply typing `confindr.py` on the command line, no matter what directory you're currently in.


### Usage
- Program takes a folder with paired or single-ended fastq files as input. Files can be uncompressed, or compressed with gzip/bzip2.
- outputs results to a csv file which you name when calling the script.
- New_Detector.py is what is responsible for running analysis.

#### Example Usages

Detect contamination on any fastq files within fastq folder, outputs results to outputname.csv, and uses files in the 
databases folder for contamination detection - this folder should be the databases folder in this repository.

`confindr.py Fastq_Folder outputname databases`

#### Options

```
usage: confindr.py [-h] [-t THREADS] [-n NUMBER_SUBSAMPLES] [-k KMER_SIZE]
                       [-s SUBSAMPLE_DEPTH] [-c KMER_CUTOFF]
                       fastq_directory output_name databases

positional arguments:
  fastq_directory       Folder that contains fastq files you want to check for
                        contamination. Will find any fastq file that contains
                        .fq or .fastq in the filename.
  output_name           Base name for output/temporary directories.
  databases             Databases folder. Should contain rMLST_combined.fasta,
                        profiles.txt, and refseq.msh as well as
                        RefSeqSketchesDefaults.msh

optional arguments:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS
                        Number of threads to run analysis with.
  -n NUMBER_SUBSAMPLES, --number_subsamples NUMBER_SUBSAMPLES
                        Number of times to subsample.
  -k KMER_SIZE, --kmer-size KMER_SIZE
                        Kmer size to use for contamination detection.
  -s SUBSAMPLE_DEPTH, --subsample_depth SUBSAMPLE_DEPTH
                        Depth to subsample to. Higher increases sensitivity,
                        but also false positive rate. Default is 20.
  -c KMER_CUTOFF, --kmer_cutoff KMER_CUTOFF
                        Number of times you need to see a kmer before it is
                        considered trustworthy. Kmers with counts below this
                        number will be discarded.
```