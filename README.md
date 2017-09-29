# KmerContam

### Program Requirements
- bbmap (>= 37.23) installed and present on your $PATH
- jellyfish (>= 2.2.6) installed and on your $PATH
- Python 3.5 (2.7 should also work, though that isn't tested.)
- NCBI BLAST+ (>=2.2.31) 
- mash (>=1.1.1 - https://github.com/marbl/Mash) installed and on your $PATH

### Python Package Requirements
- pysam >= 0.11.2.2
- OLCTools >= 0.2.3
- biopython >= 1.70

### Usage
- Program takes a folder with paired or single-ended fastq files as input. Files can be uncompressed, or compressed with gzip/bzip2.
- outputs results to a csv file which you name when calling the script.
- Currently runs through New_Detector.py

#### Options
- Threads (-t): Number of threads to run analysis on. Default is number of cores on your system.
- Number subsamples (-n): Number of times to subsample. More is generally better, although more time-consuming. Default is 5.
- Kmer size (-k): Kmer size to use. Default is 31. Other values may make results extremely unreliable. USE WITH CAUTION.
- Subsample depth (-s): Depth to subsample to. Default is 20. Higher values increase sensitivity, but also false positive rate.
- Kmer cutoff (-c): Number of times a kmer has to be seen before it's considered to be trustworthy. Default value 2.
#### Example Usages

Detect contamination on any fastq files within fastq folder, outputs results to outputname.csv, and uses database.fasta
as the database of rMLST genes.

`python3 New_Detector.py Fastq_Folder outputname database.fasta`

## Using Docker
- A docker image will be created for ease of use. It can be downloaded at (insert here eventually where you'll put it).
- Add in instructions on mounting the files you want and running the docker image here.
