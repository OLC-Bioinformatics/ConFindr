# KmerContam

### Program Requirements
- bbmap (>= 37.23) installed and present on your $PATH
- jellyfish (>= 2.2.6) installed and on your $PATH
- Python 3.5 (2.7 should also work - not fully tested)
- CLARK (>=1.2.3) installed and on your $PATH

### Python Package Requirements
- pysam >= 0.11.2.2
- OLCTools >= 0.1.16

### Usage
- Program takes a folder with paired or single-ended fastq files as input. Files can be uncompressed, or compressed with gzip/bzip2.
- outputs results to a csv file which you name when calling the script.
- Currently runs through Detector.py

#### Options
- Classify (-c): If the -c flag is added, when suspected cross-species contamination is found CLARK will be run (light version) to try to identify which species are present. Default=False.
- Threads (-t): Number of threads to run analysis on. Default is number of cores on your system.
- Trim_fastq (-tr): Performs quality trimming using bbduk. Off by default, but should probably almost always be turned on.
- Remove\_bad_reads (-x): Removes reads which have contaminating kmers in them. Still highly experimental.

#### Example
python3 Detector.py Fastq_Folder outputname.csv

## Using Docker
- A docker image has also been created for ease of use. It can be downloaded at (insert here eventually where you'll put it).
- Add in instructions on mounting the files you want and running the docker image here.
