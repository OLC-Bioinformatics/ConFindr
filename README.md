# KmerContam

### Program Requirements
- bbmap (>= 37.23) installed and present on your $PATH
- jellyfish (>= 2.2.6) installed and on your $PATH
- Python 2.7
- CLARK (>=1.2.3) installed and on your $PATH

### Python Package Requirements
- jellyfish (see https://github.com/gmarcais/Jellyfish for installation instructions)
- pysam >= 0.11.2.2
- OLCTools >= 0.1.16

### Usage
- Program takes a folder with paired fastq files as input. Files can be uncompressed, or compressed with gzip/bzip2.
- outputs results to a csv file which you name when calling the script.
- Currently runs through test.py

#### Options
- Classify (-c): If the -c flag is added, when suspected cross-species contamination is found CLARK will be run (light version) to try to identify which species are present. Default=False.
- Threads (-t): Number of threads to run analysis on. Default is number of cores on your system.
- Trim_fastq (-tr): Performs quality trimming using bbduk. Off by default, but should probably almost always be turned on.

#### Example
python test.py Fastq_Folder outputname.csv
