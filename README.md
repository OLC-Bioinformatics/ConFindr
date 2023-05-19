[![CircleCI](https://dl.circleci.com/status-badge/img/gh/OLC-Bioinformatics/ConFindr/tree/main.svg?style=svg)](https://dl.circleci.com/status-badge/redirect/gh/OLC-Bioinformatics/ConFindr/tree/main)
[![PyPI version](https://badge.fury.io/py/confindr.svg)](https://badge.fury.io/py/confindr)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/confindr/README.html)


# ConFindr

This program is designed to find bacterial intra-species contamination in raw Illumina data. It does this
 by looking for multiple alleles of core, single copy genes.

For **complete instructions on installation and usage**, please visit [the ConFindr GitHub Pages site](https://olc-bioinformatics.github.io/ConFindr/).

## Important Note

ConFindr has only been validated using rMLST databases. **Please use them if possible** (`--rmlst`). Complete installation instructions can be found [here](https://olc-bioinformatics.github.io/ConFindr/install/#downloading-confindr-databases).

## Quickstart

### Installing ConFindr

1. Follow the instructions [here](https://bioconda.github.io/) to add the Bioconda channel to your list of conda channels, if it hasn't already been added.

2. Install ConFindr into a new conda environment named 'confindr':

`conda create -n confindr -c bioconda confindr=0.8.1`

3. Activate the new conda environment:

`conda activate confindr`

### Downloading and setting up the rMLST databases

Instructions for downloading and setting up the rMLST databases can be found [here](https://olc-bioinformatics.github.io/ConFindr/install/#downloading-confindr-databases).

### Testing ConFindr

1. To obtain an example dataset, run the following command, which will create a folder named `example-data` in your current working directory: 

`wget https://ndownloader.figshare.com/files/9972709 && tar xf 9972709 && rm 9972709`

2. As of version `0.7.0` ConFindr can be run automatically on _Escherichia_, _Salmonella_, and _Listeria_ with no further 
work on your part, using core-gene databases (*experimental*). Simply run:

`confindr -i example-data -o example-out`

3. To use the *recommended* rMLST database (after installation):

`confindr -i example-data -o example-out --rmlst`

Once ConFindr finishes running, take a look at the `confindr_report.csv` file found in `example-out`—it shows that multiple
alleles were found for many sites within the genes that ConFindr examines, meaning that this sample is quite contaminated!

If you want to run ConFindr on genera other than the 3 listed above, you'll need to get access to and download the rMLST databases by following the instructions [here](https://olc-bioinformatics.github.io/ConFindr/install/#downloading-confindr-databases).

## Running ConFindr in a Python Script

If you want to run ConFindr from within a script instead of running from the command line, here's how:

```python
from confindr_src import confindr

# Find read files.
paired_reads = confindr.find_paired_reads('path_to_fastq_folder', forward_id='_R1', reverse_id='_R2')
# Run confindr. This assumes that you have already downloaded the databases. If you haven't,
# you can run confindr.check_for_databases_and_download(database_location='path/where/you/want/to/download, tmpdir='a/tmp/dir')
for pair in paired_reads:
    confindr.find_contamination(pair=pair,
                                forward_id='_R1', # change if yours is different
                                threads=4, 
                                output_folder='path/to/output',
                                databases_folder='path/to/databases')
                                
```

## Reporting Issues

If you have any problems installing or running ConFindr, or have feature request,
please open an issue here on GitHub.

## Citing ConFindr

ConFindr has been published in PeerJ—if you use it in your work, please cite the following:

```
Low AJ, Koziol AG, Manninger PA, Blais B, Carrillo CD. 2019. ConFindr: rapid detection of intraspecies and cross-species contamination in bacterial whole-genome sequence data. PeerJ 7:e6995 https://doi.org/10.7717/peerj.6995
```