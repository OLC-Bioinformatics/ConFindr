[![Build status](https://travis-ci.org/OLC-Bioinformatics/ConFindr.svg?master)](https://travis-ci.org/OLC-Bioinformatics)
[![PyPI version](https://badge.fury.io/py/confindr.svg)](https://badge.fury.io/py/confindr)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/confindr/README.html)

# ConFindr

This program is designed to find bacterial intra-species contamination in raw Illumina data. It does this
 by looking for multiple alleles of core, single copy genes.

For complete instructions on installation and usage, please visit [the ConFindr github pages site](https://olc-bioinformatics.github.io/ConFindr/).

## Quickstart

To install ConFindr, use conda (see [here](https://bioconda.github.io/) for instructions on getting bioconda set up): 

`conda install -c bioconda confindr`

To get an example dataset, use this command, which will create a folder called `example-data` in your current working directory: 

`wget https://ndownloader.figshare.com/files/9972709 && tar xf 9972709 && rm 9972709`

As of version `0.7.0` ConFindr can run automatically on _Escherichia_, _Salmonella_, and _Listeria_ with no further 
work on your part - just run:

`confindr.py -i example-data -o example-out`

Once ConFindr finishes running, take a look at the `confindr_report.csv` file found in `example-out` - it shows that multiple
alleles were found for many sites within the genes that ConFindr examines, meaning that this sample is quite contaminated!

If you want to run ConFindr on genera other than the 3 listed above, you'll need to get access to and download rMLST databases. 
[You can find instructions on how to do that here](https://olc-bioinformatics.github.io/ConFindr/install/#downloading-confindr-databases).

### Running ConFindr in a Python Script

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

ConFindr has been published in PeerJ - if you use it in your work, please cite the following:

```
Low AJ, Koziol AG, Manninger PA, Blais B, Carrillo CD. 2019. ConFindr: rapid detection of intraspecies and cross-species contamination in bacterial whole-genome sequence data. PeerJ 7:e6995 https://doi.org/10.7717/peerj.6995
```

