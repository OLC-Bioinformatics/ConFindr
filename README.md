[![Build status](https://travis-ci.org/lowandrew/ConFindr.svg?master)](https://travis-ci.org/lowandrew)
[![PyPI version](https://badge.fury.io/py/confindr.svg)](https://badge.fury.io/py/confindr)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/confindr/README.html)

# ConFindr

This program is designed to find bacterial intra-species contamination in raw Illumina data. It does this
 by looking for multiple alleles of rMLST genes, which are known to be universal across the bacterial kingdom
 and present only in single copies.

For complete instructions on installation and usage, please visit [the ConFindr github pages site](https://lowandrew.github.io/ConFindr/).

## Quickstart

To install ConFindr, use conda: 

`conda install -c bioconda confindr`

To get an example dataset, use this command: 

`wget https://ndownloader.figshare.com/files/9972709 && tar xf 9972709 && rm 9972709`

To run confindr on that dataset: 

`confindr.py -i example-data -o example-out`

This will download the databases ConFindr needs and run ConFindr on the example dataset (a intraspecies _Escherichia coli_ mix), and put the results in a folder
called `example-out` in your current working directory. Take a look at `confindr_report.csv` in that directory to see
the ConFindr results, which will show the sample is contaminated.


### Running ConFindr in a Python Script

If you want to run ConFindr from within a script instead of running from the command line, here's how:

```python
from confindr import confindr

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

### Using a cgMLST scheme instead of rMLST

As of ConFindr 0.4.4, ConFindr has the option to use a cgMLST scheme instead of an rMLST scheme for increased
contamination detection sensitivity. This hasn't been tested extensively, but looks to be working. Runtime is
increased by a factor of 2 or 3 compared to running against rMLST genes.

To use this option, you'll need a a cgMLST FASTA file - all FASTA headers should be in format >genename_allele

In order to decrease computation time, clustering the cgMLST FASTA before running is recommended. CD-HIT on default
parameters does this fairly well.

cgMLST files that are already clustered are available for _Salmonella_ and _Escherichia_. To get them:

Escherichia: `wget 'https://scist01.blob.core.windows.net/olc/Escherichia_cgmlst.fasta?sp=r&st=2018-10-25T13:45:13Z&se=2020-10-31T21:45:13Z&spr=https&sv=2017-11-09&sig=0fgvcf6R%2BSSiz7gPm7KKm5M78wpjAPjGawhzt%2BlY7iE%3D&sr=b' -O Escherichia_cgmlst.fasta`

Salmonella: `wget 'https://scist01.blob.core.windows.net/olc/Salmonella_cgmlst.fasta?sp=r&st=2018-10-25T13:35:21Z&se=2020-10-31T21:35:21Z&spr=https&sv=2017-11-09&sig=w5Pq9e4hsa6PGr458%2Bx4b5zqf4F2a6OUUL9H3ewQTNc%3D&sr=b' -O Salmonella_cgmlst.fasta`

When using a cgMLST database, ConFindr will use the provided scheme for all samples regardless of genus.

Actually calling ConFindr with cgMLST:

`confindr.py -i folder-with-Escherichia-files -o cgmlst-output -cgmlst /path/to/Escherichia_cgmlst.fasta`


## Reporting Issues

If you have any problems installing or running ConFindr, or have feature request,
please open an issue here on GitHub.