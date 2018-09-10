[![Build status](https://travis-ci.org/lowandrew/ConFindr.svg?master)](https://travis-ci.org/lowandrew)
[![PyPI version](https://badge.fury.io/py/confindr.svg)](https://badge.fury.io/py/confindr)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/confindr/README.html)

# ConFindr

This program is designed to find bacterial intra-species contamination in raw Illumina data. It does this
 by looking for multiple alleles of rMLST genes, which are known to be universal across the bacterial kingdom
 and present only in single copies.

For complete instructions on installation and usage, please visit [the ConFindr github pages site](https://lowandrew.github.io/ConFindr/).

## Quickstart

To install ConFindr, use conda: `conda install -c bioconda confindr`

To get an example dataset, use this command: `wget https://ndownloader.figshare.com/files/9972709 && tar xf 9972709 && rm 9972709`

To run confindr on that dataset: `confindr.py -i example-data -o example-out`

This will download the databases ConFindr needs and run ConFindr on the example dataset (a intraspecies _Escherichia coli_ mix), and put the results in a folder
called `example-out` in your current working directory. Take a look at `confindr_report.csv` in that directory to see
the ConFindr results, which will show the sample is contaminated.