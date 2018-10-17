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

Important note: Under current default settings, ConFindr will call a lot of false positives if you 
input high depth (>100X) samples. If your samples are high depth, add `-bf 0.05` to your command - this will make it so that
at least 5 percent of bases in a column have to support a call, as well as having at least 2 high quality bases.
This will likely be changed to a default setting in a future version of ConFindr.

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

## Reporting Issues

If you have any problems installing or running ConFindr, or have feature request,
please open an issue here on GitHub.