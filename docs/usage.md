## Example Dataset

An example dataset has been uploaded to FigShare. You can download it to your current working directory with the following command:

`wget https://ndownloader.figshare.com/files/9972709 && tar xf 9972709`

This example dataset contains two different serotypes of _Escherichia coli_ mixed together - it's about 80/20 split of O103 and O157. Contamination like this is difficult to detect
with regular tools - it's possible to pick up that it's two different strains, but it can be finicky. ConFindr, however, has no difficulty picking up the fact that this sample is contaminated.

In order to have ConFindr analyze this sample, the parameters you need to provide are:

- `-i, --input_directory`: The path to a directory containing the reads, in FASTQ format, that you want analyzed. If you're using the example dataset, you'll want to enter `example-data`
- `-o, --output_name`: The base name for your output. If you put `output` for this parameter, a folder called `output` will be created, and a file called `confindr_report.csv` with contamination
information will be created in this folder

So, if the `example-data` folder were downloaded to your current working directory and you want to have an output folder called `output`, the command to run ConFindr would be:

`confindr.py -i example-data -o output`

You can use absolute or relative paths, and trailing slashes are also acceptable for the directories specified.
If ConFindr is properly installed, you should see something similar to the following appear on your terminal:

```bash
  2019-04-02 15:06:04  Welcome to ConFindr 0.7.0! Beginning analysis of your samples... 
  2019-04-02 15:06:04  Could not find Escherichia_db_cgderived.fasta 
  2019-04-02 15:06:04  Could not find Listeria_db_cgderived.fasta 
  2019-04-02 15:06:04  Could not find Salmonella_db_cgderived.fasta 
  2019-04-02 15:06:04  Could not find refseq.msh 
  2019-04-02 15:06:04  Databases not present - downloading basic databases now... 
  2019-04-02 15:06:04  Downloading mash refseq sketch... 
  2019-04-02 15:06:07  Downloading cgMLST-derived data for Escherichia, Salmonella, and Listeria... 
  2019-04-02 15:06:12  Did not find rMLST databases, if you want to use ConFindr on genera other than Listeria, Salmonella, and Escherichia, you'll need to download them. Instructions are available at https://olc-bioinformatics.github.io/ConFindr/install/#downloading-confindr-databases
 
  2019-04-02 15:06:12  Beginning analysis of sample example... 
  2019-04-02 15:06:12  Checking for cross-species contamination... 
  2019-04-02 15:06:29  Extracting conserved core genes... 
  2019-04-02 15:06:37  Quality trimming... 
  2019-04-02 15:06:38  Detecting contamination... 
  2019-04-02 15:07:05  Done! Number of contaminating SNVs found: 214
 
  2019-04-02 15:07:05  Contamination detection complete! 
```

The run shouldn't take too long - depending on how powerful your machine is, it should be done in
one to two minutes (slightly longer if an *Escherichia* specific database has not yet been set up).
Once the run is done, you'll be able to inspect your results. Take a look at `output/confindr_report.csv`:
The `ContamStatus` column should read `True`, and the `NumContamSNVs` column should have a value of something close to 200.

In any future uses of ConFindr, databases will not need to be re-downloaded.

## Interpreting ConFindr Results

The results file that ConFindr produces is in comma-separated value (CSV) format, which can be opened by any spreadsheet application (Excel, LibreOffice, etc.) or your favorite text editor.

The file has the following headers: `Sample`, `Genus`, `NumContamSNVs`, `ContamStatus`, `PercentContam`, `PercentContamStandardDeviation`, and `BasesExamined`. Of these, ContamStatus is the most important - it will be `True` if a sample
is contaminated, and `False` if a sample is not contaminated. Detailed descriptions of each header follow.

- `Sample`: The name of the sample. ConFindr will take everything before the first underscore (\_) character to be the name of the sample, as done with samples coming from an Illumina MiSeq.
- `Genus`: The genus that ConFindr thinks your sample is. If ConFindr couldn't figure out what genus your sample is from, this will be NA.
If multiple genera were found, they will all be listed here, separated by a `:`
- `NumContamSNVs`: The number of times ConFindr found sites with more than one base present.
- `ContamStatus`: The most important of all! Will read `True` if contamination is present in the sample, and `False` if contamination is not present. The result will be `True` if any of the following conditions are met:
	- More than 1 contaminating SNV per 10000 base pairs examined was found.
	- There is cross contamination between genera.
- `PercentContam`: Based on the depth of the minor variant for sites with multiple bases, ConFindr guesses
at what percent of your reads come from a contaminant. The more sequencing depth you have, the more accurate this will
get. For lower levels of contamination (around 5 percent) this tends to get overestimated, but the number gets more accurate as
contamination level increases, as well as sequencing depth.
- `PercentContamStandardDeviation`: The standard deviation of the percentage contamination estimate. Very high values may
indicate something strange is going on.
- `BasesExamined`: The number of bases ConFindr examined when making the contamination call. Will usally be around 20kb for rMLST databases,
 and will vary when other databases are used.
- `DatabaseDownloadDate`: Date that rMLST databases were downloaded, if you have them. As these are curated and updated regularly,
it's a good idea to re-run `confindr_database_setup` every now and then.

ConFindr will also produce two CSV files for each sample - one called `samplename_contamination.csv`, which shows the contaminating
sites, and one called `samplename_rmlst.csv`, which shows ConFindr's guess at which allele is present for each rMLST gene.

## Using ConFindr in a Python Script

In the event you'd rather integrate ConFindr into a script than run from the command line, here's how:

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

## Using Schemes other than rMLST

As of ConFindr 0.4.4, ConFindr has the option to use a cgMLST scheme instead of an rMLST scheme for increased
contamination detection sensitivity. This hasn't been tested extensively, but looks to be working. Runtime is
increased by a factor of 2 or 3 compared to running against rMLST genes. 

To use this option, you'll need a a cgMLST FASTA file - all FASTA headers should be in format >genename_allele

In order to decrease computation time, clustering the cgMLST FASTA before running is recommended. CD-HIT on default
parameters does this fairly well.

cgMLST files that are already clustered are available for _Salmonella_ and _Escherichia_. To get them:

Escherichia: `wget 'https://scist01.blob.core.windows.net/olc/Escherichia_cgmlst.fasta?sp=r&st=2020-07-15T13:17:25Z&se=2029-07-15T21:17:25Z&spr=https&sv=2019-10-10&sr=b&sig=zF01kDCgmsmBWJ0cSpyfYi6CLldIjalwU2RgsswKmmI%3D' -O Escherichia_cgmlst.fasta`

Salmonella: `wget 'https://scist01.blob.core.windows.net/olc/Salmonella_cgmlst.fasta?sp=r&st=2020-07-15T13:13:16Z&se=2029-07-15T21:13:16Z&spr=https&sv=2019-10-10&sr=b&sig=ncG%2F5MzKt57p1BUdyFtnhJUk9Yfi6x3rFhSQWPlT2Ek%3D' -O Salmonella_cgmlst.fasta`

When using a cgMLST database, ConFindr will use the provided scheme for all samples regardless of genus.

Actually calling ConFindr with cgMLST:

`confindr.py -i folder-with-Escherichia-files -o cgmlst-output -cgmlst /path/to/Escherichia_cgmlst.fasta`


## Optional Arguments

ConFindr has a few optional arguments that allow you to modify its other parameters. Optional arguments are:

- `-t, --threads`: The number of threads to run ConFindr analysis with. The default is to use all threads available on your machine, and ConFindr scales very well with more threads, so it's recommended that this option be left at the default unless you need the computational resources for something else.
- `-d`, --databases`: Path to ConFindr databases. These will be downloaded automatically if not present.
- `-k`, --keep_files`: Set this flag to keep intermediate files. Useful if you want to do manual inspection of the BAM files
that ConFindr creates, which are deleted by default.
- `-fid, --forward_id`: The identifier for forward reads in your input FASTQ folder. By default, this is `_R1`. If you follow a different naming scheme, this is the parameter to change.
- `-rid, --reverse_id`: The identifier for reverse reads in your input FASTQ folder. By default, this is `_R2`. If you follow a different naming scheme, this is the parameter to change. 
- `-v, --version`: Display ConFindr version and exit.
- `-verbosity, --verbosity`: How much you want printed to the screen. Choose `debug` to get some extra, or `warning` to
get almost nothing. Default is `info`.
- `-b`, `--base_cutoff`: The number of high-quality bases needed to call a site as multiallelic, and therefore 
contributing to contamination. Defaults to 2, which is usually sensitive without producing false positives.
If dealing with high depth samples, adding the `-bf` parameter set to around `0.05` is likely to be helpful in reducing
false positives.
- `-bf`, `--base_fraction_cutoff`: The proportion of high-quality bases needed to call a site as multiallelic, and therefore 
contributing to contamination. Must be between 0 and 1. Not used by default.
- `-q`, `--quality_cutoff`: The phred score a base needs to have before it's considered
trustworthy enough to contribute to a site being multiallelic. Defaults to 20, which should
be suitable for most purposes. 
- `--rmlst`: By default, ConFindr will use custom core-gene derived datasets for _Escherichia_, _Listeria_, and _Salmonella_
instead of rMLST. Activate this flag to force use of rMLST genes for all genera.
- `--cross_details`: By default, when ConFindr finds cross-contaminated samples it stops analysis. Activate
this flag to have analysis of number of cSNVs continue in order to get an estimate of percentage contamination.
