# Usage

## Example Dataset

An [example dataset](https://figshare.com/articles/dataset/Minimal_dataset_for_ConFindr_testing_using_pytest/22852937) has been uploaded to FigShare for testing ConFindr.
This example dataset contains FASTQ reads from contaminated and uncontaminated genomes of _Escherichia coli_ O121:H19, _Salmonella enterica_ Heidelberg, and _Listeria monocytogenes_.
These reads were obtained by selecting samples from the originally published dataset that was used to test ConFindr's performance, and then downsampled by different factors to obtain a more minimal dataset.
The metadata for these samples can be downloaded from [here](https://figshare.com/ndownloader/files/40599716), with additional metadata available in the Supplemental Information of the [original ConFindr publication](https://peerj.com/articles/6995).
The dataset may take a few minutes to download depending on your internet connection speed.

If you wish to download the dataset to analyze it yourself using `confindr.py`, you can run the following command to download and extract the files within your current working directory:

```bash
wget https://figshare.com/ndownloader/files/41228577 -O test_samples.tar.gz && \
    tar -xzvf test_samples.tar.gz && \
    rm test_samples.tar.gz
```

ConFindr can also be tested with [pytest](https://docs.pytest.org/en/7.3.x/) automated testing.
This is the recommended testing strategy for those who wish to contribute to ConFindr on GitHub.
We recommend installing ConFindr within a virtual environment using Conda (see instructions [here](install.md#installing-using-conda-recommended)), cloning the ConFindr GitHub repository, and then installing ConFindr in development mode using the following commands:

```bash
conda activate confindr
git clone https://github.com/OLC-Bioinformatics/ConFindr
cd ConFindr/
pip install -e .
wget https://figshare.com/ndownloader/files/41228577 -O tests/test_samples.tar.gz && \
  tar -xzvf tests/test_samples.tar.gz -C tests/ && \
  rm tests/test_samples.tar.gz
```

To run the automated tests with pytest:

```bash
conda activate confindr
cd ConFindr/
python -m pytest tests/ -vvv
```

You should see output similar to below if all tests have passed:

```
================================================ test session starts =================================================
platform linux -- Python 3.7.12, pytest-4.5.0, py-1.11.0, pluggy-0.11.0 -- /home/liam/miniconda3/envs/confindr-8.0/bin/python
cachedir: .pytest_cache
rootdir: /media/liam/Cerulean/projects/PRJ-LPB-00001/bin/confindr-8.0/ConFindr
collected 38 items                                                                                                   

tests/test_confindr.py::test_integration PASSED                                                                [  2%]
tests/test_confindr.py::test_present_dependency PASSED                                                         [  5%]
tests/test_confindr.py::test_nonexistent_dependency PASSED                                                     [  7%]
tests/test_confindr.py::test_r1_fastqs PASSED                                                                  [ 10%]
tests/test_confindr.py::test_1_fastqs PASSED                                                                   [ 13%]
tests/test_confindr.py::test_empty_fastqs PASSED                                                               [ 15%]
tests/test_confindr.py::test_unpaired_fastq PASSED                                                             [ 18%]
tests/test_confindr.py::test_run_cmd_success PASSED                                                            [ 21%]
tests/test_confindr.py::test_run_cmd_failure_exit_code PASSED                                                  [ 23%]
tests/test_confindr.py::test_two_hq_bases_above_threshold PASSED                                               [ 26%]
tests/test_confindr.py::test_just_one_hq_bases_above_threshold PASSED                                          [ 28%]
tests/test_confindr.py::test_two_hq_bases_above_threshold_custom_params PASSED                                 [ 31%]
tests/test_confindr.py::test_just_one_hq_base_above_threshold_custom_params PASSED                             [ 34%]
tests/test_confindr.py::test_three_hq_bases_above_threshold PASSED                                             [ 36%]
tests/test_confindr.py::test_two_out_of_three_hq_bases_above_threshold PASSED                                  [ 39%]
tests/test_confindr.py::test_two_hq_bases_above_fraction_threshold PASSED                                      [ 42%]
tests/test_confindr.py::test_two_hq_bases_above_fraction_threshold_low_coverage PASSED                         [ 44%]
tests/test_confindr.py::test_two_hq_bases_above_fraction_threshold_low_coverage_one_base_counts PASSED         [ 47%]
tests/test_confindr.py::test_just_one_hq_bases_above_fraction_threshold PASSED                                 [ 50%]
tests/test_confindr.py::test_three_hq_bases_above_fraction_threshold PASSED                                    [ 52%]
tests/test_confindr.py::test_two_out_of_three_hq_bases_above_fraction_threshold PASSED                         [ 55%]
tests/test_confindr.py::test_valid_base_fraction_none PASSED                                                   [ 57%]
tests/test_confindr.py::test_valid_base_fraction_zero PASSED                                                   [ 60%]
tests/test_confindr.py::test_valid_base_fraction_one PASSED                                                    [ 63%]
tests/test_confindr.py::test_valid_base_fraction_between_zero_one PASSED                                       [ 65%]
tests/test_confindr.py::test_invalid_base_fraction PASSED                                                      [ 68%]
tests/test_confindr.py::test_total_length_fasta PASSED                                                         [ 71%]
tests/test_confindr.py::test_write_output_creates_file_if_does_not_exist PASSED                                [ 73%]
tests/test_confindr.py::test_write_output_appends_if_file_does_exist PASSED                                    [ 76%]
tests/test_confindr.py::test_base_dict_to_string_two_base_descending PASSED                                    [ 78%]
tests/test_confindr.py::test_base_dict_to_string_two_base_ascending PASSED                                     [ 81%]
tests/test_confindr.py::test_base_dict_to_string_three_bases PASSED                                            [ 84%]
tests/test_confindr.py::test_valid_xmx_string_gigabytes PASSED                                                 [ 86%]
tests/test_confindr.py::test_valid_xmx_string_megabytes PASSED                                                 [ 89%]
tests/test_confindr.py::test_valid_xmx_string_kilobytes PASSED                                                 [ 92%]
tests/test_confindr.py::test_invalid_xmx_bad_suffix PASSED                                                     [ 94%]
tests/test_confindr.py::test_invalid_xmx_float PASSED                                                          [ 97%]
tests/test_confindr.py::test_invalid_xmx_not_an_integer PASSED                                                 [100%]

============================================= 38 passed in 57.70 seconds =============================================
```

## Running ConFindr using included core-gene databases

**ConFindr has only been validated using rMLST databases. Please use them if possible (`--rmlst`, see instructions below).**

In order to analyze samples using ConFindr with the pre-included core-gene databases for _Escherichia_, _Listeria_, and _Salmonella_, the parameters you need to provide are:

- `-i, --input_directory`: The path to a directory containing the reads, in FASTQ format, that you want analyzed. If you're using the example dataset, you'll want to enter `example-data`
- `-o, --output_name`: The base name for your output. If you put `output` for this parameter, a folder called `output` will be created, and a file called `confindr_report.csv` with contamination information will be created in this folder

So, if the `test_samples` folder were downloaded to your current working directory and you want to have an output folder called `output`, the command to run ConFindr would be:

`confindr -i test_samples -o output`

You can use absolute or relative paths, and trailing slashes are also acceptable for the directories specified.
If ConFindr is properly installed, you should see something similar to the following appear on your terminal:

```bash
$ confindr -i test_samples -o output
  2023-05-16 15:25:51  Welcome to ConFindr 0.8.1! Beginning analysis of your samples... 
  2023-05-16 15:25:51  Could not find Escherichia_db_cgderived.fasta 
  2023-05-16 15:25:51  Could not find Listeria_db_cgderived.fasta 
  2023-05-16 15:25:51  Could not find Salmonella_db_cgderived.fasta 
  2023-05-16 15:25:51  Could not find refseq.msh 
  2023-05-16 15:25:51  Databases not present - downloading basic databases now... 
  2023-05-16 15:25:51  Downloading mash refseq sketch... 
  2023-05-16 15:25:53  Downloading cgMLST-derived data for Escherichia, Salmonella, and Listeria... 
  2023-05-16 15:25:55  Since this is the first time you are using this database, it needs to be indexed by KMA. This might take a while 
  2023-05-16 15:25:57  Since this is the first time you are using this database, it needs to be indexed by KMA. This might take a while 
  2023-05-16 15:25:58  Since this is the first time you are using this database, it needs to be indexed by KMA. This might take a while 
  2023-05-16 15:25:59  Did not find rMLST databases, if you want to use ConFindr on genera other than Listeria, Salmonella, and Escherichia, you'll need to download them. Instructions are available at https://olc-bioinformatics.github.io/ConFindr/install/#downloading-confindr-databases
 
  2023-05-16 15:25:59  Beginning analysis of sample SRX5084910_SRR8268082... 
  2023-05-16 15:25:59  Checking for cross-species contamination... 
  2023-05-16 15:26:00  Extracting conserved core genes... 
  2023-05-16 15:26:01  Quality trimming... 
  2023-05-16 15:26:01  Detecting contamination... 
  2023-05-16 15:26:05  Done! Number of contaminating SNVs found: 208

...

  2023-05-16 15:26:51  Contamination detection complete! 
```

The run shouldn't take too long—depending on how powerful your machine is, it should be done in one to two minutes.
Once the run is complete, you'll be able to inspect your results.
Take a look at `output/confindr_report.csv`.

If your files are paired-end reads that don't following the `_R1` for forward reads and `_R2` for reverse reads naming convention, then you'll need to specify the pattern to search for using `-fid` and `-rid`.

### Running ConFindr using rMLST databases

The rMLST databases must be downloaded and set up prior to use. You can find the instructions for this [here](install.md#downloading-confindr-databases).

To analyze the test samples using ConFindr with the downloaded rMLST alleles databases, you need to provide an `--input_directory` and an  `--output_name`, as well as the `--rmlst` option on the command line. The command to analyze the test sample in this case would be:

`confindr -i test_samples -o output --rmlst`

In any future uses of ConFindr, databases will not need to be re-downloaded.

## Interpreting ConFindr Results

The results file that ConFindr produces is in comma-separated value (CSV) format, which can be opened by any spreadsheet application (Excel, LibreOffice, etc.) or your favorite text editor.

The file has the following headers: `Sample`, `Genus`, `NumContamSNVs`, `ContamStatus`, and `BasesExamined`. Of these, **ContamStatus is the most important**—it will be `True` if a sample
is contaminated, and `False` if a sample is not contaminated. Detailed descriptions of each header follow.

- `Sample`: The name of the sample. ConFindr will take everything before the first underscore (\_) character to be the name of the sample, as done with samples coming from an Illumina MiSeq.
- `Genus`: The genus that ConFindr thinks your sample is. If ConFindr couldn't figure out what genus your sample is from, this will be NA.
If multiple genera were found, they will all be listed here, separated by a `:`
- `NumContamSNVs`: The number of times ConFindr found sites with more than one base present.
- `ContamStatus`: The most important of all! Will read `True` if contamination is present in the sample, and `False` if contamination is not present. The result will be `True` if any of the following conditions are met:
	- More than 1 contaminating SNV per 10000 base pairs examined was found.
	- There is cross contamination between genera.
- `BasesExamined`: The number of bases ConFindr examined when making the contamination call. Will usally be around 20kb for rMLST databases,
 and will vary when other databases are used.
- `DatabaseDownloadDate`: Date that rMLST databases were downloaded, if you have them. As these are curated and updated regularly,
it's a good idea to re-run `confindr_database_setup` every now and then.

ConFindr will also produce two CSV files for each sample: one called `samplename_contamination.csv`, which shows the contaminating
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

To use this option, you'll need a a cgMLST FASTA file—all FASTA headers should be in format >genename_allele

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
- `-d`, `--databases`: Path to ConFindr databases. These will be downloaded automatically if not present.
- `-k`, `--keep_files`: Set this flag to keep intermediate files. Useful if you want to do manual inspection of the BAM files
that ConFindr creates, which are deleted by default.
- `-fid, --forward_id`: The identifier for forward reads in your input FASTQ folder. By default, this is `_R1`. If you follow a different naming scheme, this is the parameter to change.
- `-rid, --reverse_id`: The identifier for reverse reads in your input FASTQ folder. By default, this is `_R2`. If you follow a different naming scheme, this is the parameter to change. 
- `-v, --version`: Display ConFindr version and exit.
- `-verbosity, --verbosity`: How much you want printed to the screen. Choose `debug` to get some extra, or `warning` to
get almost nothing. Default is `info`.
- `-b`, `--base_cutoff`: The number of high-quality bases needed to call a site as multiallelic, and therefore 
contributing to contamination. Defaults to 3, which is usually sensitive without producing false positives.
If dealing with high depth samples, adding the `-bf` parameter set to around `0.05` is likely to be helpful in reducing
false positives.
- `-bf`, `--base_fraction_cutoff`: The proportion of high-quality bases needed to call a site as multiallelic, and therefore 
contributing to contamination. Must be between 0 and 1. Not used by default.
- `-e`, `--error_cutoff`: Value to use for the calculated error cutoff when setting the `--base_cutoff` value. Default is 1.0%.
- `-q`, `--quality_cutoff`: The phred score a base needs to have before it's considered
trustworthy enough to contribute to a site being multiallelic. Defaults to 20, which should
be suitable for most purposes. 
- `--rmlst`: By default, ConFindr will use custom core-gene derived datasets for _Escherichia_, _Listeria_, and _Salmonella_
instead of rMLST. Activate this flag to force use of rMLST genes for all genera.
- `-tmp`, `--tmp`: If your ConFindr databases are in a location that you don't have write access to, you can enter this option to specify a temporary directory to save genus-specific databases to.
- `-Xmx`, `--Xmx`: Very occasionally, parts of the pipeline that use the BBMap suite will have their memory reservation fail and request insufficient, or sometimes negative, memory. If this happens to you, you can use this flag to override automatic memory reservation and use an amount of memory requested by you. `-Xmx 20g` will specify 20 gigabytes of RAM, and `-Xmx 800m` will specify 800 megabytes.
- `-cgmlst`, `--cgmlst`: Path to a cgMLST database to use for contamination detection instead of using the default rMLST database. Sequences in this file should have headers in format `>genename_allelenumber`. To speed up ConFindr runs, clustering the cgMLST database with CD-HIT before running ConFindr is recommended. This is highly experimental, and results should be interpreted with great care.
- `--fasta`: If activated, ConFindr will look for FASTA files instead of FASTQ for unpaired reads.
- `-m`, `--min_matching_hashes`: Minimum number of matching hashes in a MASH screen in order for a genus to be considered present in a sample. Default is 150.