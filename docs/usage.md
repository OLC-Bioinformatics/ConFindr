# Usage with a Pip/Conda Install

If you used pip or conda to install, all you need to do is activate your ConFindr virtualenv/type `confindr.py` on the command line. Doing so without providing parameters will give a message telling you what parameters you do need to provide.


## Example Dataset

An example dataset has been uploaded to FigShare. You can download it to your current working directory with the following command:

`wget https://ndownloader.figshare.com/files/9972709 && tar xf 9972709`

This example dataset contains two different serotypes of _Escherichia coli_ mixed together - it's about 80/20 split of O103 and O157. Contamination like this is difficult to detect
with regular tools - it's possible to pick up that it's two different strains, but it can be finicky. ConFindr, however, has no difficulty picking up the fact that this sample is contaminated.

In order to have ConFindr analyze this sample, the parameters you need to provide are:

- `-i, --input_directory`: The path to a directory containing the reads, in FASTQ format, that you want analyzed. If you're using the example dataset, you'll want to enter `example-data`
- `-o, --output_name`: The base name for your output. If you put `output` for this parameter, a folder called `output` will be created, and a file called `confindr_report.csv` with contamination
information will be created in this folder
- `-d, ---databases`: The path to the databases directory obtained when downloading the `databases` folder from FigShare.

So, if the `databases` and `example-data` folders were downloaded to your current working directory and you want to have an output folder called `output`, the command to run ConFindr would be:

`confindr.py -i example-data -o output -d databases`

You can use absolute or relative paths, and trailing slashes are also acceptable for the directories specified. If ConFindr is properly installed, you should see something similar to the following appear on your terminal:

```bash
  2018-06-07 12:41:11  Welcome to ConFindr 0.3.1! Beginning analysis of your samples...
  2018-06-07 12:41:11  Beginning analysis of sample example...
  2018-06-07 12:41:11  Checking for cross-species contamination...
  2018-06-07 12:41:23  Extracting rMLST genes...
  2018-06-07 12:41:27  Quality trimming...
  2018-06-07 12:41:27  Beginning 3 cycles of contamination detection...
  2018-06-07 12:41:27  Working on cycle 1 of 3...
  2018-06-07 12:41:35  Working on cycle 2 of 3...
  2018-06-07 12:41:43  Working on cycle 3 of 3...
  2018-06-07 12:41:51  Finished analysis of sample example!
  2018-06-07 12:41:51  Contamination detection complete!
```

The run shouldn't take too long - depending on how powerful your machine is, it should be done in
one to two minutes (slightly longer if an *Escherichia* specific database has not yet been set up.
Once the run is done, you'll be able to inspect your results.
The `ContamStatus` column should read `True`, and the `NumContamSNVs` column should have a value somewhere between 30 and 40.

# Usage with a Docker Install

If you used Docker to install ConFindr, usage will be slightly different. Assuming you have the `databases` from the [Installation](install.md) step and `example-data` from above in your current working
directory, the command to run ConFindr would be:

`docker run -it -v /path/to/current/directory/:/data olcbioinformatics/confindr confindr.py -i /data/example-data -o /data/output -d /data/databases`

You should see the same output to the terminal that was mentioned above, and have the same output files in a folder called `output` in your current working directory.

## Interpreting ConFindr Results

The results file that ConFindr produces is in comma-separated value (CSV) format, which can be opened by any spreadsheet application (Excel, LibreOffice, etc.) or your favorite text editor.

The file has the following headers: Sample, Genus, NumContamSNVs, NumUniqueKmers, CrossContamination, and ContamStatus. Of these, ContamStatus is the most important - it will be `True` if a sample
is contaminated, and `False` if a sample is not contaminated. Detailed descriptions of each header follow.

- `Sample`: The name of the sample. ConFindr will take everything before the first underscore (\_) character to be the name of the sample, as done with samples coming from an Illumina MiSeq.
- `Genus`: The genus that ConFindr thinks your sample is. If ConFindr couldn't figure out what genus your sample is from, this will be NA.
If multiple genera were found, they will all be listed here, separated by a `:`
- `NumContamSNVs`: The number of times ConFindr found a kmer that had a mismatch, indicating the potential for multiple alleles of one gene being present. Completely clean samples should have a value of 0.
- `NumUniqueKmers`: The number of unique kmers found by ConFindr for a sample. Numbers substantially above the total length of the rMLST genes (~35000 base pairs) can indicate contamination.
- `ContamStatus`: The most important of all! Will read `True` if contamination is present in the sample, and `False` if contamination is not present. The result will be `True` if any of the following conditions are met:
	- 3 or more contaminating SNVs are found. 
	- More than 45000 unique rMLST kmers are found.
	- There is cross contamination between genera.


## Optional Arguments

ConFindr has a few optional arguments that allow you to modify its other parameters. Optional arguments are:

- `-t, --threads`: The number of threads to run ConFindr analysis with. The default is to use all threads available on your machine, and ConFindr scales very well with more threads, so it's recommended that this option be left at the default unless you need the computational resources for something else.
- `-n, --number_subsamples`: The number of times you want ConFindr to sample your rMLST reads to try to detect contamination. By default this is set to 3.
- `-k, --kmer-size`: The kmer size ConFindr uses to try to detect contamination. Default is 31. Usage with other values may produce very unreliable results and is _not_ recommended.
- `-s, --subsample_depth`: Coverage depth to subsample to. Default value is 20, which provides a good tradeoff between sensitvity and specificity. Going any lower will make it very difficult to detect contamination, and going higher will increase the false positive rate.
- `-c, --kmer_cutoff`: The cutoff for the number of times a kmer must be seen before it is considered trustworthy and is included in the analysis. By default set to 2. Setting any lower this this essentially guarantees that your analysis will be overrun by false positives called by sequencing errors.
- `-fid, --forward_id`: The identifier for forward reads in your input FASTQ folder. By default, this is `_R1`. If you follow a different naming scheme, this is the parameter to change.
- `-rid, --reverse_id`: The identifier for reverse reads in your input FASTQ folder. By default, this is `_R2`. If you follow a different naming scheme, this is the parameter to change. 
- `-v, --version`: Display ConFindr version and exit.

Generally speaking, none of these parameters should be changed; ConFindr has been tested extensively with its default parameters and been found to work very well. 

