# Usage with a Pip Install

If you used pip to install, all you need to do is type `confindr.py` on the command line. Doing so without providing parameters will give a message telling you what parameters you do need to provide.

The parameters you need to provide are (in this order!):

- `fastq_directory`: The path to a directory containing the reads, in FASTQ format, that you want analyzed.
- `output_directory`: The base name for your output. For example, putting `/home/user/confindr_output` would create an output file called `/home/user/confindr_output.csv`
- `databases`: The path to the databases directory obtained when cloning the ConFindr GitHub repository.

So, if the ConFindr repository was cloned to `/home/user` and the data to be analyzed is in `/home/user/example-data`, the command to run ConFindr would be:

`confindr.py /home/user/example-data /home/user/confindr_output /home/user/ConFindr/databases/`

This command will create a file called `/home/user/confindr_output.csv` with results, and a file called `/home/user/confindr_output.log` which shows the log of what STDOUT and STDERR from programs called
in the ConFindr pipeline.

# Usage with Docker

Usage with docker is very similar to usage with a pip install. In order to have the output from the docker container accessible by you, you'll need to mount some folders - it's recommended to put the output in your input folder in this case. In this case, if the data to be analyzed is in `/home/user/example-data`, the command to run ConFindr would be:

`docker run -v /home/user/example-data/:/home/user/example-data/ olcbioinformatics/confindr confindr.py /home/user/example-data/ /home/user/example-data/confindr_output /home/user/ConFindr/databases/` 

In this case, the output file `confindr_output.csv` will be found in `/home/user/example-data/`.

## Optional Arguments

ConFindr has a few optional arguments that allow you to modify its parameters. These should be placed after the positional arguments described above. Optional arguments are:

- `-t, --threads`: The number of threads to run ConFindr analysis with. The default is to use all threads available on your machine, and ConFindr scales very well with more threads, so it's recommended that this option be left at the default unless you need the computational resources for somethin else.
- `-n, --number_subsamples`: The number of times you want ConFindr to sample your rMLST reads to try to detect contamination. By default this is set to 5.
- `-k, --kmer-size`: The kmer size ConFindr uses to try to detect contamination. Default is 31. Usage with other values may produce very unreliable results and is _not_ recommended.
- `-s, --subsample_depth`: Coverage depth to subsample to. Default value is 20, which provides a good tradeoff between sensitvity and specificity. Going any lower will make it very difficult to detect contamination, and going higher will increase the false positive rate.
- `-c, --kmer_cutoff`: The cutoff for the number of times a kmer must be seen before it is considered trustworthy and is included in the analysis. By default set to 2. Setting any lower this this essentially guarantees that your analysis will be overrun by false positives called by sequencing errors.

Generally speaking, none of these parameters should be changed; ConFindr has been tested extensively with its default parameters and been found to work very well. 

