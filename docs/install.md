# Installation

## System Requirements

ConFindr has been tested with Debian-based Linux systems, 
but should in principle work on any flavour of Linux, as well as MacOSX. 
Windows isn't supported, but it may very well be installable via bioconda. If you have any success running
ConFindr on a Windows machine, let me know!

ConFindr should run on any regular desktop/laptop with 8 GB or RAM or more.

## Downloading ConFindr Databases

As of ConFindr 0.7.0, databases for detecting contamination in _Escherichia_, _Listeria_, and _Salmonella_ derived from
core-gene schemes are freely available and will be automatically downloaded by ConFindr when it runs. If you only want 
to run ConFindr on these three genera, nothing further is necessary. If you want to run ConFindr on any other genera, keep
reading on how to get access to the necessary databases.

ConFindr uses the ribosomal multi-locus sequence typing (rMLST) scheme to detect contamination in genera other 
than the ones listed above. These databases are available free for academic use, but you will need to jump through a few 
hoops before you can get access to them due to an associated [licence agreement](https://pubmlst.org/rmlst/rMLST_licence.pdf).
Non-academic use will require a commercial licence.

Here are the steps to getting databases downloaded:

1. Register for a PubMLST account if you do not already have one. Link to register is [here](https://pubmlst.org/bigsdb). 
Click on `Register for a site-wide account.`

2. Login to your account at [https://pubmlst.org/bigsdb](https://pubmlst.org/bigsdb) and request access to `Ribosomal MLST genomes (pubmlst_rmlst_isolates)` and `Ribosomal MLST typing (pubmlst_rmlst_seqdef)` under 'Database registrations'. Additionally, create a PubMLST API key on the same page under 'API keys'. The generated client ID (consumer key) and client secret (consumer secret) will enable you to access the database programatically.

3. Once you've gotten your consumer key and consumer secret, put them into a text file
with the key on the first line and the secret on the second. It should look something like the below
snippet:

```
efKXmqp2D0EBlMBkZaGC2lPf
F$M+fQ2AFFB2YBDfF9fpHF^qSWJdmmN%L4Fxf5Gur3
```

4. Install ConFindr as shown in the next section.

5. With ConFindr installed, use the command `confindr_database_setup` to have ConFindr download the latest version
of the rMLST databases. This script takes two arguments - a `-s` where you give the path to the text file containing your consumer 
key and secret, and a `-o` to specify where you want the sequences downloaded. Only the `-s` is mandatory. If your output
directory is not specified, ConFindr will first search for an environmental variable called `CONFINDR_DB`, and if it can't
find that it will automatically download to a folder called `.confindr_db` in your home directory.

## Installing Using Conda (Recommended)

1. Follow the instructions [here](https://bioconda.github.io/) to add the Bioconda channel to your list of conda channels, if it hasn't already been added.

2. Install ConFindr into a new conda environment named 'confindr':

`conda create -n confindr -c bioconda confindr=0.8.2`

3. Activate the new conda environment:

`conda activate confindr`

Typing `confindr -h` into the command-line will show the help menu for the program. See the [Usage](usage.md) section for instructions on how to use ConFindr, including a ConFindr run on an example dataset.

## Manual Installation

### Executable

ConFindr can also be installed using `pip`. Use of a virtual environment for ConFindr is highly recommended. To create a virtualenv:

1. Create an empty directory (i.e. `mkdir ~/Virtual_Environments/ConFindr`).
2. Virtualenv that directory (`virtualenv -p /usr/bin/python3 ~/Virtual_Environments/ConFindr`).
3. Activate the virtualenv (`source ~/Virtual_Environments/ConFindr/bin/activate`).
4. Install ConFindr—this should also install any packages that ConFindr depends upon (`pip install confindr`).

With this done, you'll need to make sure that any necessary dependencies are installed.

### Dependencies

Before using ConFindr when installed using `pip`, you'll need to download and add the following programs to your $PATH:

- [BBMap (>=39.01)](https://jgi.doe.gov/data-and-tools/bbtools/)
- [Mash (>=2.3)](https://github.com/marbl/Mash/releases)
- [KMA (>=1.4.9)](https://bitbucket.org/genomicepidemiology/kma)
- [Python (>=3.9.15)](https://www.python.org/downloads/)
- [SAMtools (>=1.17)](https://github.com/samtools/samtools)
- [pysam (>=0.21.0)](https://pypi.org/project/pysam/)

If you want to run ConFindr in Nanopore mode (`-dt Nanopore`), you'll also need to install [minimap2](https://github.com/lh3/minimap2).

Instructions on adding programs to your $PATH can be found [here](https://stackoverflow.com/questions/14637979/how-to-permanently-set-path-on-linux-unix).

If ConFindr can't find these dependencies when you try to run it, you will see an error message and the program will quit.



