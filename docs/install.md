# System Requirements

ConFindr has a fair number of dependencies.
Testing of ConFindr has been done with Ubuntu 16.04 and Linux Mint. Other variants of linux should have no issues,
and MacOS systems should also work. Windows is not supported at this time, but Windows users should be able to use the docker version.


## Installing Using Pip

#### Executable

ConFindr can also be installed using pip. Use of a virtual environment for ConFindr is highly recommended. To create a virtualenv:

- Create an empty directory (i.e. `mkdir ~/Virtual_Environments/ConFindr`)
- Virtualenv that directory (`virtualenv -p /usr/bin/python3 ~/Virtual_Environments/ConFindr`)
- Activate the virtualenv (`source ~/Virtual_Environments/ConFindr/bin/activate`)
- Install ConFindr - this should also install any packages that ConFindr depends on (`pip install confindr`)

With this done, you'll need to make sure that any necessary dependencies are installed.

#### Dependencies

Before using ConFindr, you'll need to download and add the following programs to your $PATH:

- [BBTools (>=37.23)](https://jgi.doe.gov/data-and-tools/bbtools/)
- [Jellyfish (>= 2.2.6)](https://github.com/gmarcais/Jellyfish/releases)
- [NCBI BLAST+ (>=2.2.31](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [Mash (>=2.0)](https://github.com/marbl/Mash/releases)
- [Python (>=3.5)](https://www.python.org/downloads/)

Instructions on adding programs to your $PATH can be found [here](https://stackoverflow.com/questions/14637979/how-to-permanently-set-path-on-linux-unix).

If ConFindr can't find these dependencies when you try to run it, you will see an error message. ConFindr will continue to attempt to run, but will likely crash at some point in the process.


#### Databases

The databases necessary for making ConFindr run are available for download from FigShare.

Navigate to the place you would like to download the database, and use the following commands to download and uncompress the folder:

`wget https://ndownloader.figshare.com/files/9827251
tar xf 9827251`

These commands should create a folder called `databases` in your current working directory. This folder contains everything you need to run ConFindr - it's what will be specified with the `-d` option.
