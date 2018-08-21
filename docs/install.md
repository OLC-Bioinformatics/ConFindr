# System Requirements

ConFindr has been tested with Debian-based Linux systems, but should in principle work on any flavour of Linux, as well as MacOSX. Windows is not supported at this time,
but Windows users may be able to use the ConFindr Docker image in order to run ConFindr.

To run ConFindr, your computer should have a minimum of 8GB of RAM, and at least 10GB of disk space. Any number of processors will work, with more generally being better.

## Downloading ConFindr Databases

#### Automatic Download

As of ConFindr 0.3.4, ConFindr databases will automatically downloaded (by default to ~/.confindr_db - this location can be changed
by setting the environmental variable $CONFINDR_DB). You may still manually download the databases and supply the path
as done for previous versions.

#### Manual Download
The databases necessary for making ConFindr run are available for download from FigShare.

Navigate to the place you would like to download the database, and use the following commands to download and uncompress the folder:

`wget https://ndownloader.figshare.com/files/11864267 && tar xf 11864267 && rm 11864267`

These commands should create a folder called `databases` in your current working directory. This folder contains everything you need to run ConFindr - it's what will be specified with the `-d` option.

## Installing Using Conda (Recommended)

ConFindr is available within bioconda - to get bioconda installed and running see instructions [here](https://bioconda.github.io/).

With bioconda running, you can install ConFindr with the following command:

`conda install -c bioconda confindr`

With that done, typing `confindr.py` will bring access the ConFindr pipeline. See the [Usage](usage.md) section for instructions on how to use ConFindr, including a ConFindr run on an example dataset.

## Installing Using Docker

ConFindr can also be used via Docker, which will take care of installing all of ConFindr's dependencies for you. If you do not have Docker installed, instructions on how to install it
can be found [here](https://docs.docker.com/engine/installation/).

With Docker installed, all you have to do to install ConFindr is enter the following command, which will pull the ConFindr image from the Docker Hub and put it on your machine:

`docker pull olcbioinformatics/confindr`

You can verify that the pull was successful by entering the command `docker images`. You should see `olcbioinformatics/confindr` in the list. For instructions on using the image, see 
the [Usage](usage.md) section.

## Manual Install

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

If ConFindr can't find these dependencies when you try to run it, you will see an error message and the program will quit.



