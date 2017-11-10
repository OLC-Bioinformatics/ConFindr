# System Requirements

ConFindr has a fair number of dependencies. The easiest way to install ConFindr is by using docker, but other methods are also possible. 
Testing of ConFindr has been done with Ubuntu 16.04 and Linux Mint. Other variants of linux should have no issues,
and MacOS systems should also work. Windows is not supported at this time, but Windows users should be able to use the docker version.

## Installing Using Docker

To install using docker: `docker pull olcbioinformatics/confindr`

This should download the newest image, complete with the databases that ConFindr needs. For usage with docker, see the Usage page (HYPERLINK ME!)

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

If ConFindr can't find these dependencies when you try to run it, an _ImportError_ will be raised with a list of unfindable dependences.


#### Databases

Once you have the executable and dependencies installed, you'll just need to download the databases that ConFindr depends on.

To do this, you'll need to install Git LFS (instructions [here](https://git-lfs.github.com/)). 

Then, clone the ConFindr Git repository (`git clone https://github.com/lowandrew/ConFindr.git`). The `databases` folder is the important one - you'll need it for calling ConFindr, as seen in the `Usage` section. 
