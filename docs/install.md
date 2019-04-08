# System Requirements

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
than the ones listed above. These databases are 
freely available, but you will need to jump through a few hoops before you can get access to them due to an 
associated [licence agreement](https://pubmlst.org/rmlst/rMLST_licence.pdf).

Here are the steps to getting databases downloaded:

- Register for a PubMLST account if you do not already have one. Link to register is [here](https://pubmlst.org/bigsdb). 
Click on `Register for a site-wide account.`

- Login to your account at [https://pubmlst.org/bigsdb](https://pubmlst.org/bigsdb) and request access to 
Ribosomal MLST genome and Ribosomal MLST locus/sequence definitions under `Registrations`. Additionally, email Keith Jolley
(keith.jolley@zoo.ox.ac.uk) and request a consumer key and consumer secret so that you'll be able
to access the database programatically.

- Once you've gotten your consumer key and consumer secret from Keith, put them into a text file
with the key on the first line and the secret on the second. It should look something like the below
snippet:

```
efKXmqp2D0EBlMBkZaGC2lPf
F$M+fQ2AFFB2YBDfF9fpHF^qSWJdmmN%L4Fxf5Gur3
```

- Install ConFindr as shown in the next section.

- With ConFindr installed, use the command `confindr_database_setup` to have ConFindr download the latest version
of the rMLST databases. This script takes two arguments - a `-s` where you give the path to the text file containing your consumer 
key and secret, and a `-o` to specify where you want the sequences downloaded. Only the `-s` is mandatory. If your output
directory is not specified, ConFindr will first search for an environmental variable called `CONFINDR_DB`, and if it can't
find that it will automatically download to a folder called `.confindr_db` in your home directory.

## Installing Using Conda (Recommended)

ConFindr is available within bioconda - to get bioconda installed and running see instructions [here](https://bioconda.github.io/).

With bioconda running, you can install ConFindr with the following command:

`conda install -c bioconda confindr`

With that done, typing `confindr.py` will bring access the ConFindr pipeline. See the [Usage](usage.md) section for instructions on how to use ConFindr, including a ConFindr run on an example dataset.

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
- [Mash (>=2.0)](https://github.com/marbl/Mash/releases)
- [KMA (>=1.2.0)](https://bitbucket.org/genomicepidemiology/kma)
- [Python (>=3.5)](https://www.python.org/downloads/)

If you want to run in Nanopore mode, you'll also need to get [minimap2](https://github.com/lh3/minimap2).

Instructions on adding programs to your $PATH can be found [here](https://stackoverflow.com/questions/14637979/how-to-permanently-set-path-on-linux-unix).

If ConFindr can't find these dependencies when you try to run it, you will see an error message and the program will quit.



