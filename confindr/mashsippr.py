#!/usr/bin/env python
import subprocess
import time
from confindr.mash import *
from sipprCommon.objectprep import Objectprep
from accessoryFunctions.accessoryFunctions import *

__author__ = 'adamkoziol'


class MashSippr(object):

    def runner(self):
        """
        Run the necessary methods in the correct order
        """
        printtime('Starting mashsippr analysis pipeline', self.starttime)
        if not self.pipeline:
            # Create the objects to be used in the analyses
            objects = Objectprep(self)
            objects.objectprep()
            self.runmetadata = objects.samples
        # Run the analyses
        Mash(self, self.analysistype)

    def __init__(self, args, pipelinecommit, startingtime, scriptpath):
        """
        :param args: command line arguments
        :param pipelinecommit: pipeline commit or version
        :param startingtime: time the script was started
        :param scriptpath: home path of the script
        """
        import multiprocessing
        # Initialise variables
        self.commit = str(pipelinecommit)
        self.starttime = startingtime
        self.homepath = scriptpath
        # Define variables based on supplied arguments
        if args.path.endswith('/'):
            self.path = args.path
        else:
            self.path = os.path.join(args.path, '')
        assert os.path.isdir(self.path), u'Supplied path is not a valid directory {0!r:s}'.format(self.path)
        if args.sequencepath.endswith('/'):
            self.sequencepath = args.sequencepath
        else:
            self.sequencepath = os.path.join(args.sequencepath, '')
        assert os.path.isdir(self.sequencepath), u'Sequence path  is not a valid directory {0!r:s}' \
            .format(self.sequencepath)
        if args.targetpath.endswith('/'):
            self.targetpath = args.targetpath
        else:
            self.targetpath = os.path.join(args.targetpath, '')
        self.reportpath = os.path.join(self.path, 'reports')
        assert os.path.isdir(self.targetpath), u'Target path is not a valid directory {0!r:s}' \
            .format(self.targetpath)
        # ref file path is used to work with sub module code with a different naming scheme
        self.reffilepath = self.targetpath
        self.bcltofastq = args.bcltofastq
        self.miseqpath = args.miseqpath
        self.miseqfolder = args.miseqfolder
        self.fastqdestination = args.fastqdestination
        self.forwardlength = args.forwardlength
        self.reverselength = args.reverselength
        self.numreads = 2 if self.reverselength != 0 else 1
        self.customsamplesheet = args.customsamplesheet
        # Use the argument for the number of threads to use, or default to the number of cpus in the system
        self.cpus = int(args.cpus if args.cpus else multiprocessing.cpu_count())
        self.pipeline = args.pipeline
        if self.pipeline:
            self.runmetadata = args.runmetadata
        else:
            self.runmetadata = MetadataObject()
        self.analysistype = 'mash'
        self.copy = False
        # Run the analyses
        self.runner()

if __name__ == '__main__':
    # Argument parser for user-inputted values, and a nifty help menu
    from argparse import ArgumentParser
    # Get the current commit of the pipeline from git
    # Extract the path of the current script from the full path + file name
    homepath = os.path.split(os.path.abspath(__file__))[0]
    # Find the commit of the script by running a command to change to the directory containing the script and run
    # a git command to return the short version of the commit hash
    commit = subprocess.Popen('cd {} && git rev-parse --short HEAD'.format(homepath),
                              shell=True, stdout=subprocess.PIPE).communicate()[0].rstrip()
    # Parser for arguments
    parser = ArgumentParser(description='Perform modelling of parameters for GeneSipping')
    parser.add_argument('path',
                        help='Specify input directory')
    parser.add_argument('-s', '--sequencepath',
                        required=True,
                        help='Path of .fastq(.gz) files to process.')
    parser.add_argument('-t', '--targetpath',
                        required=True,
                        help='Path of target files to process.')
    parser.add_argument('-n', '--cpus',
                        help='Number of threads. Default is the number of cores in the system')
    parser.add_argument('-b', '--bcltofastq',
                        action='store_true',
                        help='Optionally run bcl2fastq on an in-progress Illumina MiSeq run. Must include:'
                             'miseqpath, and miseqfolder arguments, and optionally readlengthforward, '
                             'readlengthreverse, and projectName arguments.')
    parser.add_argument('-m', '--miseqpath',
                        help='Path of the folder containing MiSeq run data folder')
    parser.add_argument('-f', '--miseqfolder',
                        help='Name of the folder containing MiSeq run data')
    parser.add_argument('-d', '--fastqdestination',
                        help='Optional folder path to store .fastq files created using the fastqCreation module. '
                             'Defaults to path/miseqfolder')
    parser.add_argument('-r1', '--forwardlength',
                        default='full',
                        help='Length of forward reads to use. Can specify "full" to take the full length of '
                             'forward reads specified on the SampleSheet')
    parser.add_argument('-r2', '--reverselength',
                        default='full',
                        help='Length of reverse reads to use. Can specify "full" to take the full length of '
                             'reverse reads specified on the SampleSheet')
    parser.add_argument('-c', '--customsamplesheet',
                        help='Path of folder containing a custom sample sheet (still must be named "SampleSheet.csv")')
    parser.add_argument('-P', '--projectName',
                        help='A name for the analyses. If nothing is provided, then the "Sample_Project" field '
                             'in the provided sample sheet will be used. Please note that bcl2fastq creates '
                             'subfolders using the project name, so if multiple names are provided, the results '
                             'will be split as into multiple projects')
    # Get the arguments into an object
    arguments = parser.parse_args()
    arguments.pipeline = False
    # Define the start time
    start = time.time()

    # Run the script
    MashSippr(arguments, commit, start, homepath)

    # Print a bold, green exit statement
    print('\033[92m' + '\033[1m' + "\nElapsed Time: %0.2f seconds" % (time.time() - start) + '\033[0m')


class PipelineInit(object):

    def __init__(self, inputobject):
        inputobject.pipeline = True
        MashSippr(inputobject, inputobject.commit, inputobject.starttime, inputobject.homepath)
