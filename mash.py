#!/usr/bin/env python
from subprocess import call
from threading import Thread
from accessoryFunctions.accessoryFunctions import *

__author__ = 'adamkoziol'


class Mash(object):
    def sketching(self):
        printtime('Indexing files for {} analysis'.format(self.analysistype), self.starttime)
        # Create the threads for the analysis
        for i in range(self.cpus):
            threads = Thread(target=self.sketch, args=())
            threads.setDaemon(True)
            threads.start()
        # Populate threads for each gene, genome combination
        for sample in self.metadata:
            # Create the analysis type-specific GenObject
            setattr(sample, self.analysistype, GenObject())
            # Set attributes
            sample[self.analysistype].reportdir = os.path.join(sample.general.outputdirectory, self.analysistype)
            make_path(sample[self.analysistype].reportdir)
            sample[self.analysistype].targetpath = self.referencefilepath if not self.pipeline else os.path.join(
                self.referencefilepath, self.analysistype)
            sample[self.analysistype].refseqsketch = \
                sample[self.analysistype].targetpath + '/RefSeqSketchesDefaults.msh'
            sample[self.analysistype].sketchfilenoext = '{}/{}'.format(sample[self.analysistype].reportdir,
                                                                       sample.name)
            sample[self.analysistype].sketchfile = sample[self.analysistype].sketchfilenoext + '.msh'
            # Make the mash output directory if necessary
            make_path(sample[self.analysistype].reportdir)
            # Create a file containing the path/name of the filtered, corrected fastq files
            sample[self.analysistype].filelist = '{}/{}_fastqfiles.txt'.format(sample[self.analysistype].reportdir,
                                                                               sample.name)
            with open(sample[self.analysistype].filelist, 'w') as filelist:
                filelist.write('\n'.join(sample.general.trimmedcorrectedfastqfiles))

            # Create the system call
            sample.commands.sketch = 'mash sketch -m 2 -p {} -l {} -o {}' \
                .format(self.cpus, sample[self.analysistype].filelist, sample[self.analysistype].sketchfilenoext)
            # Add each sample to the threads
            try:
                self.sketchqueue.put(sample)
            except (KeyboardInterrupt, SystemExit):
                printtime('Received keyboard interrupt, quitting threads', self.starttime)
                quit()
        # Join the threads
        self.sketchqueue.join()
        self.mashing()

    def sketch(self):
        while True:
            sample = self.sketchqueue.get()
            if not os.path.isfile(sample[self.analysistype].sketchfile):
                call(sample.commands.sketch, shell=True, stdout=self.fnull, stderr=self.fnull)
            self.sketchqueue.task_done()

    def mashing(self):
        printtime('Performing {} analyses'.format(self.analysistype), self.starttime)
        # Create the threads for the analysis
        for i in range(self.cpus):
                threads = Thread(target=self.mash, args=())
                threads.setDaemon(True)
                threads.start()
        # Populate threads for each gene, genome combination
        for sample in self.metadata:
            sample[self.analysistype].mashresults = '{}/{}.tab'.format(sample[self.analysistype].reportdir,
                                                                       sample.name)

            sample.commands.mash = \
                'mash dist -p {} {} {} | sort -gk3 > {}'.format(self.threads,
                                                                sample[self.analysistype].refseqsketch,
                                                                sample[self.analysistype].sketchfile,
                                                                sample[self.analysistype].mashresults)
            try:
                self.mashqueue.put(sample)
            except (KeyboardInterrupt, SystemExit):
                printtime('Received keyboard interrupt, quitting threads', self.starttime)
                quit()
        # Join the threads
        self.mashqueue.join()
        self.parse()

    def mash(self):
        while True:
            sample = self.mashqueue.get()
            if not os.path.isfile(sample[self.analysistype].mashresults):
                # call(sample.commands.mash, shell=True, stdout=self.fnull, stderr=self.fnull)
                call(sample.commands.mash, shell=True)
            self.mashqueue.task_done()

    def parse(self):
        import re
        printtime('Determining closest refseq genome', self.starttime)
        for sample in self.metadata:
            # Open the results and extract the first line of data
            mashdata = open(sample[self.analysistype].mashresults).readline().rstrip()
            # Split on tabs
            data = mashdata.split('\t')
            try:
                referenceid, queryid, sample[self.analysistype].mashdistance, sample[self.analysistype]. \
                    pvalue, sample[self.analysistype].nummatches = data
                # The database is formatted such that the reference file name is usually preceded by '-.-'
                # e.g. refseq-NZ-1005511-PRJNA224116-SAMN00794588-GCF_000303935.1-.-Escherichia_coli_PA45.fna
                #      refseq-NZ-1639-PRJNA224116-SAMN03349770-GCF_000951975.1-p3KSM-Listeria_monocytogenes.fna
                sample[self.analysistype].closestrefseq = re.findall(r'.+-(.+)\.fna', referenceid)[0]
                sample[self.analysistype].closestrefseqgenus = sample[self.analysistype].closestrefseq.split('_')[0]
                sample[self.analysistype].closestrefseqspecies = sample[self.analysistype].closestrefseq.split('_')[1]
                # Set the closest refseq genus - will be used for all typing that requires the genus to be known
                sample.general.referencegenus = sample[self.analysistype].closestrefseqgenus
            except ValueError:
                # Populate the attribute with negative results
                sample[self.analysistype].closestrefseqgenus = 'NA'
        # Create the report
        self.reporter()

    def reporter(self):
        make_path(self.reportpath)
        header = 'Strain,ReferenceGenus,ReferenceFile,ReferenceGenomeMashDistance,Pvalue,NumMatchingHashes\n'
        data = ''
        for sample in self.metadata:
            try:
                data += '{},{},{},{},{},{}\n'.format(sample.name,
                                                     sample[self.analysistype].closestrefseqgenus,
                                                     sample[self.analysistype].closestrefseq,
                                                     sample[self.analysistype].mashdistance,
                                                     sample[self.analysistype].pvalue,
                                                     sample[self.analysistype].nummatches)
            except AttributeError:
                data += '{}\n'.format(sample.name)
        # Create the report file
        reportfile = '{}/mash.csv'.format(self.reportpath)
        with open(reportfile, 'w') as report:
            report.write(header)
            report.write(data)

    def __init__(self, inputobject, analysistype):
        from queue import Queue
        self.metadata = inputobject.runmetadata.samples
        self.referencefilepath = inputobject.reffilepath
        self.starttime = inputobject.starttime
        self.reportpath = inputobject.reportpath
        self.cpus = inputobject.cpus
        self.threads = int(self.cpus / len(self.metadata)) if self.cpus / len(self.metadata) > 1 else 1
        self.sketchqueue = Queue(maxsize=self.cpus)
        self.mashqueue = Queue(maxsize=4)
        self.analysistype = analysistype
        self.pipeline = inputobject.pipeline
        self.fnull = open(os.devnull, 'w')  # define /dev/null
        self.sketching()
