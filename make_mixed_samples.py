import os
import glob

pairs = [['2014-SEQ-0157', '2014-SEQ-1305'], ['2014-SEQ-0448', '2015-SEQ-1646'], ['2014-SEQ-0671', '2016-SEQ-1048'],
         ['2014-SEQ-1305', '2014-SEQ-1310'], ['2015-SEQ-1046', '2015-SEQ-1155'], ['2015-SEQ-1213', '2016-SEQ-0761'],
         ['2015-SEQ-1213', '2016-SEQ-0766'], ['2015-SEQ-1256', '2017-SEQ-0284']]

coverage = dict()
coverage['2014-SEQ-0157'] = 35
coverage['2014-SEQ-0448'] = 82
coverage['2014-SEQ-0671'] = 93
coverage['2014-SEQ-1305'] = 58
coverage['2015-SEQ-1046'] = 34
coverage['2015-SEQ-1213'] = 96
coverage['2015-SEQ-1256'] = 134

levels = ['0.01', '0.02', '0.05', '0.1', '0.15', '0.2']

for pair in pairs:
    forward_contam = glob.glob(pair[1] + '*R1*')
    reverse_contam = glob.glob(pair[1] + '*R2*')
    forward_source = glob.glob(pair[0] + '*R1*')
    reverse_source = glob.glob(pair[0] + '*R2*')
    for level in levels:
        desired_bases = int(5200000.0 * float(coverage[pair[0]]) * float(level))
        cmd = 'reformat.sh in1={} in2={} out1={} out2={} samplebasestarget={} overwrite'.format(forward_contam[0],
                                                                                                reverse_contam[0],
                                                                                                'contam_R1.fastq.gz',
                                                                                                'contam_R2.fastq.gz',
                                                                                                str(desired_bases))
        os.system(cmd)
        cmd = 'zcat {} {} > {}'.format(forward_source[0], 'contam_R1.fastq.gz', pair[0] + '_' + pair[1] + '_R1_' + level + '.fastq')
        os.system(cmd)
        cmd = 'zcat {} {} > {}'.format(reverse_source[0], 'contam_R2.fastq.gz', pair[0] + '_' + pair[1] + '_R2_' + level + '.fastq')
        os.system(cmd)
