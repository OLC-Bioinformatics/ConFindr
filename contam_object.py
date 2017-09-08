
class ContamObject:
    def __init__(self, forward_reads):
        self.forward_reads = forward_reads
        self.reverse_reads = 'NA'
        self.unique_kmers = -1
        self.snv_count = -1
        self.trimmed_forward_reads = 'NA'
        self.trimmed_reverse_reads = 'NA'
        self.forward_rmlst_reads = 'NA'
        self.reverse_rmlst_reads = 'NA'
        self.subsampled_forward = 'NA'
        self.subsampled_reverse = 'NA'
        self.jellyfish_file = 'NA'
        self.mer_fasta = 'NA'
        self.samfile = 'NA'
