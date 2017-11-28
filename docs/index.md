# What is ConFindr?

ConFindr is a pipeline that can detect contamination in bacterial NGS data, both between and within species. It can do this with exceptional sensitivity - two samples mixed together with as few as 500 SNPs between them (> 99.9 percent identity!) can be identified. This allows for stringent quality control of NGS samples.

# How Does ConFindr Work?

### Brief Overview

ADD IMAGE HERE EVENTUALLY

ConFindr works by looking at rMLST genes. These 53 genes are known to be single copy and conserved across all bacteria, making them excellent markers. As they are known to be single copy (with some caveats), any sample that has multiple alleles of one or more rMLST gene is likely to be contaminated. To identify the presence of multiple alleles in a sample, the following workflow is followed:

1. Determine the genus of each sample so that genus-specific rMLST databases can be constructed.
2. Perform stringent quality trimming and bait out reads that contain rMLST gene sequence, using BBDuk.
3. Subsample rMLST reads to a depth of approximately 20X.
4. Split reads into kmers using jellyfish - typically to a size of k=31. 
5. Compare all of the kmers found to all other kmers found, looking for pairs of kmers that differ only by one substitution - these are assumed to represent multiple alleles of the same gene.
6. Repeat steps 3-5 a few times (typically 5) and take the median number of contaminating kmers from these repetitions.
7. Use Mash to do a quick check for cross-species contamination.

### Detailed Overview & FAQ

##### Why Subsample Reads?

In our testing, the coverage depth of samples we tested varied wildly. Some samples had as little coverage as 20-30X, while others were well over 200X. We found that samples with very high coverage (over 100X) frequently had false positive results coming up, as some kmers were present often enough due to sequencing error to look like another allele was present. Enforcing a higher cutoff for number of kmers present on the high depth samples caused drastic cuts to sensitivity. Subsampling reads to 20X coverage prevents sequencing errors while maintaining good sensitivity.

##### Why Not Just Map Reads to the rMLST Database?

The obvious way to approach this problem would be to map reads to the rMLST database and look for hits to multiple alleles of the same gene. This approach was attempted in early iterations of the pipeline, but was found to work fairly poorly. The issue here was that many rMLST alleles are similar enough that the same read will map to many of them. Because of this, we would often only have one or two reads that mapped unambiguously to rMLST genes across a sample - not enough to reliably call contamination.


