# What is ConFindr?

ConFindr is a pipeline that can detect contamination in bacterial NGS data, both between and within species. It can do this with exceptional sensitivity - two samples mixed together with as few as
500 SNPs between them (> 99.9 percent identity!) can be identified. This allows for stringent quality control of NGS samples.

## How Does ConFindr Work?

ConFindr works by looking at rMLST genes. These 53 genes are known to be single copy and conserved across all bacteria, making them excellent markers. As they are known to be single copy (with some caveats), any sample that has multiple alleles of one or more rMLST gene is likely to be contaminated. To identify the presence of multiple alleles in a sample, the following workflow is followed:

1. Use Mash to determine the genus of each sample so that genus-specific rMLST databases can be constructed
and check for interspecies contamination.
2. Perform stringent quality trimming and bait out reads that contain rMLST gene sequence, using BBDuk.
3. Subsample rMLST reads to a depth of approximately 20X.
4. Split reads into kmers using jellyfish - typically to a size of k=31. 
5. Compare all of the kmers found to all other kmers found, looking for pairs of kmers that differ only by one substitution
 - these are assumed to represent multiple alleles of the same gene, and therefore contamination, since
 these genes are only single copy.
6. Repeat steps 3-5 a few times (typically 3) and take the median number of contaminating kmers from these repetitions,
in order to avoid any issues with an incredibly unlucky subsample.

## ConFindr Intra-species Contamination Detection Performance

The below graph shows the magnitude of contamination detected for several synthetic datasets in *Escherichia coli*.
Strains were mixed together that were either identical (and so should have no contamination), not identical but
have the same rMLST type (and so are contaminated, but beyond the limit of detection), two strains with the same serotype,
and therefore very closely related, or two strains of differing serotypes.

The black line on the graph represents our cutoff for contaminated samples - any sample with a magnitude above that
can reliably be called contaminated. As can be seen, two different serotypes from the same species are reliably detected
at contamination levels of 5 percent or higher, two strains of the same serotype are often detected at 10 percent
contamination and almost always at 20 percent contamination, and two strains that have the same rMLST or are identical
never have contamination detected. Results are very similar for other species.

![alt text](performance.png "ConFindr Performance")

## ConFindr Inter-species Contamination Detection

ConFindr seems to reliably be able to detect interspecies contamination at levels of 5 percent or above, but its checks
for interspecies contamination are not particularly rigorous. If you're very worried about interspecies contamination,
it would be a good idea to put your samples through some sort of metagenomics software (such as Kraken or Kaiju).