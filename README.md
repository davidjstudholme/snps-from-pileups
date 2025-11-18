Note: If you are looking for the code to accompany this paper: **Reyes-Herrera et al. (2023)**. Genome sequence data reveal at least two distinct incursions of the Tropical Race 4 variant of Fusarium wilt into South America. _Phytopathology_ 113: 90–97, then please go to https://github.com/davidjstudholme/fusarium_TR4_genomics.

# SNPsFromPileups

The purpose of this script is as part of a simple pipeline for parsing a set of
[SAMtools mpileup](http://samtools.sourceforge.net/pileup.shtml) files to identify SNPs and tabulate
the allele-frequencies at each SNP across a set of biological samples.

The scenario is that we have performed genomic re-sequencing on a set of several biological samples.
The genomic sequence reads have been aligned against a reference genome sequence using a tool such as [BWA](https://github.com/lh3/bwa) to generate a ```.bam``` file.
We want to identify single-nucleotide positions in the genome that show variation between samples (i.e. SNPs) and we want to estimate the allele frequencies in each sample at each of those SNP sites.
* The samples could be individuals, for example each could be an individual plant or animal. If these individuals are diploid, then we expect that the allele frequencies will always be 0, 0.5 or 1. If the individuals are triploid then we would expect 0, 0.33, 0.67, or 1. And so on for other levels of ploidy. There is a discrete set of possible values for allele frequency.
* Alternatively, the samples could be populations rather than individuals; for example, they could be microbial cultures containing many millions of individuals of various different genotypes. In this case, the allele frequency depends on the relative abundance of each genotype within the population. There is a continuous set of possible values for allele frequency.

So, let's assume that we have one ```.bam``` file for each sample (genomic sequence reads aligned against reference genome ```genome.fasta```).
The first step is to identify candidate SNPs using [SAMtools](http://www.htslib.org/) version 1.6 and
[BCFtools](https://samtools.github.io/bcftools/bcftools.html) version 1.6: 

```
for alignmentFile in *.bam; do samtools mpileup -u -f genome.fasta $alignmentFile > $alignmentFile.bcf; done

for alignmentFile in *.bam; do bcftools call -m -v -Ov $alignmentFile.bcf > $alignmentFile.vcf; done

```

So, now we have a set of ```.vcf``` files containing candidate SNPs. Let's again use version 1.6 of BCFtools to filter these to keep only high-confidence ones:
```
for alignmentFile in *.bam; do bcftools filter --SnpGap 100 --include '(REF="A" | REF="C" | REF="G" | REF="T") & %QUAL>=35 & MIN(DP)>=5 & INDEL=0' $alignmentFile.vcf > $alignmentFile.filtered.vcf; done
```

This ```bcftools filter``` step eliminates indels with low-confidence single-nucleotide variant calls.
It also eliminates candidate SNVs within 10 base pairs of an indel, since alignment artefacts are relatively common in the close vicinity of indels.

Our ```filtered.vcf``` files now contain a catalogue of relatively high-confidence SNPs.

We also need the alignments in SAMtools' mpileup format. The ```SNPsFromPileups.pl``` script will extract the allele frequencies from these files. These files are potentially very large and parsing every line is potentially slow. Therefore, the ```filtered.vcf``` files are used to reduce the search space. That is, ```SNPsFromPileups.pl``` ignores any genomic positions that are not listed in at least one of the ```filtered.vcf``` files. (Note that we could have used the ```-l``` option in ```samtools mpileup``` if the ```filtered.vcf``` were converted into ```.bed``` format).

To generate the mpileup files (```.pileup```):

```
for alignmentFile in *.bam; do samtools mpileup -f genome.fasta $alignmentFile > $alignmentFile.pileup; done
```

Now we are ready to generate the table of allele frequencies (considering only sites where read-coverage is at least 10x):

```
perl get_snps_from_pileups.pl 10 *.filtered.vcf *.pileup > snps.csv
```

The [get_snps_from_pileups.pl](get_snps_from_pileups.pl) script attempts to be efficient in its use of memory and compute and can handle genome szies of hundred of Mbp.
However, if you are working with much smaller genome sizes, e.g. bacterial genomes of a few Mbp, then you can instead use [get_snps_from_pileups_small_genome.pl](get_snps_from_pileups_small_genome.pl):

```
perl get_snps_from_pileups_small_genome.pl 10 *.filtered.vcf *.pileup > snps.csv
```

Finally, you can use [get_haplotypes_and_aligned_fasta_from_csv.pl](get_haplotypes_and_aligned_fasta_from_csv.pl)
to convert the SNPs table (snps.csv) into a Nexus-formatted file that can serve as input into a tools such as PopArt:
 
```
perl get_haplotypes_and_aligned_fasta_from_csv.pl snps.csv
```
