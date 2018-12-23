# arcnumt
Scripts used to detect archaic NUMTs in NGS data
================================================

Description
-----------

A pipeline including a collection of scripts used to analyse NUMTs discovered in whole genome paired read data.
flanking_region_analysis.py is a script used to calculate match ratios with archaic genomes of a specific genomic region.
numt_stats.py is a script to calculate various statistics for discovered NUMTs.
mito_variance.py is a script to calculate pariwise differences between all sequences of an alignment.

Required resources
------------------

This workflow is based on the output of dinumt (https://github.com/mills-lab/dinumt) including the supplementary files obtained with the option --output_support.

For some steps third party software is required. Here is a list of those I used, but they can be replaced by other software doing the same:
* bam2fastx: https://github.com/PacificBiosciences/bam2fastx
* bwa mem: http://bio-bwa.sourceforge.net/
* samtools: http://www.htslib.org/doc/samtools.html
* bam-rewrap: https://bitbucket.org/ustenzel/biohazard-tools
* GATK 4.0 HaplotypeCaller: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php
* GATK 3.8 FastaAlternativeReferenceMaker: https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_fasta_FastaAlternateReferenceMaker.php


python modules needed:
* pysam
* argparse
* collections
* multiprocessing
* vcf
* Bio


In addition you will need:
* RSRS reference sequence: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3322232/
* alignment file with chimp and human mitochondrial genomes used for phylogenetic analysis. Must contain RSRS as a reference

Workflow
--------
Split the NUMT read files into single bam-files for each NUMT:
~~~
split_sam.py -s sample_support.sam \
~~~
Further process each individual NUMT to obtain its sequence and combine it with corresponding sequences from various mitochondrial genomes:
~~~
bam2fastx -q -A -o numt.fq -N numt.bam \
bwa mem RSRS.fasta numt.fq | samtools view -b | bam-rewrap RSRS:16569 | samtools sort > numt.sorted.bam; samtools index numt.sorted.bam \
fixbam.py -s numt.sorted.bam -o numt.sorted.fixed.bam; samtools index numt.sorted.fixed.bam \
getbed.py -s numt.sorted.fixed.bam -o numt.bed \

gatk4.0 HaplotypeCaller -L numt.bed -R RSRS.fasta -I numt.sorted.fixed.bam -O numt.vcf \
java -jar gatk3.8 -T FalstaAlternateReferenceMaker -R RSRS.fasta -o numt.fasta -L numt.bed -V numt.vcf \

extract_mito.py -n numt.fasta -a aligned_mt_genomes.fasta -b numt.bed -o numt_mito.fasta
~~~
