# arcnumt
Scripts used to detect archaic NUMTs in NGS dataDetection of archaic NUMTs in NGS data.
=======================================================================================

Description
-----------

A pipeline including a collection of scripts used to analyse NUMTs discovered in whole genome paired read data.

Required resources
------------------

This workflow is based on the output of dinumt (https://github.com/mills-lab/dinumt) including the supplementary files obtained with the option --output_support.

For some steps third party software is required. Here is a list of those I used with its function, but they can be replaced by other software doing the same:
* bam2fastx: https://github.com/PacificBiosciences/bam2fastx
* bwa mem: http://bio-bwa.sourceforge.net/
* samtools: http://www.htslib.org/doc/samtools.html
* bam-rewrap: https://bitbucket.org/ustenzel/biohazard-tools


python modules needed:
* pysam
* argparse
* collections
*


In addition you will need:
* RSRS reference sequence: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3322232/
* 
To be uploaded soon.

Workflow
--------

~~~
split_sam.py -s sample_support.sam \


~~~

