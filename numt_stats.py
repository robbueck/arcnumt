#! /usr/bin/env python3

import glob
import pysam
from collections import defaultdict
import numpy as np
import os.path
import re
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--info-file', action='store', help='''file with
                    sample names in first and sample coverage in 10th
                    column''')
parser.add_argument('-d', '--directory-input', action='store', help='''
                    directory with alll files from numt analysis ending with
                    /''')

args = parser.parse_args()

NUMT_STAT = defaultdict(list)
numt_id = 0
SAMPLE_COVERAGE = {}
with open(args.info_file, "r") as handle:
    for line in handle.readlines():
        sample, cov = line.split()[0, 9]
        SAMPLE_COVERAGE[sample] = cov


numtfiles = args.directory_input + '*.fixed.rg.bam'
for f in glob.glob(numtfiles):
    print("Analyze numt nr. " + str(numt_id))
    numt_id += 1
    sample_cov = "bla"
    seq_len = 0
    samfile = pysam.AlignmentFile(f, "rb")
    numt_sample = re.split('_', os.path.basename(f))[0]
    numt = "_".join(os.path.basename(f).split("_")[-3:-1])
    numt_chrom = numt.split("_")[0]
    samfile = pysam.AlignmentFile(f, "rb")
    fasta_file = f.split(".")[0] + ".numt.fa"

    # sequence might be split in multiple segments
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq_len += len(record.seq)
    READ_GC = []
    READ_LEN = []
    for read in samfile.fetch():
        sequence = read.query_sequence
        READ_GC.append(sequence.count("G") + sequence.count("C"))
        READ_LEN.append(len(sequence))
    if len(READ_LEN) == 0:
        continue

    gc_content = (sum(READ_GC) / sum(READ_LEN))
    COVERAGE = [pileupcolumn.n for pileupcolumn in samfile.pileup()]
    n = len(COVERAGE)
    min_cov = min(COVERAGE)
    am_min = round(100 * COVERAGE.count(min_cov) / n, 1)
    mean = round(np.average(COVERAGE), 2)
    std = round(np.std(COVERAGE), 2)
    am_below_5 = round(100 * len([a for a in COVERAGE if a <= 5]) / n, 1)
    if numt_chrom in ["X", "Y"]:
        sample_cov = "?"
        per_sample_cov = "?"
    else:
        numt_sample_chr = numt_sample + "_" + str(numt_chrom)
        sample_cov = SAMPLE_COVERAGE[numt_sample_chr]
        per_sample_cov = 100 * (float(mean) / float(sample_cov))

    NUMT_STAT[numt_id] = [numt_sample, numt, min_cov, am_min,
                          mean, std, am_below_5,
                          gc_content, seq_len, sample_cov, per_sample_cov]


outfile = open("/mnt/scratch/rob/output/ind_results_dinumt/"
               "ind_numt_stats_corr_5.txt", "w+")
outfile.write("SAMPLE\tNUMT\tMIN_COV\tCOUNT_MIN_COV"
              "\tMEAN_COV\tSTD\tN_BELOW_5\tGC_CONTENT\t"
              "SEQ_LEN\tSAMPLE_COV\tPER_SAMPLE_COV\n")
for key in NUMT_STAT.keys():
    outfile.write("\t".join([str(a) for a in NUMT_STAT[key]])
                  + "\n")
