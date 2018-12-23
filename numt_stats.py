#! /usr/bin/env python3

import glob
import pysam
from collections import defaultdict
import numpy as np
import os.path
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--info-file', action='store', help='''file with
                    sample names in first and sample coverage in 10th
                    column''')
parser.add_argument('-d', '--directory-input', action='store', help='''
                    directory with alll files from numt analysis ending with
                    /''')
parser.add_argument('-o', '--outfile', action='store')

args = parser.parse_args()

NUMT_STAT = defaultdict(list)
numt_id = 0
SAMPLE_COVERAGE = {}
with open(args.info_file, "r") as handle:
    for line in handle.readlines():
        sample = line.split()[0]
        cov = line.split()[11]
        SAMPLE_COVERAGE[sample] = cov


numtfiles = args.directory_input + '*.fixed.rg.bam'
count = 0
for f in glob.glob(numtfiles):
    print("Analyze numt nr. " + str(numt_id))
    numt_id += 1
    seq_len = 0
    gc_count = 0
    samfile = pysam.AlignmentFile(f, "rb")
    name = os.path.basename(f).split('_')
    numt_sample = '_'.join(name[:-3])
    numt = "_".join(name[-3:-1])
    numt_chrom = numt.split("_")[0]
    fasta_file = f.split(".")[0] + ".numt.fa"

    # sequence might be split in multiple segments
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq_len += len(record.seq)
            gc_count += record.seq.count('G')
            gc_count += record.seq.count('C')
    if seq_len > 0:
        gc_content = gc_count / seq_len
    else:
        gc_content = 0
    READ_GC = []
    READ_LEN = []
    for read in samfile.fetch():
        sequence = read.query_sequence
        READ_GC.append(sequence.count("G") + sequence.count("C"))
        READ_LEN.append(len(sequence))
    if len(READ_LEN) == 0:
        continue

    gc_content2 = (sum(READ_GC) / sum(READ_LEN))
    COVERAGE = [pileupcolumn.n for pileupcolumn in samfile.pileup()]
    n = len(COVERAGE)
    mean = round(np.average(COVERAGE), 2)
    std = round(np.std(COVERAGE), 2)
    am_below_5 = round(100 * len([a for a in COVERAGE if a <= 5]) / n, 1)
    # if numt_chrom in ["X", "Y"]:
    #     sample_cov = "?"
    #     per_sample_cov = "?"
    sample_cov = SAMPLE_COVERAGE[numt_sample]
    per_sample_cov = 100 * (float(mean) / float(sample_cov))

    NUMT_STAT[numt_id] = [numt_sample, numt, mean, std, am_below_5,
                          gc_content, seq_len, sample_cov, per_sample_cov]


outfile = open(args.outfile, "w+")
outfile.write("SAMPLE\tNUMT\tMEAN_COV\tSTD\tN_BELOW_5\tGC_CONTENT\t"
              "SEQ_LEN\tSAMPLE_COV\tPER_SAMPLE_COV\n")
for key in NUMT_STAT.keys():
    outfile.write("\t".join([str(a) for a in NUMT_STAT[key]])
                  + "\n")
