#! /usr/bin/env python

# splits a supplementary file from dinumt containing mitochondrial reads of
# NUMTs into files for individual NUMTs according to the mapping of their mate
# on the reference

import pysam
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--samfile", action="store",
                    help="sam file")
args = parser.parse_args()


samfile = pysam.AlignmentFile(args.samfile)

iter = samfile.fetch()
MT_reads = {}
NUMTS = defaultdict(list)
count = 0

# index the reads
for read in iter:
    if read.reference_name == "MT":
        count += 1
        MT_reads[count] = read

# merge reads belonging to one numt
for entry in MT_reads.keys():
    x = 0
    read = MT_reads[entry]
    start = read.next_reference_start
    # numts are more than 5000 bp away from each other on the reference
    end = read.next_reference_start + 5000
    for a in NUMTS.keys():
        if (a.split("_")[0] == read.next_reference_name and
            int(a.split("_")[1]) <= end and
            int(a.split("_")[2]) >= start):
            if int(a.split("_")[1]) < start:
                start = int(a.split("_")[1])
            if int(a.split("_")[2]) > end:
                end = int(a.split("_")[2])
            numt_name = str(read.next_reference_name)\
                + "_" + str(start) + "_" + str(end)
            NUMTS[numt_name] = NUMTS.pop(a)
            NUMTS[numt_name].append(entry)
            x = 1
    if x == 0:
        numt_name = str(read.next_reference_name)\
            + "_" + str(start)\
            + "_" + str(end)
        NUMTS[numt_name].append(entry)

for numt in NUMTS.keys():
    outfile_name = args.samfile.split(".")[0] + "_" + numt.split('_')[0] + '_'
    + numt.split('_')[1][:4] + ".bam"
    outfile = pysam.AlignmentFile(outfile_name, "wb", template=samfile)
    for read in NUMTS[numt]:
        outfile.write(MT_reads[read])
print(args.samfile + " processed")
