#! /usr/bin/env python
# creates a bed file for the regions where reads map

import pysam
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--samfile", action="store",
                    default="/mnt/scratch/rob/output/results_dinumt/"
                    "pap_merged_refmitos_mapped.sam",
                    help="sam file")
parser.add_argument("-o", "--outfile", action="store")
args = parser.parse_args()


samfile = pysam.AlignmentFile(args.samfile, "rb")

ref_start = {}
positions = []
start = []
end = []

for pileupcolumn in samfile.pileup():
    if pileupcolumn.n > 4:
        positions.append(pileupcolumn.pos)

for f in positions:
    if (f - 1) not in positions:
        start.append(f)
    if (f + 1) not in positions:
        end.append(f + 1)

numt = 0
for a in start:
    ref_start[a] = end[numt]
    numt += 1

bedfile = open(args.outfile, "w+")

start = list(ref_start.keys())
name = 0
start.sort()
for a in start:
    name += 1
    bedfile.write("RSRS"
                  + "\t" + str(a) + "\t" + str(ref_start[a]) + "\t" +
                  "numt_" + str(name) + "\n")
print('regions written to bedfile')
