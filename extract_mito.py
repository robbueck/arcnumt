#! /usr/bin/env python

from Bio import SeqIO
import argparse
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

parser = argparse.ArgumentParser()
parser.add_argument("-n", "--numt", action="store",
                    help="fasta file with numt")
parser.add_argument("-a", "--alignment", action="store",
                    default="/home/robert_buecking/References/mitorefs/"
                    "aligned_mito_101humans_RSRS_10ancienthumans_17nea_4den"
                    "_sima_pan.fa",
                    help="alignment of numts in fasta")
parser.add_argument("-b", "--bed-file", action="store",
                    help="bed file with mitochondrial coordinates")
parser.add_argument("-o", "--outfile", action="store")
args = parser.parse_args()


# get numts and corresponding mito seqs and write them to file


# update positions on reference
def update_pos(pos):
    start = pos
    change = 0.0001
    while change > 0:
        oldstart = start
        start = pos + ref_seq[:start].count("-")
        change = start - oldstart
    return start


# read files and prepare lists
reference = list(SeqIO.parse(args.alignment, "fasta"))
numts = list(SeqIO.parse(args.numt, "fasta"))
with open(args.bed_file, "r") as bed_file:
    coordinates = bed_file.readlines()

for record in reference:
    if record.id == "RSRS":
        ref_seq = record.seq

nr = 0
NUMT = {}

if len(coordinates) > 0:
    for n in coordinates:
        nr += 1
        NUMT[nr] = n.split("\t")
        NUMT[nr][1] = update_pos(int(NUMT[nr][1]))
        NUMT[nr][2] = update_pos(int(NUMT[nr][2]))
    if NUMT[1][1] <= 1 and NUMT[nr][2] >= 16568:
        numts.append(numts.pop(0))
        NUMT[nr + 1] = NUMT[1]
        del NUMT[1]
    numt1 = numts[-1]
    numt1.seq = sum([nu.seq for nu in numts], Seq("", generic_dna))
    outfile = (args.outfile)
    with open(outfile, "w+") as handle:
        SeqIO.write(numt1, handle, "fasta")
        for record in SeqIO.parse(args.alignment, "fasta"):
            RECORDS = []
            for key in sorted(NUMT.keys()):
                start = NUMT[key][1]
                end = NUMT[key][2]
                RECORDS.append(record.seq[start:end])
            record.seq = sum(RECORDS, Seq("", generic_dna))
            SeqIO.write(record, handle, "fasta")
print("seqs written to " + args.numt.split(".")[0] + "_mito.fasta")
