#! /usr/bin/env python3
'''
calculates pairwise differences between all sequences of an alignment file
sequences are grouped in Denisova, Neanderthal and Modern human. Den and Nea
seqs are identified by keywords in their id, all other seqs are treated as
modern humans. NUMT seq name must be a number between 0 an 19 to be identified.
'''
from Bio import SeqIO
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--alignment-file", action="store")
parser.add_argument('-o', '--outfile', action='store')
args = parser.parse_args()

################################################################
# defining functions
################################################################


def which_human(human_seq):
    if "denisova" in human_seq.id.lower():
        group = "DEN"
    elif "neanderthal" in human_seq.id.lower():
        group = "NEA"
    elif human_seq.id in [str(a) for a in range(0, 20, 1)]:
        group = "NUM"
    # elif "heidelbergensis" in human_seq.id.lower():
    #     group = "SIM"
    # elif 'pan' in human_seq.id.lower():
    #     group = 'PAN'
    else:
        group = "SAP"
    return group


def seq_diff(seq1, seq2):
    differences = 0
    N = range(0, len(seq1))
    for a in N:
        if seq1[a] != seq2[a]:
            differences += 1
    return differences


def histogram(l, m_dif):
    hist = []
    for a in range(0, m_dif + 1):
        hist.append(str(len([b for b in l if b == a])))
    return hist


##############################################################
# run
##############################################################

# remove nonhuman or reference seqs
TRASH = ["RSRS", "rCRS_NC_012920116569/116569",
         "NC_002083.1_PongoA1602616499/116499",
         "NC_011120.1_GorillaA1592316412/116412",
         'NC_001643.1_PanA1598616554/116554',
         'KF683087.1_Homo_heidelbergensis']
alignment = [a for a in list(SeqIO.parse(args.alignment_file, "fasta"))
             if a.id not in TRASH]
DIFFERENCES = defaultdict(list)
HISTOGRAM = {}
max_diff = 0

for record1 in alignment:
    alignment = [a for a in alignment if a.id != record1.id]
    hum1 = which_human(record1)
    for record2 in alignment:
        hum2 = which_human(record2)
        CATS = [hum1, hum2]
        CATS.sort()
        category = "_".join(CATS)
        difference = seq_diff(record1.seq, record2.seq)
        DIFFERENCES[category].append(difference)
        if difference > max_diff:
            max_diff = difference

for key in DIFFERENCES.keys():
    HISTOGRAM[key] = histogram(DIFFERENCES[key], max_diff)
CAT_LIST = list(HISTOGRAM.keys())
CAT_LIST.sort()


# write output file
outfile = open(args.outfile, "w+")
outfile.write('DIFF\t' + '\t'.join(CAT_LIST) + '\n')
for d in range(0, max_diff + 1):
    outfile.write(str(d) + '\t' +
                  '\t'.join([HISTOGRAM[a][d] for a in CAT_LIST]) + '\n')
