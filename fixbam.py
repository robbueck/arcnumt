#! /usr/bin/env python

import pysam
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--samfile", action="store",
                    default="/mnt/scratch/rob/output/results_dinumt/"
                    "pap_merged_refmitos_mapped.sam",
                    help="sam file")
parser.add_argument('-o', '--outfile', action='store')
args = parser.parse_args()


head = {'SQ': [{'LN': 16569, 'SN': 'RSRS'}],
        'RG': [{'PU': 'unit1', 'LB': 'lib1', 'ID': '4',
                'PL': 'illumina', 'SM': '20'}],
        'PG': [{'PN': 'bam-rewrap', 'ID': 'bam-rewrap',
                'CL': 'RSRS:16569\tVN:0.2\tPP:bwa'},
               {'PN': 'bwa', 'ID': 'bwa', 'VN': '0.7.17-r1188',
                'CL': '/home/robert_buecking/Programs/bwa/bwa mem -t 3 '
                '/home/robert_buecking/References/mitorefs/RSRS/'
                'circular_RSRS.fasta UV043_3_13848303_13858818.fq'}],
        'HD': {'SO': 'coordinate', 'VN': '1.5'}}
samfile = pysam.AlignmentFile(args.samfile, "rb")
outfile = pysam.AlignmentFile(args.outfile, "wb", header=head)
iter = samfile.fetch()

for record in iter:
    record.set_tag('RG', '4', 'Z')
    if record.next_reference_id == -1:
        record.next_reference_start = -1
    outfile.write(record)
print("Finished")
