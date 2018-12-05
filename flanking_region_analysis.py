#! /usr/bin/env python3
# calculates a match ratio between different genomes and denisova
# to check, if match ratio for indiiduals with insertion is higher than for
# individuals without.
import vcf
import numpy as np
import multiprocessing as mp
import argparse
import os.path
# import matplotlib.pyplot as plt

###############################################################
# get input files
##############################################################

parser = argparse.ArgumentParser()
parser.add_argument('-w', '--worldwide-freq', nargs='+', action='store',
                    help='''list of phased vcfs from populations
                    without introgression''')
parser.add_argument('-m', '--worldwide-match', nargs='*',
                    help='''vcf files to calculate match ratios for not
                    introgressed populations. large vcfs can be split by
                    samples for faster parallel processing. if not specified,
                    files given with -w are used''')
parser.add_argument('-s', '--samples-vcf', nargs='+', action='store',
                    help='list of phased vcfs from samples to analyse')
parser.add_argument('-a', '--archaic-vcf', action='store',
                    help='vcf of putative introgressed archaic human')
parser.add_argument('-p', '--numt_phase', action='store', help='''file with
                    phase information for archaic numt. first column sample
                    name, secon column phase''')
parser.add_argument('-t', '--threads', action='store')
parser.add_argument('-r', '--region', action='store', help='''chromosomal
                    position, e.g. 3:1384800-1385800''')

args = parser.parse_args()


output = mp.Queue()

FREQ_FILES = {}
MATCH_FILES = []
SAMPLES_FILE = []
for input_file in args.wordwide_freq:
    name = os.path.basename(input_file).split('.')[0]
    FREQ_FILES[name] = vcf.Reader(filename=input_file)

if args.wordwide_match:
    for input_file in args.wordwide_freq:
        MATCH_FILES.add(vcf.Reader(filename=input_file))
else:
    MATCH_FILES = FREQ_FILES

for input_file in args.samples_vcf:
    SAMPLES_FILE.add(vcf.Reader(filename=input_file))

den_file = vcf.Reader(filename=args.archaic_vcf)

PHASE = {}
phasefile = open(args.numt_phase, 'r')
for line in phasefile.readlines():
    sample, phase = line.rstrip().split('\t')
    PHASE[sample] = phase

STUDIES = MATCH_FILES + SAMPLES_FILE
for stud in STUDIES:
    for sample in stud.samples:
        if sample not in PHASE.keys():
            PHASE[sample] = 3

chrom, region = args.region.split(':')
start, end = region.split('-')
##############################################################
# defenition of functions
##############################################################


def flat_list(list_of_lists):
    # take all values of list of list and put them in one list
    try:
        flat_list = [item for sublist in list_of_lists for item in sublist]
    except TypeError:
        print("Not a list of lists")
        flat_list = list_of_lists
    return flat_list


def check_match(den_snps, hum_snps, p, no_phase):
    """
    checks if denisovan SNP also found in sample
    for samples with numt insertion: just phase with insertion
    """
    stats = {}
    POSITIONS = [a for a in den_snps.keys() if den_snps[a].gt_type != 0]
    # exclude reference alleles
    POSITIONS += [a for a in hum_snps.keys() if
                  hum_snps[a].data[0].split('|')[p] == 1]
    for pos in POSITIONS:
        if pos not in den_snps.keys():
            continue
            # exclude unsequenced sites in DEN
        if pos in no_phase:
            continue
            # exclude unphased sites
        if pos not in hum_snps.keys():
            stats[pos] = 0
            continue
        hum = hum_snps[pos].gt_bases.split('|')
        if hum_snps[pos].data[0].split('|')[p] == 0:
            stats[pos] = 0
            continue
            # human is reference
        if den_snps[pos].gt_type == 0:
            stats[pos] = 0
            continue
        if hum[p] in den_snps[pos].gt_bases.split('/'):
            stats[pos] = 1
        else:
            stats[pos] = 0
    if not stats:
        return {0: 0}
    else:
        return stats


def get_sample_stats(n, vcf_hum):
    """
    creates dictionary for all sites between a modern humand and a denisovan
    with 0: no denisova match and 1: denisova match
    """
    # declaring variables
    DEN_VS = {}
    unphased = []
    for record in den_file.fetch(chrom, start, end):
        DEN_VS[record.POS] = record.genotype('Denisova')
    temp_name = '/home/robert_buecking/flanking_region/temp_file' + str(n)
    temp = open(temp_name, 'w+')
    print('analzing file: ' + vcf_hum.filename)
    for sample in vcf_hum.samples:
        HUM_VS = {}
        for record in vcf_hum.fetch(chrom, start, end):
            if record.genotype(sample).phased:
                # record is phased and not homozygous reference
                HUM_VS[record.POS] = (record.genotype(sample))
            else:
                unphased.append(record.POS)
        stats0 = check_match(DEN_VS, HUM_VS, 0, unphased)
        stats1 = check_match(DEN_VS, HUM_VS, 1, unphased)
        temp.write(sample + '\t0\t'
                   + ';'.join([str(a) + ':' + str(b) for a, b in
                              stats0.items()]) + '\n'
                   + sample + '\t1\t'
                   + ';'.join([str(a) + ':' + str(b) for a, b in
                               stats1.items()]) + '\n')
    temp.close()
    print('finished analyzing file: ' + vcf_hum.filename)
    output.put((n, temp_name))


def allel_freq(vcf_file, positions):
    '''takes a vcf_file and the positions where modern humans and dens match.
    for alleles where modern humans
    and denisova match, report the allele frequency in vcf-file'''
    FREQ = {}
    for pos in set(positions):
        FREQ[pos] = 0
        den = 'jq'
        for record in den_file.fetch(3, pos - 1, pos):
            if record.POS == pos:
                den = record.ALT
                if len(den) > 1:
                    raise Exception('too many alt-dens')
        for record in vcf_file.fetch(3, pos - 1, pos):
            if record.POS == pos:
                if den[0] in record.ALT:
                    x = record.ALT.index(den[0])
                    FREQ[pos] = record.aaf[x]
    return FREQ


def low_freqs(sample_matching_sites):
    '''counts how many matching alleles have a frequency below 0.05 in each
    study file'''
    sample = sample_matching_sites[0]
    matching_sites = sample_matching_sites[1]
    freq_count_sgdp = allel_freq(sgdp_file, matching_sites)
    freq_count_gp = allel_freq(gp_file, matching_sites)
    FREQ_GP = [a for a in freq_count_gp.keys()
               if freq_count_gp[a] <= 0.05 and a <= 13851000]
    FREQ_SGDP = [a for a in freq_count_sgdp.keys()
                 if freq_count_sgdp[a] <= 0.05 and a <= 1385100]
    # output_gp.put(FREQ)
    # output_sgdp.put(FREQ)
    print('finished low_freqs')
    return((sample, FREQ_GP, FREQ_SGDP))


################################################################
# run
###############################################################

# dictionaries for results for output
NUMT_SUM_STATS = {}
MATCH_RATIO = {}
GP_ALLEL_FREQS_SAMPLE = {}
PAP_ALLEL_FREQS_SAMPLE = {}
IND_ALLEL_FREQS_SAMPLE = {}
SGDP_ALLEL_FREQS_SAMPLE = {}


# pool = mp.Pool(processes=2)
# results = [pool.apply(get_sample_stats, args=(x, ))
#            for x in range(2)]
# print(results)

print('start analyzing match ratios')

processes = [mp.Process(target=get_sample_stats, args=(x, STUDIES[x]))
             for x in range(7)]
for p in processes:
    p.start()

for p in processes:
    p.join()

results = [output.get() for p in processes]
# if program fails after here, just comment lines above and uncomment lines
# below to prevent starting analysis all over again. also comment os.remove
# below so temp files don't get removed
# results = [(0, '/home/robert_buecking/flanking_region/temp_file0'),
#            (1, '/home/robert_buecking/flanking_region/temp_file1'),
#            (2, '/home/robert_buecking/flanking_region/temp_file2'),
#            (3, '/home/robert_buecking/flanking_region/temp_file3'),
#            (4, '/home/robert_buecking/flanking_region/temp_file4'),
#            (5, '/home/robert_buecking/flanking_region/temp_file5'),
#            (6, '/home/robert_buecking/flanking_region/temp_file6')]
# results.sort()
# SNPS = {k: v for d in results for k, v in d.items()}
# positions.update(set(a for a in SNPS[sample].keys()))

# for sample in SNPS.keys():
#     print('llll')
#     MATCH_RATIO_UP[sample] = np.average(list(
#         SNPS[sample][x] for x in SNPS[sample].keys() if x <= 13848625)
#     )
#     MATCH_RATIO_DONW[sample] = np.average(list(
#         SNPS[sample][x] for x in SNPS[sample].keys() if x > 13848625)
#     )

TEMP_FILES = [r[1] for r in results]
MATCHING_SITES = []
MATCHING_SITES_PER_SAMPLE = {}
INFO_SAMPLE = {}
for study in TEMP_FILES:
    f = open(study, 'r')
    for l in f.readlines():
        line = l.rstrip()
        line = line.split('\t')
        sample = line[0]
        phase = line[1]
        if int(phase) == PHASE[sample]:
            sample = 'N_' + sample
        sample = sample + '_' + str(phase)
        INFO_SAMPLE[sample] = line[0]
        matches = [int(x.split(':')[1]) for x in line[3].split(';')]
        MATCH_RATIO[sample] = (str(np.average(matches)) + '\t'
                               + str(len(matches)))
        sites = line[3].split(';')
        MATCHING_SITES_PER_SAMPLE[sample] = list()
        for x in sites:
            if int(x.split(':')[1]) == 1:
                MATCHING_SITES.append(int(x.split(':')[0]))
                MATCHING_SITES_PER_SAMPLE[sample].append(
                    int(x.split(':')[0]))
    f.close()


print('finished match ratio analysis, start analyzing matching sites')


# counting alleles matching denisova with low frequency in global studies

print('computing low freqs in 1000GP and SGDP')

# LOW_FREQ_GP = {}
# LOW_FREQ_SGDP = {}
# pool = mp.Pool(processes=40)
# results = pool.imap_unordered(low_freqs,
#                               list(MATCHING_SITES_PER_SAMPLE.items()))
# pool.close()
#
# while (True):
#     completed = results._index
#     if (completed == len(MATCHING_SITES_PER_SAMPLE)):
#         break
#     print("Waiting for", len(MATCHING_SITES_PER_SAMPLE) - completed,
#           "tasks to complete...")
#     time.sleep(50)
#
# for x in results:
#     LOW_FREQ_GP[x[0]] = x[1]
#     LOW_FREQ_SGDP[x[0]] = x[2]

LOW_SAMPLE_ALLEL_FREQS = {}
SAMPLE_ALLEL_FREQS = {}
LOW_ALLEL_FREQS = {}
ALLEL_FREQS = {}
current = 0
for sample in MATCHING_SITES_PER_SAMPLE.keys():
    ALLEL_FREQS = {}
    current += 1
    if current in range(0, 6000, 200):
        print('finished ' + str(current) + ' Samples')
    for f in FREQ_FILES.keys():
        ALLEL_FREQS[f] = allel_freq(FREQ_FILES[f], MATCHING_SITES)
        LOW_ALLEL_FREQS[f] = len([a for a in ALLEL_FREQS[f].values()
                                  if a <= 0.05])
    SAMPLE_ALLEL_FREQS[sample] = ALLEL_FREQS
    LOW_SAMPLE_ALLEL_FREQS[sample] = LOW_ALLEL_FREQS


#############################################################
# write to files
print("writing data to files")

SAMPLE_INFO = {}
with open("/home/robert_buecking/References/all_studies_samples.txt",
          'r') as f:
    for l in f:
        line = l.split("\t")
        SAMPLE_INFO[line[0]] = line[1:]

with open("/home/robert_buecking/flanking_region/match_ratios.txt", "w+") as f:
    f.write('SAMPLE\tMATCH_RATIO\tSITES\t\tLOW_FREQS_' +
            ('\tLOW_FREQS_').join(sorted(FREQ_FILES.keys())) +
            '\tPOPULATION\tREGION\tCOUNTRY\tSTUDY\n')
    for sample in MATCH_RATIO.keys():
        LOWS = [LOW_SAMPLE_ALLEL_FREQS[sample][f] for f in
                sorted(FREQ_FILES.keys())]
        f.write(sample + "\t" + str(MATCH_RATIO[sample]) + "\t"
                + LOWS.join('\t') + '\t'
                + ("\t").join(SAMPLE_INFO[INFO_SAMPLE[sample]]))

with open('/home/robert_buecking/flanking_region/den_freqs.txt', 'w+') as f:
    for sample in MATCHING_SITES_PER_SAMPLE.keys():
        positions = sorted(SAMPLE_ALLEL_FREQS.keys())
        f.write('\n' + sample + '\nSTUDY\t' + '\t'.join(map(str, positions))
                + '\n')
        for study in sorted(FREQ_FILES.keys()):
            f.write(study + '\t' +
                    '\t'.join(map(str,
                                  [SAMPLE_ALLEL_FREQS[sample][study][a] for
                                   a in positions]))
                    + '\n')

for study in TEMP_FILES:
    os.remove(study)
