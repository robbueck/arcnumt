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
parser.add_argument('-p', '--numt-phase', action='store', help='''file with
                    phase information for archaic numt. first column sample
                    name, secon column phase''')
parser.add_argument('-t', '--threads', action='store', type=int)
parser.add_argument('-r', '--region', action='store', help='''chromosomal
                    position, e.g. 3:1384800-1385800''')
parser.add_argument('-o', '--outfiles-prefix', action='store')
parser.add_argument('-i', '--info-file', action='store', help='''file with
                    sample informations in columns. first row name of
                    categories, first column sample names''')
parser.add_argument('-T', '--temp-dir', action='store')

args = parser.parse_args()


output = mp.Queue()

FREQ_FILES = {}
MATCH_FILES = []
SAMPLES_FILE = []
for input_file in args.worldwide_freq:
    name = os.path.basename(input_file).split('.')[0]
    # FREQ_FILES[name] = vcf.Reader(filename=input_file)
    FREQ_FILES[name] = input_file

if args.worldwide_match:
    for input_file in args.worldwide_match:
        # MATCH_FILES.append(vcf.Reader(filename=input_file))
        MATCH_FILES.append(input_file)
else:
    MATCH_FILES = list(FREQ_FILES.values())

for input_file in args.samples_vcf:
    # SAMPLES_FILE.append(vcf.Reader(filename=input_file))
    SAMPLES_FILE.append(input_file)

den_file = vcf.Reader(filename=args.archaic_vcf)

PHASE = {}
phasefile = open(args.numt_phase, 'r')
for line in phasefile.readlines():
    sample, phase = line.rstrip().split('\t')
    PHASE[sample] = int(phase)

STUDIES = MATCH_FILES + SAMPLES_FILE
# for stud in STUDIES:
#     for sample in stud.samples:
#         if sample not in PHASE.keys():
#             PHASE[sample] = 3

chrom, region = args.region.split(':')
chrom = int(chrom)
start, end = region.split('-')
start = int(start)
end = int(end)
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


def get_sample_stats(intake):
    """
    creates dictionary for all sites between a modern humand and a denisovan
    with 0: no denisova match and 1: denisova match
    """
    # declaring variables
    n = intake[0]
    vcf_hum = vcf.Reader(filename=intake[1])
    DEN_VS = {}
    unphased = []
    for record in den_file.fetch(chrom, start, end):
        DEN_VS[record.POS] = record.genotype('Denisova')
    temp_name = args.temp_dir + 'temp_file' + str(n)
    temp = open(temp_name, 'w+')
    print('analzing SAMPLES_FILE: ' + vcf_hum.filename)
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
    # output.put((n, temp_name))
    return (n, temp_name)


def allel_freq(vcf_file_name, positions):
    '''takes a vcf_file and the positions where modern humans and dens match.
    for alleles where modern humans
    and denisova match, report the allele frequency in vcf-file'''
    vcf_file = vcf.Reader(filename=vcf_file_name)
    FREQ = {}
    for pos in set(positions):
        FREQ[pos] = 0
        den = 'jq'
        for record in den_file.fetch(chrom, pos - 1, pos):
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
    FREQS = {}
    sample = sample_matching_sites[0]
    matching_sites = sample_matching_sites[1]
    for study in FREQ_FILES:
        freq_count = allel_freq(study, matching_sites)
        FREQS[study] = [a for a in freq_count.keys() if freq_count[a] <= 0.05]
    print('finished low_freqs')
    return(sample, FREQS)


################################################################
# run
###############################################################

# dictionaries for results for output
MATCH_RATIO = {}


print('start analyzing match ratios')
# get matching positions
pool = mp.Pool(processes=args.threads)
results = pool.map(get_sample_stats, [(x, STUDIES[x]) for x in
                                      range(len(STUDIES))])

# processes = [mp.Process(target=get_sample_stats, args=(x, STUDIES[x]))
#              for x in range(len(STUDIES))]
# for p in processes:
#     p.start()
#
# for p in processes:
#     p.join()
#
# results = [output.get() for p in processes]


# if program fails after here, just comment lines above and uncomment lines
# below to prevent starting analysis all over again. change names of tmep files
# also comment os.remove
# at the end of the script so temp files don't get removed
# results = [(0, '/home/robert_buecking/flanking_region/temp_file0'),
#            (1, '/home/robert_buecking/flanking_region/temp_file1'),
#            (2, '/home/robert_buecking/flanking_region/temp_file2'),
#            (3, '/home/robert_buecking/flanking_region/temp_file3'),
#            (4, '/home/robert_buecking/flanking_region/temp_file4'),
#            (5, '/home/robert_buecking/flanking_region/temp_file5')]
results.sort()

# calculate match ratios
TEMP_FILES = [r[1] for r in results]
MATCHING_SITES = []
MATCHING_SITES_PER_SAMPLE = {}
INFO_SAMPLE = {}
for study in TEMP_FILES:
    f = open(study, 'r')
    for l in f.readlines():
        line = l.rstrip().split('\t')
        # line = line.split('\t')
        sample, phase, sites = line
        sites = sites.split(';')
        if sample in PHASE.keys() and int(phase) == PHASE[sample]:
            sample = 'X_' + sample
        else:
            sample = 'O_' + sample
        sample = sample + '_' + str(phase)

        matches = [int(x.split(':')[1]) for x in sites]
        MATCH_RATIO[sample] = (str(np.average(matches)) + '\t'
                               + str(len(matches)))

        MATCHING_SITES_PER_SAMPLE[sample] = list()
        for x in sites:
            x = x.split(':')
            if int(x[1]) == 1:
                MATCHING_SITES.append(int(x[0]))
                MATCHING_SITES_PER_SAMPLE[sample].append(int(x[0]))
    f.close()


print('finished match ratio analysis, start analyzing matching sites')


# counting alleles matching denisova with low frequency in global studies

print('computing low freqs in ' + (', ').join(FREQ_FILES.keys()))

LOW_SAMPLE_ALLEL_FREQS = {}
SAMPLE_ALLEL_FREQS = {}
current = 0
samps = MATCHING_SITES_PER_SAMPLE.keys()
lowsamples = []
lsamp = 0
for sample in samps:
    ALLEL_FREQS = {}
    LOW_ALLEL_FREQS = {}
    if current in range(0, len(samps), 200):
        print(str(len(samps) - current) + ' Samples left to analyse')
    for f in FREQ_FILES.keys():
        ALLEL_FREQS[f] = allel_freq(FREQ_FILES[f],
                                    MATCHING_SITES_PER_SAMPLE[sample])
        LOW_ALLEL_FREQS[f] = len([a for a in ALLEL_FREQS[f].values()
                                  if a <= 0.1])
    SAMPLE_ALLEL_FREQS[sample] = ALLEL_FREQS
    LOW_SAMPLE_ALLEL_FREQS[sample] = LOW_ALLEL_FREQS
    current += 1
    # if LOW_ALLEL_FREQS[list(FREQ_FILES.keys())[0]] > 0:
    #     if lsamp == 0:
    #         lsamp = sample
    #     print(f + '\t' + sample)
    #     print(lsamp)
    #     print(LOW_SAMPLE_ALLEL_FREQS[lsamp])
    #     lowsamples.append(sample)
    #     print((LOW_SAMPLE_ALLEL_FREQS))
# for stu in FREQ_FILES.keys():
#     for sample in samps:
#         low = LOW_SAMPLE_ALLEL_FREQS[sample][stu]
#         if low > 0:
#             print(stu + sample + str(low))
#         if sample in lowsamples:
#             print('lowsamples: ' + sample)

#############################################################
# write to files
print("writing data to files")

SAMPLE_INFO = {}
with open(args.info_file, 'r') as f:
    info_fields = f.readline().split('\t')[1:]
    for l in f:
        line = l.split("\t")
        SAMPLE_INFO[line[0]] = line[1:]

match_filename = args.outfiles_prefix + '_match_ratios.txt'
with open(match_filename, "w+") as f:
    f.write('SAMPLE\tMATCH_RATIO\tSITES\tLOW_FREQS_' +
            ('\tLOW_FREQS_').join(sorted(FREQ_FILES.keys())) + '\t'
            + ('\t').join(info_fields))
    for sample in samps:
        LOWS = [str(LOW_SAMPLE_ALLEL_FREQS[sample][s]) for s in
                sorted(FREQ_FILES.keys())]
        if sample in lowsamples:
            print(LOWS)
            print(LOW_SAMPLE_ALLEL_FREQS[sample])
        f.write(sample + "\t" + str(MATCH_RATIO[sample]) + "\t"
                + ('\t').join(LOWS) + '\t'
                + ("\t").join(SAMPLE_INFO[sample[2:-2]]))

den_freq_filename = args.outfiles_prefix + '_den_freqs.txt'
with open(den_freq_filename, 'w+') as f:
    for sample in MATCHING_SITES_PER_SAMPLE.keys():
        positions = sorted(MATCHING_SITES_PER_SAMPLE[sample])
        f.write('\n' + sample + '\nSTUDY\t' + '\t'.join(map(str, positions))
                + '\n')
        for study in sorted(FREQ_FILES.keys()):
            f.write(study + '\t' +
                    '\t'.join(map(str,
                                  [SAMPLE_ALLEL_FREQS[sample][study][a] for
                                   a in positions]))
                    + '\n')

# for study in TEMP_FILES:
#     os.remove(study)
