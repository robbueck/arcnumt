#! /usr/bin/env python3
# calculates a match ratio between different genomes and denisova
# to check, if match ratio for indiiduals with insertion is higher than for
# individuals without.
import vcf
import numpy as np
import multiprocessing as mp
from scipy.stats import binom_test
# import time
# import matplotlib.pyplot as plt

###############################################################
# get input files
##############################################################

output = mp.Queue()

ind_file = vcf.Reader(filename='/home/robert_buecking/References/'
                      'phased_flanking/ind_chr3_138_snp_only.vcf.gz')
pap_file = vcf.Reader(filename='/home/robert_buecking/References/'
                      'phased_flanking/pap_chr3_138_snp_only.vcf.gz')
gp_file = vcf.Reader(filename='/home/robert_buecking/References/'
                     'phased_flanking/gp_chr3_138_snp_only.vcf.gz')
gp_file_0 = vcf.Reader(filename='/home/robert_buecking/References/'
                       'phased_flanking/gp_chr3_138_snp_only_samples.vcf.gz')
gp_file_1 = vcf.Reader(filename='/home/robert_buecking/References/'
                       'phased_flanking/gp_chr3_138_snp_only_samples1.vcf.gz')
gp_file_2 = vcf.Reader(filename='/home/robert_buecking/References/'
                       'phased_flanking/gp_chr3_138_snp_only_samples2.vcf.gz')
gp_file_3 = vcf.Reader(filename='/home/robert_buecking/References/'
                       'phased_flanking/gp_chr3_138_snp_only_samples3.vcf.gz')
sgdp_file = vcf.Reader(filename='/home/robert_buecking/References/'
                       'phased_flanking/sgdp_chr3_138_snp_only.vcf.gz')
den_file = vcf.Reader(filename='/mnt/454/Vindija/high_cov/genotypes/'
                      'Denisova/chr3_mq25_mapab100.vcf.gz')


PHASE = {'UV043': 1, 'UV1260': 1, 'UV925': 1, 'UV927': 1, 'UV956': 0,
         'UV958': 0, 'UV971': 1, 'UV979': 1, 'S_Papuan-5': 0, 'S_Papuan-10': 0,
         "ALR06": 1, "ALR11": 0, "FAN-KEI-024": 0, "NG88-F": 1, "UV030": 1}
# PHASE = {'UV043': 1, 'UV1260': 1, 'UV925': 1, 'UV927': 1, 'UV956': 0,
#          'UV958': 0, 'UV971': 1, 'UV979': 1}
IND_SAMPLES = ind_file.samples
PAP_SAMPLES = pap_file.samples
SGDP_SAMPLES = sgdp_file.samples
GP_SAMPLES = gp_file.samples
ALL_SAMPLES = PAP_SAMPLES + SGDP_SAMPLES + GP_SAMPLES + IND_SAMPLES
for sample in ALL_SAMPLES:
    if sample not in PHASE.keys():
        PHASE[sample] = 3
STUDIES = [pap_file, sgdp_file, gp_file_0, gp_file_1, gp_file_2, gp_file_3,
           ind_file]

STUDY_DICT = {}
POP_DICT = {}
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

    # match = 0
    # if p == 3:
    #     den_set = set()
    #     for alt in den_snp.ALT:
    #         den_set.update(str(alt))
    #     hum_set = set(hum_snp.split("|"))
    #     if (den_set & hum_set):
    #         match = 1
    # else:
    #     if hum_snp.split("|")[p] in den_snp.ALT:
    #         match = 1


def get_sample_stats(n, vcf_hum):
    """
    creates dictionary for all sites between a modern humand and a denisovan
    with 0: no denisova match and 1: denisova match
    """
    # declaring variables
    DEN_VS = {}
    unphased = []
    chrom = 3
    start = 13851000 - 20000
    end = 13851000 + 20000
    if vcf_hum == ind_file:
        start = 13848625 - 10000
        end = 13848625 + 10000
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
        # if PHASE[sample] == 3:
        stats0 = check_match(DEN_VS, HUM_VS, 0, unphased)
        stats1 = check_match(DEN_VS, HUM_VS, 1, unphased)
        up0 = [stats0[x] for x in stats0.keys() if x >= 13851000]
        up1 = [stats1[x] for x in stats1.keys() if x >= 13851000]
        # print(up1)
        p_binom0 = binom_test(up0.count(1), len(up0), np.average(up1))
        p_binom1 = binom_test(up1.count(1), len(up1), np.average(up0))
        #   if p_binom0 > p_binom1:
        #       p_binom0 = p_binom1
        #   if np.average(up0) > np.average(up1):
        #       sample_stats[sample] = stats0
        #   else:
        #       sample_stats[sample] = stats1
        # else:
        #     sample_stats[sample] = check_match(DEN_VS, HUM_VS, PHASE[sample])
        #     p_binom0 = 0
        temp.write(sample + '\t0\t' + str(p_binom0) + '\t'
                   + ';'.join([str(a) + ':' + str(b) for a, b in
                              stats0.items()]) + '\n'
                   + sample + '\t1\t' + str(p_binom1) + '\t'
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
        # current += 1
        # if current in range(0, 6000, 200):
        #     print('finished ' + str(current) + ' Samples')
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
NUMT_ACC_UP = {}
NUMT_ACC_DOWN = {}
NUMT_JACK_UP = {}
NUMT_JACK_DOWN = {}
MATCH_RATIO_UP = {}
MATCH_RATIO_DOWN = {}
ALLEL_FREQS = {}
GP_ALLEL_FREQS_SAMPLE = {}
# PAP_ALLEL_FREQS_SAMPLE = {}
IND_ALLEL_FREQS_SAMPLE = {}
SGDP_ALLEL_FREQS_SAMPLE = {}

# positions = set()


# pool = mp.Pool(processes=2)
# results = [pool.apply(get_sample_stats, args=(x, ))
#            for x in range(2)]
# print(results)

print('start analyzing match ratios')

# processes = [mp.Process(target=get_sample_stats, args=(x, STUDIES[x]))
#              for x in range(7)]
# for p in processes:
#     p.start()
#
# for p in processes:
#     p.join()
#
# results = [output.get() for p in processes]
# if program fails after here, just comment lines above and uncomment lines
# below to prevent starting analysis all over again. also comment os.remove
# below so temp files don't get removed
results = [(0, '/home/robert_buecking/flanking_region/temp_file0'),
           (1, '/home/robert_buecking/flanking_region/temp_file1'),
           (2, '/home/robert_buecking/flanking_region/temp_file2'),
           (3, '/home/robert_buecking/flanking_region/temp_file3'),
           (4, '/home/robert_buecking/flanking_region/temp_file4'),
           (5, '/home/robert_buecking/flanking_region/temp_file5'),
           (6, '/home/robert_buecking/flanking_region/temp_file6')]
results.sort()
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
P_BINOM = {}
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
        P_BINOM[sample] = line[2]
        matches_up = [int(x.split(':')[1]) for x in line[3].split(';') if
                      int(x.split(':')[0]) <= 13851000]
        matches_down = [int(x.split(':')[1]) for x in line[3].split(';') if
                        int(x.split(':')[0]) > 13851000]
        MATCH_RATIO_UP[sample] = (str(np.average(matches_up)) + '\t'
                                  + str(len(matches_up)))
        MATCH_RATIO_DOWN[sample] = (str(np.average(matches_down)) + '\t'
                                    + str(len(matches_down)))
        # if sample in ['N_UV043_1', 'UV043_0', 'N_UV1260_1', 'N_UV1260_0',
        #               'UV925_0', 'N_UV925_1', 'UV927_0', 'N_UV927_1',
        #               'N_UV956_0', 'N_UV958_0', 'UV897_1', 'DNG09_1',
        #               'N_UV971_1', 'N_UV979_1', 'N_S_Papuan-5_0',
        #               'N_S_Papuan-10_0', 'N_ALR06_1', 'HG04099_0',
        #             'N_ALR11_0', 'N_FAN-KEI-0240', 'N_NG88-F_1', 'N_UV030_1',
        #               'NA19309', 'HG02651_1', 'HG00357_0', 'S_Igorot-2_1',
        #               'HG00119_1', 'HG04099_0', 'HG02541_1']:
        sites = line[3].split(';')
        MATCHING_SITES_PER_SAMPLE[sample] = list()
        for x in sites:
            if int(x.split(':')[1]) == 1:
                MATCHING_SITES.append(int(x.split(':')[0]))
                MATCHING_SITES_PER_SAMPLE[sample].append(
                    int(x.split(':')[0]))
    f.close()


print('finished match ratio analysis, start analyzing matching sites')

# SITES = sorted(set(MATCHING_SITES))
# ALLEL_FREQS['numt'] = {pos: (MATCHING_SITES.count(pos) / 15)
#                       for pos in SITES if pos >= 13838634 or pos <= 13858483}
# ALLEL_FREQS['numt'].update({pos: (MATCHING_SITES.count(pos) / 10)
#                             for pos in SITES if pos < 13838634
#                             and pos > 13858483})
# ALLEL_FREQS['sgdp'] = allel_freq(sgdp_file, MATCHING_SITES)
# # ALLEL_FREQS['pap'] = allel_freq(pap_file, MATCHING_SITES)
# ALLEL_FREQS['gp'] = allel_freq(gp_file, MATCHING_SITES)
# # ALLEL_FREQS['ind'] = allel_freq(ind_file, MATCHING_SITES)

# for sample in MATCHING_SITES_PER_SAMPLE.keys():
#     GP_ALLEL_FREQS_SAMPLE[sample] = allel_freq(
#         gp_file, MATCHING_SITES_PER_SAMPLE[sample])
#     SGDP_ALLEL_FREQS_SAMPLE[sample] = allel_freq(
GP_ALLEL_FREQS_SAMPLE = allel_freq(gp_file, MATCHING_SITES)
SGDP_ALLEL_FREQS_SAMPLE = allel_freq(sgdp_file, MATCHING_SITES)
#       sgdp_file, MATCHING_SITES_PER_SAMPLE[sample])
#     PAP_ALLEL_FREQS_SAMPLE[sample] = allel_freq(
#         pap_file, MATCHING_SITES_PER_SAMPLE[sample])
#     IND_ALLEL_FREQS_SAMPLE[sample] = allel_freq(
#         ind_file, MATCHING_SITES_PER_SAMPLE[sample])

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


LOW_FREQ_GP = {}
LOW_FREQ_SGDP = {}
current = 0
for sample in MATCHING_SITES_PER_SAMPLE.keys():
    current += 1
    if current in range(0, 6000, 200):
        print('finished ' + str(current) + ' Samples')
    freq_count_gp = allel_freq(gp_file, MATCHING_SITES_PER_SAMPLE[sample])
    freq_count_sgdp = allel_freq(sgdp_file, MATCHING_SITES_PER_SAMPLE[sample])
    LOW_FREQ_GP[sample] = [a for a in freq_count_gp.keys() if
                           freq_count_gp[a] <= 0.05 and a <= 13851000]
    LOW_FREQ_SGDP[sample] = [a for a in freq_count_sgdp.keys()
                             if freq_count_gp[a] <= 0.05 and a <= 13851000]


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
    f.write('SAMPLE\tMATCH_RATIO_UP\tSITES_UP\tMATCH_RATIO_DOWN\tSITES_DOWN\t'
            + 'P_BINOM_PHASE\tLOW_FREQ_GP\tLOW_FREQ_SGDP\tPOPULATION\tREGION'
            + '\tCOUNTRY\tSTUDY\n')
    for sample in MATCH_RATIO_UP.keys():
        f.write(sample + "\t" + str(MATCH_RATIO_UP[sample]) + "\t"
                + str(MATCH_RATIO_DOWN[sample]) + "\t" + P_BINOM[sample] + '\t'
                + str(len(LOW_FREQ_GP[sample])) + '\t'
                + str(len(LOW_FREQ_SGDP[sample])) + '\t'
                + ("\t").join(SAMPLE_INFO[INFO_SAMPLE[sample]]))

with open('/home/robert_buecking/flanking_region/den_freqs.txt', 'w+') as f:
    # f.write('STUDY' + '\t' + '\t'.join(map(str, SITES)) + '\n')
    # for study in ALLEL_FREQS.keys():
    #     f.write(study + '\t' + '\t'.join(map(str, [ALLEL_FREQS[study][a]
    #                                                for a in SITES])) + '\n')
    positions = sorted(list(GP_ALLEL_FREQS_SAMPLE.keys()))
    f.write('\n' + sample + '\n' + 'STUDY' + '\t'
            + '\t'.join(map(str, positions)) + '\n')
#         f.write('pap\t' + '\t'.join(map(str,
#                                         [PAP_ALLEL_FREQS_SAMPLE[sample][a]
#                                          for a in positions])) + '\n')
    f.write('sgdp\t' + '\t'.join(map(str,
                                     [SGDP_ALLEL_FREQS_SAMPLE[a]
                                      for a in positions])) + '\n')
    f.write('gp\t' + '\t'.join(map(str, [GP_ALLEL_FREQS_SAMPLE[a]
                                         for a in positions])) + '\n')
#         f.write('ind\t' + '\t'.join(map(str,
#                                         [IND_ALLEL_FREQS_SAMPLE[sample][a]
#                                          for a in positions])) + '\n')
# for study in TEMP_FILES:
#     os.remove(study)
