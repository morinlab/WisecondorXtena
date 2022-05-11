# WisecondorX

import os

import numpy as np

from wisecondorX.overall_tools import exec_R, get_z_score, get_median_segment_variance, get_cpa

'''
Writes plots.
'''


def exec_write_plots(rem_input, results):
    json_plot_dir = os.path.abspath(rem_input['args'].outid + '_plot_tmp')
    json_dict = {
        'R_script': str('{}/include/plotter.R'.format(rem_input['wd'])),
        'ref_gender': str(rem_input['ref_gender']),
        'beta': str(rem_input['args'].beta),
        'zscore': str(rem_input['args'].zscore),
        'binsize': str(rem_input['binsize']),
        'n_reads': str(rem_input['n_reads']),
        'cairo': str(rem_input['args'].cairo),
        'results_r': results['results_r'],
        'results_w': results['results_w'],
        'results_c': results['results_c'],
        'ylim': str(rem_input['args'].ylim),
        'infile': str('{}.json'.format(json_plot_dir)),
        'out_dir': str('{}.plots'.format(rem_input['args'].outid)),
    }
    exec_R(json_dict)


'''
Calculates zz-scores, marks aberrations and
writes tables.
'''


def generate_output_tables(rem_input, results):
    _generate_bins_bed(rem_input, results)
    _generate_segments_and_aberrations_bed(rem_input, results)
    _generate_seg_file(rem_input, results)
    _generate_chr_statistics_file(rem_input, results)


def _generate_bins_bed(rem_input, results):
    bins_file = open('{}_bins.bed'.format(rem_input['args'].outid), 'w')
    bins_file.write('chr\tstart\tend\tid\tratio\ttypeR\tzscore\ttypeZ\n')

    results_r = results['results_r']
    results_z = results['results_z']
    binsize = rem_input['binsize']

    for chr in range(len(results_r)):
        chr_name = str(chr + 1)
        if chr_name == '23':
            chr_name = 'X'
        if chr_name == '24':
            chr_name = 'Y'
        feat = 1
        for i in range(len(results_r[chr])):
            r = results_r[chr][i]
            tR = results_r[chr][i]
            z = results_z[chr][i]
            tZ = results_z[chr][i]
            if r == 0:
                r = 'nan'
            else:
                if tR < 0:
                    tR = 'loss'
                else:
                    tR = 'gain'
            if z == 0:
                z = 'nan'
            else:
                if tZ < 0:
                    tZ = 'loss'
                else:
                    tZ = 'gain'
            feat_str = '{}:{}-{}'.format(chr_name, str(feat), str(feat + binsize - 1))
            row = [chr_name, feat, feat + binsize - 1, feat_str, r, tR, z, tZ]
            bins_file.write('{}\n'.format('\t'.join([str(x) for x in row])))
            feat += binsize

    bins_file.close()


def _generate_segments_and_aberrations_bed(rem_input, results):
    segments_file = open('{}_segments.bed'.format(rem_input['args'].outid), 'w')
    abberations_file = open('{}_aberrations.bed'.format(rem_input['args'].outid), 'w')
    segments_file.write('chr\tstart\tend\tratio\tzscore\n')
    abberations_file.write('chr\tstart\tend\tratio\tzscore\ttype\tcopy.number\n')

    for segment in results['results_c']:
        chr_name = str(segment[0] + 1)
        if chr_name == '23':
            chr_name = 'X'
        if chr_name == '24':
            chr_name = 'Y'
        row = [chr_name,
               int(segment[1] * rem_input['binsize'] + 1),
               int(segment[2] * rem_input['binsize']),
               segment[4], segment[3]]
        segments_file.write('{}\n'.format('\t'.join([str(x) for x in row])))
        ploidy = 2

        if (chr_name == 'X' or chr_name == 'Y') and rem_input['actual_gender'] == 'M':
            ploidy = 1
        if rem_input['args'].beta is not None:
            if __get_aberration_cutoff(rem_input['args'].beta, ploidy)[1] < float(segment[4]) < __get_aberration_cutoff(rem_input['args'].beta, ploidy)[2]:
                abberations_file.write('{}\tgain\t3\n'.format('\t'.join([str(x) for x in row])))
            elif float(segment[4]) < __get_aberration_cutoff(rem_input['args'].beta, ploidy)[0]:
                abberations_file.write('{}\tloss\t1\n'.format('\t'.join([str(x) for x in row])))
            elif float(segment[4]) > __get_aberration_cutoff(rem_input['args'].beta, ploidy)[2]:
                abberations_file.write('{}\tamp\t4\n'.format('\t'.join([str(x) for x in row])))
        elif isinstance(segment[3], str):
            continue
        else:
            if float(segment[3]) > rem_input['args'].zscore: #and float(segment[3]) < rem_input['args'].zscore * 5:
                abberations_file.write('{}\tgain\t3\n'.format('\t'.join([str(x) for x in row])))
            elif float(segment[3]) < - rem_input['args'].zscore:
                abberations_file.write('{}\tloss\t1\n'.format('\t'.join([str(x) for x in row])))
            #elif float(segment[3]) > rem_input['args'].zscore * 5:
            #    abberations_file.write('{}\tamp\t4\n'.format('\t'.join([str(x) for x in row])))

    segments_file.close()
    abberations_file.close()


def __get_aberration_cutoff(beta, ploidy):
    loss_cutoff = np.log2((ploidy - (beta / 2)) / ploidy)
    gain_cutoff = np.log2((ploidy + (beta / 2)) / ploidy)
    amp_cutoff = np.log2((ploidy + ((beta + 1)/2)) / ploidy)
    return loss_cutoff, gain_cutoff, amp_cutoff

def _generate_seg_file(rem_input, results):
    seg_file = open('{}.seg'.format(rem_input['args'].outid), 'w')
    seg_file.write('ID\tchr\tstart\tend\tLOH_flag\ttype\tcopy.number\tlog.ratio\n')

    for segment in results['results_c']:
        ID = os.path.basename(rem_input['args'].outid)
        LOH_flag = 0
        #type = {'loss': 'loss', 'gain': 'gain', 'neutral': 'neutral'} -> tried making a dictionary and then incorporating that into the fstrings below but it didnt work properly so I had to put each variable separately as seen below. 
        loss = 'loss'
        gain = 'gain'
        neutral = 'neut'
        amplification = 'amp'

        #copy_number = {'one': 1, 'two': 2, 'three': 3}
        one = 1
        two = 2
        three = 3
        four = 4
        start = int(segment[1] * rem_input['binsize'] + 1)
        end = int(segment[2] * rem_input['binsize'])

        chr_name = str(segment[0] + 1)
        if chr_name == '23':
            chr_name = 'X'
        if chr_name == '24':
            chr_name = 'Y'

        row = [ID, chr_name,
               start,
               end,
               int(LOH_flag),
               segment[4]]
        #seg_file.write('{}\n'.format('\t'.join([str(x) for x in row])))
        # commented this out because it causes the lines to duplicate. ^
        ploidy = 2


        #the original seg_file.write code commented out below here was used before incorporating f strings, but these do not allow log ratio  to be placed at the end as a column, so f strings were used and each column had to be made into a variable (as seen above)
        #added these fstrings on July 14, 2021
        if (chr_name == 'X' or chr_name == 'Y') and rem_input['actual_gender'] == 'M':
            ploidy = 1
        if rem_input['args'].beta is not None:
            if __get_aberration_cutoff(rem_input['args'].beta, ploidy)[1] < float(segment[4]) < __get_aberration_cutoff(rem_input['args'].beta, ploidy)[2]:
                #seg_file.write('{}\tgain\t3\n'.format('\t'.join([str(x) for x in row])))
                seg_file.write(f"{ID}\t{chr_name}\t{start}\t{end}\t{LOH_flag}\t{gain}\t{three}\t{segment[4]}\n")
            elif float(segment[4]) < __get_aberration_cutoff(rem_input['args'].beta, ploidy)[0]:
                #seg_file.write('{}\tloss\t1\n'.format('\t'.join([str(x) for x in row])))
                seg_file.write(f"{ID}\t{chr_name}\t{start}\t{end}\t{LOH_flag}\t{loss}\t{one}\t{segment[4]}\n")
            elif float(segment[4]) > __get_aberration_cutoff(rem_input['args'].beta, ploidy)[2]:
                #seg_file.write('{}\tamp\t4\n'.format('\t'.join([str(x) for x in row])))
                seg_file.write(f"{ID}\t{chr_name}\t{start}\t{end}\t{LOH_flag}\t{amplification}\t{four}\t{segment[4]}\n")
            else:
                #seg_file.write('{}\tneut\t2\n'.format('\t'.join([str(x) for x in row]))) 
                seg_file.write(f"{ID}\t{chr_name}\t{start}\t{end}\t{LOH_flag}\t{neutral}\t{two}\t{segment[4]}\n")
        elif isinstance(segment[3], str):
            continue
        else:
            if float(segment[3]) > rem_input['args'].zscore: #and float(segment[3]) < rem_input['args'].zscore * 5:
                #seg_file.write('{}\tgain\t3\n'.format('\t'.join([str(x) for x in row])))
                seg_file.write(f"{ID}\t{chr_name}\t{start}\t{end}\t{LOH_flag}\t{gain}\t{three}\t{segment[4]}\n")
            elif float(segment[3]) < - rem_input['args'].zscore:
                #seg_file.write('{}\tloss\t1\n'.format('\t'.join([str(x) for x in row])))
                seg_file.write(f"{ID}\t{chr_name}\t{start}\t{end}\t{LOH_flag}\t{loss}\t{one}\t{segment[4]}\n")
            #elif float(segment[3]) > rem_input['args'].zscore * 5:
            #    seg_file.write('{}\tamp\n'.format('\t'.join([str(x) for x in row])))
            else:
               #seg_file.write('{}\tneut\t2\n'.format('\t'.join([str(x) for x in row]))) 
               seg_file.write(f"{ID}\t{chr_name}\t{start}\t{end}\t{LOH_flag}\t{neutral}\t{two}\t{segment[4]}\n")

    seg_file.close()

def _generate_chr_statistics_file(rem_input, results):
    stats_file = open('{}_statistics.txt'.format(rem_input['args'].outid), 'w')
    stats_file.write('chr\tratio.mean\tratio.median\tzscore\n')
    chr_ratio_means = [np.ma.average(results['results_r'][chr], weights=results['results_w'][chr])
                       for chr in range(len(results['results_r']))]
    chr_ratio_medians = [np.median([x for x in results['results_r'][chr] if x != 0])
                         for chr in range(len(results['results_r']))]

    results_c_chr = [[x, 0, rem_input['bins_per_chr'][x] - 1, chr_ratio_means[x]]
                     for x in range(len(results['results_r']))]

    msv = round(get_median_segment_variance(results['results_c'], results['results_r']), 5)
    cpa = round(get_cpa(results['results_c'], rem_input['binsize']), 5)
    chr_z_scores = get_z_score(results_c_chr, results)

    for chr in range(len(results['results_r'])):

        chr_name = str(chr + 1)
        if chr_name == '23':
            chr_name = 'X'
        if chr_name == '24':
            chr_name = 'Y'

        row = [chr_name,
               chr_ratio_means[chr],
               chr_ratio_medians[chr],
               chr_z_scores[chr]]

        stats_file.write('\t'.join([str(x) for x in row]) + '\n')

    stats_file.write('Gender based on --yfrac (or manually overridden by --gender): {}\n'
                     .format(str(rem_input['actual_gender'])))

    stats_file.write('Number of reads: {}\n'
                     .format(str(rem_input['n_reads'])))

    stats_file.write('Standard deviation of the ratios per chromosome: {}\n'
                     .format(str(round(float(np.nanstd(chr_ratio_means)), 5))))

    stats_file.write('Median segment variance per bin (doi: 10.1093/nar/gky1263): {}\n'
                     .format(str(msv)))

    stats_file.write('Copy number profile abnormality (CPA) score (doi: 10.1186/s13073-020-00735-4): {}\n'
                     .format(str(cpa)))

    stats_file.close()
