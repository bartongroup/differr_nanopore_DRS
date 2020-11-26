import numpy as np


MISMATCH_COUNTS_COLUMNS = ['base', 'kd_mm', 'kd_m', 'cntrl_mm', 'cntrl_m']


def get_reference_base(fasta, chrom, pos):
    '''
    Get the reference base for a position, can be used to
    generate odds ratio of mismatch/match change in kd vs cntrl
    '''
    ref_base = fasta.fetch(chrom, pos, pos + 1)
    return ref_base


def calculate_mismatch_odds_ratio(ref_base, kd_counts, cntrl_counts):
    '''
    Calculate the log2 odds ratio with haldane correction
    (0.5 pseudocount)
    '''
    try:
        kd_m = float(kd_counts[ref_base])
    except KeyError:
        kd_m = 0.0
    kd_mm = float(kd_counts.sum() - kd_m)

    try:
        cntrl_m = float(cntrl_counts[ref_base])
    except KeyError:
        cntrl_m = 0.0
    cntrl_mm = float(cntrl_counts.sum() - cntrl_m)
    odds_ratio = np.log2(((cntrl_mm + 0.5) / (cntrl_m + 0.5)) / 
                         ((kd_mm + 0.5) / (kd_m + 0.5)))
    return odds_ratio, [ref_base, kd_mm, kd_m, cntrl_mm, cntrl_m]