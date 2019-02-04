import numpy as np


MISMATCH_COUNTS_COLUMNS = ['base', 'kd_mm', 'kd_m', 'cntrl_mm', 'cntrl_m']
RC = str.maketrans('ACGTSWYRN', 'TGCASWRYN')


def get_reference_base(fasta, chrom, pos, strand):
    '''
    Get the reference base for a position, can be used to
    generate odds ratio of mismatch/match change in kd vs cntrl
    '''
    ref_base = fasta.fetch(chrom, pos, pos + 1)
    if strand == '-':
        ref_base = ref_base.translate(RC)
    return ref_base


def calculate_mismatch_odds_ratio(ref_base, kd_counts, cntrl_counts):
    '''
    Calculate the log2 odds ratio with haldane correction
    (0.5 pseudocount)
    '''
    try:
        kd_mm = float(kd_counts[ref_base])
    except KeyError:
        kd_mm = 0.0
    kd_m = float(kd_counts.sum() - kd_mm)

    try:
        cntrl_mm = float(cntrl_counts[ref_base])
    except KeyError:
        cntrl_mm = 0.0
    cntrl_m = float(cntrl_counts.sum() - cntrl_mm)
    odds_ratio = np.log2(((kd_mm + 0.5) / (kd_m + 0.5)) / 
                         ((cntrl_mm + 0.5) / (cntrl_m + 0.5)))
    return odds_ratio, [ref_base, kd_mm, kd_m, cntrl_mm, cntrl_m]