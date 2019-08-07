import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency

from .mismatch import (
    get_reference_base,
    calculate_mismatch_odds_ratio,
    MISMATCH_COUNTS_COLUMNS,
)   


RESULT_COLUMNS = [
    'chrom', 'pos', 'strand', 'mismatch_odds_ratio',
    'homog_G_cntrl', 'homog_G_cntrl_pval',
    'homog_G_kd', 'homog_G_kd_pval',
    'hetero_G', 'hetero_G_pval' # FDR is calculated later
]

def filter_empty_cols(c):
    '''Remove columns for bases which are not represented in the data'''
    c_t = c.T
    c_t = c_t[c_t.any(1)]
    return c_t.T


def test_significant_der_sites(counts, fdr_threshold, fasta):
    '''Conduct G-test on counts to identify sites with differential error profiles'''
    res = []
    mismatch_counts = []
    for (chrom, pos, strand), c in counts.iterrows():
        c = filter_empty_cols(c.unstack(-1))
        cntrl = c.loc['cntrl']
        kd = c.loc['kd']
        # pool replicates and conduct G test of kd vs cntrl
        cntrl_pooled = cntrl.sum(0)
        kd_pooled = kd.sum(0)
        kd_vs_cntrl_g, kd_vs_cntrl_p_val, *_ = chi2_contingency(
            [cntrl_pooled, kd_pooled], lambda_='log-likelihood')
        # if result is likely to be significant we need to test
        # homogeneity of cntrl and kd replicates.
        # If it isn't we can save some time by skipping these tests!
        # We haven't calculated FDR yet but it cannot be smaller than P
        if kd_vs_cntrl_p_val < fdr_threshold:
            # homogeneity test of cntrl
            cntrl_hom_g, cntrl_hom_p_val, *_ = chi2_contingency(filter_empty_cols(cntrl),
                                                               lambda_='log-likelihood')
            # homogeneity test of kd
            kd_hom_g, kd_hom_p_val, *_ = chi2_contingency(filter_empty_cols(kd),
                                                        lambda_='log-likelihood')
        else:
            # Otherwise set NaNs
            cntrl_hom_g, cntrl_hom_p_val = np.nan, np.nan
            kd_hom_g, kd_hom_p_val = np.nan, np.nan

        # use reference base to generate odds ratio for mismatch compared to the reference
        # (kd_mm / kd_m) / (cntrl_m / cntrl_mm)
        ref_base = get_reference_base(fasta, chrom, pos)
        odds_ratio, mm_counts = calculate_mismatch_odds_ratio(
            ref_base, kd_pooled, cntrl_pooled
        )
        res.append([chrom, pos, strand, odds_ratio,
                    cntrl_hom_g, cntrl_hom_p_val,
                    kd_hom_g, kd_hom_p_val,
                    kd_vs_cntrl_g, kd_vs_cntrl_p_val])
        mismatch_counts.append(mm_counts)
    res = pd.DataFrame(res, columns=RESULT_COLUMNS)
    mismatch_counts = pd.DataFrame(mismatch_counts, columns=MISMATCH_COUNTS_COLUMNS)
    return res, mismatch_counts