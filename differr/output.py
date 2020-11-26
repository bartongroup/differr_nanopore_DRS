import pandas as pd
import numpy as np


class ResultsHDF5(object):
    '''
    Wrapper for writing all results (not just the significant ones in the bed file) to HDF5
    '''

    def __init__(self, h5_fn, norm_factors, references):
        self.output_h5 = pd.HDFStore(h5_fn, mode='w')
        self.output_h5.put('norm_factors', pd.Series(norm_factors))
        self.pos_data_minsizes = {
            'chrom': max([len(x) for x in references]),
            'strand': 1,
            'base': 1
        }

    def write(self, counts, mismatches, statistics):
        assert len(counts) == len(mismatches) == len(statistics)
        # extract positional data from statistics df
        pos = statistics.iloc[:, :3]
        pos['base'] = mismatches.pop('base')
        statistics = statistics.iloc[:, 3:]
        self.output_h5.append('positional_data',
                              pos,
                              min_itemsize=self.pos_data_minsizes)
        self.output_h5.append('counts', counts.reset_index(drop=True))
        self.output_h5.append('mismatch_counts', mismatches)
        self.output_h5.append('statistics', statistics)

    def close(self):
        self.output_h5.close()
        

def nlog10(val):
    return np.negative(np.log10(val))


BED_RECORD = (
    '{chrom}\t{start:d}\t{end:d}\t'
    'DER_site\t{score:.0f}\t{strand}\t'
    '{odds_ratio:.2f}\t'
    '{G_stat:.2f}\t{p_val:.2f}\t{q_val:.2f}\t'
    '{G_stat_A:.2f}\t{G_stat_A_p_val:.2f}\t'
    '{G_stat_B:.2f}\t{G_stat_B_p_val:.2f}\n'
)


def write_stats_to_bed(bed_fn, results):
    '''
    Write the results to a bed file (with a bunch of extra columns for G-test results)
    '''
    with open(bed_fn, 'w') as bed:
        for _, chrom, pos, strand, odds_ratio, *G_test_res in results.itertuples():
            G_A, G_A_p, G_B, G_B_p, G, G_p, G_fdr = G_test_res
            bed.write(BED_RECORD.format(
                chrom=chrom,
                start=pos,
                end=pos + 1,
                score=nlog10(G_fdr),
                strand=strand,
                odds_ratio=odds_ratio,
                G_stat=G,
                p_val=nlog10(G_p),
                q_val=nlog10(G_fdr),
                G_stat_A=G_A,
                G_stat_A_p_val=G_A_p,
                G_stat_B=G_B,
                G_stat_B_p_val=G_B_p
            ))