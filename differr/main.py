from itertools import chain
from multiprocessing import cpu_count

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from joblib import Parallel, delayed
import pysam
import click

from .pileups import (
    get_norm_factors,
    get_references_and_lengths,
    process_bam_chunk
)
from .statistics import test_significant_der_sites
from .output import ResultsHDF5, write_stats_to_bed



def filter_by_expression(counts, median_threshold, min_threshold):
    '''
    Filter the counts matrix by median expression of replicates with
    a condition, and a minimum expression level across all samples
    '''
    if counts.empty:
        return counts
    # If all reps are 100% matches we don't need to do any tests
    is_variable = counts.groupby(level=2, axis=1).any().sum(1)
    var_filt = is_variable > 1
    per_rep_expr = counts.groupby(level=(0, 1), axis=1).sum()
    per_cond_median_expr = per_rep_expr.groupby(level=0, axis=1).median()
    # The median expression in both conditions must be greater than this value
    med_filt = (per_cond_median_expr >= median_threshold).all(1)
    # The expression in all replicates of all conditions must be greater than this value
    # Currently there must just be at least one read in all conditions
    # i.e. threshold is zero
    min_filt = (per_rep_expr >= min_threshold).all(1)
    return counts[var_filt & med_filt & min_filt]


def _parallel_differr(kd_bam_fns, cntrl_bam_fns,
                      fasta_fn, query,
                      norm_factors,
                      median_expr_threshold=10,
                      min_expr_threshold=0,
                      fdr_threshold=0.05,
                      max_depth=10_000_000):
    kd_counts = [
        process_bam_chunk(bam_fn, query, norm_factors[bam_fn], max_depth)
        for bam_fn in kd_bam_fns
    ]
    cntrl_counts = [
        process_bam_chunk(bam_fn, query, norm_factors[bam_fn], max_depth)
        for bam_fn in cntrl_bam_fns
    ]
    counts = {
        'kd': pd.concat(kd_counts, axis=1, keys=range(len(kd_bam_fns))),
        'cntrl': pd.concat(cntrl_counts, axis=1, keys=range(len(cntrl_bam_fns))),
    }
    counts = pd.concat(counts, axis=1, join='inner').fillna(0)
    counts = filter_by_expression(counts, median_expr_threshold, min_expr_threshold)
    with pysam.FastaFile(fasta_fn) as fasta:
        stats, mismatch_counts = test_significant_der_sites(counts, fdr_threshold, fasta)
    return counts, mismatch_counts, stats


def run_differr_analysis(kd_bam_fns, cntrl_bam_fns, fasta_fn,
                         res_hdf5_fn=None,
                         median_expr_threshold=10,
                         min_expr_threshold=0,
                         fdr_threshold=0.05,
                         processes=6,
                         max_depth=10_000_000,
                         normalise=True,
                         query=None):
    '''
    run the whole analysis on a set of bam files. Returns a
    dataframe filtered by fdr.
    '''
    norm_factors = get_norm_factors(chain(kd_bam_fns, cntrl_bam_fns), normalise)
    references = get_references_and_lengths(cntrl_bam_fns[0], query)
    results = []
    if res_hdf5_fn is not None:
        res_hdf5 = ResultsHDF5(res_hdf5_fn, norm_factors, references)
    parallel_args = {
        'norm_factors': norm_factors,
        'median_expr_threshold': median_expr_threshold,
        'min_expr_threshold': min_expr_threshold,
        'fdr_threshold': fdr_threshold,
        'max_depth': max_depth,
    }
    querys = []
    for ref_name, (ref_start, ref_len) in references.items():
        splits = np.floor(np.linspace(ref_start, ref_len, processes + 1)).astype(int)
        querys += [
            (ref_name, i, j)
            for i, j in zip(splits[:-1], splits[1:])
        ]
    parallel_res = Parallel(n_jobs=processes)(
        delayed(_parallel_differr)(
            kd_bam_fns, cntrl_bam_fns, fasta_fn, q,
            **parallel_args
        ) for q in querys
    )
    for counts, mismatch_counts, stats in parallel_res:
        results.append(stats)
        if res_hdf5_fn is not None:
            res_hdf5.write(counts, mismatch_counts, stats)
    results = pd.concat(results, axis=0)
    _, results['hetero_G_fdr'], *_ = multipletests(results.hetero_G_pval.fillna(1), method='fdr_bh')
    # Filter by FDR threshold:
    results = results[results.hetero_G_fdr < fdr_threshold]
    # Filter results where the sum of homogeneity G statistics is greater
    # than G statistic for cntrl vs kd
    # only if we have replicates in both conds
    if (len(kd_bam_fns) > 1) and (len(cntrl_bam_fns) > 1):
        results = results[(results.homog_G_cntrl + results.homog_G_kd) < results.hetero_G]
    if res_hdf5_fn is not None:
        res_hdf5.close()
    return results
            

@click.command()
@click.option(
    '-a', '--cond-A-bams', required=True, multiple=True,
    help=('BAM files for low modification sample. To include multiple '
          'BAMs you need to use the -a flag multiple times.')
)
@click.option(
    '-b', '--cond-B-bams', required=True, multiple=True,
    help=('BAM files for wild type/complemented line with normal '
          'modification levels')
)
@click.option('-r', '--reference-fasta', required=True)
@click.option('-o', '--output-bed', required=True)
@click.option('-c', '--raw-counts-hdf', required=False, default=None)
@click.option('-f', '--fdr-threshold', default=0.05)
@click.option('-p', '--processes', default=1)
@click.option(
    '-m', '--max-depth', default=10_000_000,
    help='Maximum depth for pysam pileup, default is a very large number (basically samples all reads)'
)
@click.option(
    '--normalise/--no-normalise',
    default=True,
    help='Whether to normalise counts by sample depth (counts per million) before performing stats.'
)
@click.option(
    '--median-expr-threshold', required=False, default=10,
    help='The minimum value for the median number of reads per condition after normalisation. Default is 10.'
)
@click.option(
    '--min-expr-threshold', required=False, default=1,
    help='The minimum number of reads that ALL replicates should have after normalisation. Default is 1.'
)
@click.option('-q', '--query', required=False, default=None)
def differr(cond_a_bams, cond_b_bams,
            reference_fasta, output_bed,
            raw_counts_hdf,
            fdr_threshold, processes,
            max_depth, normalise,
            median_expr_threshold,
            min_expr_threshold,
            query):
    '''
    A script for detecting differential error rates in aligned Nanopore data
    '''
    results = run_differr_analysis(
        cond_a_bams,
        cond_b_bams,
        reference_fasta,
        raw_counts_hdf,
        median_expr_threshold=median_expr_threshold,
        min_expr_threshold=min_expr_threshold,
        fdr_threshold=fdr_threshold,
        processes=processes,
        max_depth=max_depth,
        normalise=normalise,
        query=query,
    )
    write_stats_to_bed(output_bed, results)


if __name__ == '__main__':
    differr()