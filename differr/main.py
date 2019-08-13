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


BATCH_SIZE = 1_000_000


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


def run_differr_analysis(kd_bam_fns, cntrl_bam_fns, fasta_fn,
                         res_hdf5_fn=None,
                         batch_size=1_000_000,
                         median_expr_threshold=10,
                         min_expr_threshold=0,
                         fdr_threshold=0.05,
                         processes=6,
                         max_depth=10_000_000,
                         normalise=True):
    '''
    run the whole analysis on a set of bam files. Returns a
    dataframe filtered by fdr.
    '''
    norm_factors = get_norm_factors(chain(kd_bam_fns, cntrl_bam_fns), normalise)
    references = get_references_and_lengths(cntrl_bam_fns[0])
    results = []
    if res_hdf5_fn is not None:
        res_hdf5 = ResultsHDF5(res_hdf5_fn, norm_factors, references)
    kd_has_replicates = len(kd_bam_fns) > 1
    cntrl_has_replicates = len(cntrl_bam_fns) > 1
    with Parallel(n_jobs=processes) as pool, pysam.FastaFile(fasta_fn) as fasta:
        for ref_name, ref_len in references.items():
            for i in range(0, ref_len, batch_size):
                query = (ref_name, i, min(ref_len, i + batch_size))
                # use joblib to process queries for all the bam files in parallel.
                # means we have to reopen the bam file each time (as they aren't
                # pickle-able) but its worth it for the speedup if batch_size is
                # big enough, probably.
                counts = pool(
                    delayed(process_bam_chunk)(
                        bam_fn, query, norm_factors[bam_fn], max_depth)
                    for bam_fn in chain(kd_bam_fns, cntrl_bam_fns)
                )
                # convert the list of dataframes to one big dataframe with multiindex
                # columns - one level for sample number and another for condition
                # luckily joblib.Parallel preserves the order of the input :)
                counts = {
                    'kd': pd.concat(counts[:len(kd_bam_fns)],
                                    axis=1, keys=range(len(kd_bam_fns))),
                    'cntrl': pd.concat(counts[len(kd_bam_fns):],
                                       axis=1, keys=range(len(cntrl_bam_fns))),
                }
                counts = pd.concat(counts, axis=1, join='inner').fillna(0)
                counts = filter_by_expression(counts, median_expr_threshold, min_expr_threshold)
                stats, mismatch_counts = test_significant_der_sites(counts, fdr_threshold, fasta)
                results.append(stats)
                if res_hdf5_fn is not None:
                    res_hdf5.write(counts, mismatch_counts, stats)
                click.echo('Chrom {}: processed {:.1f} Megabases'.format(ref_name, i / 1_000_000))
    results = pd.concat(results, axis=0)
    _, results['hetero_G_fdr'], *_ = multipletests(results.hetero_G_pval.fillna(1), method='fdr_bh')
    # Filter by FDR threshold:
    results = results[results.hetero_G_fdr < fdr_threshold]
    # Filter results where the sum of homogeneity G statistics is greater
    # than G statistic for cntrl vs kd
    # only if we have replicates in both conds
    if kd_has_replicates and cntrl_has_replicates:
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
@click.option('-p', '--processes', default=-1)
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
def differr(cond_a_bams, cond_b_bams,
            reference_fasta, output_bed,
            raw_counts_hdf,
            fdr_threshold, processes,
            max_depth, normalise,
            median_expr_threshold,
            min_expr_threshold):
    '''
    A script for detecting differential error rates in aligned Nanopore data
    '''
    if processes == -1:
        processes = min(len(cond_a_bams) + len(cond_b_bams), cpu_count())
    results = run_differr_analysis(
        cond_a_bams,
        cond_b_bams,
        reference_fasta,
        raw_counts_hdf,
        batch_size=BATCH_SIZE,
        median_expr_threshold=median_expr_threshold,
        min_expr_threshold=min_expr_threshold,
        fdr_threshold=fdr_threshold,
        processes=processes,
        max_depth=max_depth,
        normalise=normalise,
    )
    write_stats_to_bed(output_bed, results)


if __name__ == '__main__':
    differr()