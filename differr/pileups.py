from collections import Counter
import pandas as pd
import pysam


COLUMNS = ['A', 'C', 'G', 'T', 'indel']


def calculate_norm_factor(bam):
    '''
    Get the sequencing depth normalisation factor
    to apply to read counts to make them more
    comparable between libraries
    '''
    return bam.mapped / 1_000_000


def get_norm_factors(bam_list, normalise=True):
    '''
    Open all the bam files in a list and calculate a normalisation
    factor for the file - this is generally just the number of mapped
    reads divided by a million
    '''
    norm_factors = {}
    if normalise:
        for bam_fn in bam_list:
            with pysam.AlignmentFile(bam_fn) as bam:
                norm_factors[bam_fn] = calculate_norm_factor(bam)
    else:
        norm_factors = {bam_fn: 1 for bam_fn in bam_list}
    return norm_factors


def get_references_and_lengths(bam_fn):
    '''Open a bam file and get the length of all the references from the header'''
    with pysam.AlignmentFile(bam_fn) as bam:
        reference_lengths = {ref: bam.get_reference_length(ref) for ref in bam.references}
    return reference_lengths


def get_reference_pos(pileupcol):
    '''
    Get relationship to reference, i.e. chrom and pos, from bam
    (we don't know what the reference base actually is though)
    '''
    ref_name = pileupcol.reference_name
    ref_pos = pileupcol.reference_pos
    return ref_name, ref_pos


def get_single_query_seq(pileupread):
    '''
    Get the query sequence for a pysam PileUpRead
    '''
    pos = pileupread.query_position
    if pos is None:
        return ''
    else:
        seq = pileupread.alignment.query_sequence[pos]
        if len(seq) > 1:
            return ''
        else:
            return seq


def get_query_seqs(pileupcol):
    # Use pysam to get the query sequences, or if this fails fall
    # back on slower but more stable method of fetching one at a time
    try:
        query_seqs = pileupcol.get_query_sequences()
    except AssertionError:
        query_seqs = [get_single_query_seq(read) for read in pileupcol.pileups]
    query_seqs = [s.upper() if s else 'indel' for s in query_seqs]
    return query_seqs


def iter_pileupreads(pileupcol):
    '''
    Generator yielding strand and query sequence for reads in pileupcol
    Reads which are refskip (i.e. intronic) are skipped
    '''
    query_seqs = get_query_seqs(pileupcol)
    for read, seq in zip(pileupcol.pileups, query_seqs):
        if not read.is_refskip:
            is_reverse = read.alignment.is_reverse
            yield is_reverse, seq
    

def count_mismatches(col):
    '''
    For a pysam PileupColumn, computes the number of reads
    which contain each base or indels and returns a df
    '''
    ref_name, ref_pos = get_reference_pos(col)
    refs = [(ref_name, ref_pos, '+'),
            (ref_name, ref_pos, '-')]
    counts = {r: Counter() for r in refs}
    for is_reverse, query_seq_type in iter_pileupreads(col):
            strand = refs[is_reverse]
            counts[strand][query_seq_type] += 1
    return counts


def process_bam_chunk(bam_fn, query, norm_factor=1, max_depth=10_000_000):
    '''
    Open a bam file and get the base counts for a query.
    Designed to be called in parallel on a number of
    bams at once (hence why we have to reopen the bam
    for every query chunk we process).
    '''
    with pysam.AlignmentFile(bam_fn) as bam:
        chunk_res = {}
        for col in bam.pileup(
                *query, max_depth=max_depth,
                truncate=True, min_base_quality=0):

            chunk_res.update(count_mismatches(col))
        chunk_res = pd.DataFrame.from_dict(chunk_res,
                                           orient='index',
                                           columns=COLUMNS)
    return chunk_res.fillna(0) / norm_factor