#!/usr/bin/env python3

from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, Gapped
from Bio.Align import MultipleSeqAlignment
from Bio.SubsMat import MatrixInfo

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

def load_msa(msa_path, msa_format, ref_id):
    """
    The given multiple sequence alignment (MSA) should contain the reference 
    sequence.  The reference will be removed from the alignment and returned 
    separately.  The alignment will also be changed such that "." is used to 
    indicate terminal deletions while "-" is used to indicate internal 
    deletions.
    """
    msa_with_ref = AlignIO.read(msa_path, msa_format)

    ref = None
    msa = MultipleSeqAlignment([])

    for record in msa_with_ref:
        # Use "." to indicate terminal mismatches, and "-" to indicate internal 
        # mismatches.

        to_dots = lambda m: '.' * (m.end() - m.start())
        record.seq = Seq(
                re.sub('^-+|-+$', to_dots, str(record.seq)),
                record.seq.alphabet,
        )

        if record.id == ref_id:
            ref = record
        else:
            msa.append(record)

    return ref, msa

def weight_alignments(ref, msa):
    """
    Weight each alignment by percent identity to the reference.
    """
    def percent_identity(a, b):
        assert len(a) == len(b)
        N = len(b.strip('.'))
        return sum(aa == bb and aa not in '-.' for aa, bb in zip(a, b)) / N

    for record in msa:
        record.weight = percent_identity(ref.seq, record.seq)

    msa.sort(key=lambda x: x.weight, reverse=True)


def choose_deletions(ref, msa):
    scores = calc_deletion_scores(ref, msa)
    dfs = []

    for threshold in scores:
        runs = find_runs(scores >= threshold)
        df = pd.DataFrame(runs, columns=['del_start', 'del_end'])
        dfs.append(df)

    dels = pd.concat(dfs)\
            .drop_duplicates()\
            .reset_index(drop=True)
    
    def score_dels(x):
        return np.mean(scores[ x['del_start'] : x['del_end'] ])

    dels['score'] = dels.apply(score_dels, axis=1)
    dels['del_len'] = dels['del_end'] - dels['del_start']

    dels.sort_values(by='score', ascending=False, inplace=True)

    # Require that each deletion score better than the average over the whole 
    # sequence.  This also conveniently gets rid of the 0-N deletion.
    dels = dels[ dels['score'] > np.mean(scores) ]

    return dels.reset_index(drop=True)

def output_deletions(ref, dels, dels_path, dels_format='fasta'):
    unalign = str.maketrans('', '', '-.')
    wt = str(ref.seq).translate(unalign)

    records = []
    aligned_records = []

    for i, row in dels.iterrows():
        # Because all the columns in the data frame are numeric, pandas 
        # automatically coerces the ints to floats.  I think this is a bug, but 
        # it's not too hard to work around.  See pandas #28308.
        start, end = int(row['del_start']), int(row['del_end'])
        start1 = start + 1

        id = f'{ref.id}Δ{start1}' + (f'-{end}' if start1 != end else '')
        desc = f'score={row["score"]}'
        seq = Seq(
            wt[:start] + wt[end:],
            IUPAC.protein,
        )
        aligned_seq = Seq(
            wt[:start] + ('-' * (end - start)) + wt[end:],
            Gapped(IUPAC.protein, '-'),
        )

        record = SeqRecord(seq, id=id, description=desc)
        records.append(record)

        aligned_record = SeqRecord(aligned_seq, id=id, description=desc)
        aligned_records.append(aligned_record)

    SeqIO.write(records, dels_path, dels_format)

    msa = MultipleSeqAlignment(aligned_records)
    msa_path = f'{Path(dels_path).stem}.aln'
    print(msa_path)
    AlignIO.write(msa, msa_path, 'clustal')

def find_runs(a):
    # Taken from:
    # https://stackoverflow.com/questions/1066758/find-length-of-sequences-of-identical-values-in-a-numpy-array-run-length-encodi

    # Make sure each run both starts and stops.
    bounded = np.hstack(([0], a, [0]))

    # This will give 1 at run starts and -1 at run ends.
    edges = np.diff(bounded)

    # The start and end indices refer to the positions before the runs starts 
    # and ends, respectively, which would be too small by 1.  But this error is 
    # inadvertently corrected by the fact that we added an element to the 
    # beginning of the input array, effectively making every index bigger by 1.  
    # So in the end, no further corrections are necessary.
    run_starts, = np.where(edges > 0)
    run_ends, = np.where(edges < 0)

    return zip(run_starts, run_ends)


def calc_deletion_scores(ref, msa):
    """
    Sum the weight of each deletion at each position in the reference 
    sequence.
    """
    i = map_msa_indices_to_ref(ref, msa)
    deletion_scores = np.zeros(len(i))

    for record in msa:
        for msa_i in i:
            deletion_scores[i[msa_i]] += \
                    record.weight * (record.seq[msa_i] == '-')

    return deletion_scores

def calc_conservation_scores(ref, msa):
    i = map_msa_indices_to_ref(ref, msa)
    conservation_scores = np.zeros(len(i))
    subs_mat = MatrixInfo.blosum80

    def apply_subs_mat(a, b):
        if (a,b) in subs_mat:
            return subs_mat[a,b]
        if (b,a) in subs_mat:
            return subs_mat[b,a]
        return 0

    for msa_i in i:
        for record in msa:
            subs_score = apply_subs_mat(
                    ref.seq[msa_i],
                    record.seq[msa_i],
            )

            conservation_scores[i[msa_i]] += record.weight * subs_score

    return conservation_scores

def map_msa_indices_to_ref(ref, msa):
    ref_i_from_msa_j = {}

    ref_i = 0
    for msa_j, seq in enumerate(ref.seq):
        if seq == '-': continue
        ref_i_from_msa_j[msa_j] = ref_i
        ref_i += 1

    return ref_i_from_msa_j


if __name__ == '__main__':
    ref, msa = load_msa('blast_refseq_100.aln', 'clustal', 'WP_032461047.1')
    weight_alignments(ref, msa)

    # Calculate and record the deletion scores.
    deletion_scores = calc_deletion_scores(ref, msa)
    conservation_scores = calc_conservation_scores(ref, msa)

    deletion_scores.tofile('deletion_scores.dat')
    conservation_scores.tofile('conservation_scores.dat')

    #i = np.arange(len(deletion_scores))

    #conservation_scores = pd.Series(conservation_scores)
    #plt.plot(i, conservation_scores.rolling(10).mean())

    #plt.plot(i, conservation_scores)
    #plt.plot(i, deletion_scores)
    #plt.show()

    # Output a list of candidate deletions.
    dels = choose_deletions(ref, msa)
    output_deletions(ref, dels, 'blast_refseq_100_threshold_dels.fasta')


def test_find_runs():
    examples = [
            ([], []),

            ([0], []),
            ([1], [(0, 1)]),

            ([0, 0], []),
            ([0, 1], [(1, 2)]),
            ([1, 0], [(0, 1)]),
            ([1, 1], [(0, 2)]),

            ([0, 0, 0, 0], []),
            ([1, 1, 1, 1], [(0, 4)]),
            ([1, 1, 0, 0], [(0, 2)]),
            ([0, 1, 1, 0], [(1, 3)]),
            ([0, 0, 1, 1], [(2, 4)]),
            ([1, 0, 0, 1], [(0, 1), (3, 4)]),
    ]

    for example, result in examples:
        assert list(find_runs(example)) == result

