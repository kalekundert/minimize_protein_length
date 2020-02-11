#!/usr/bin/env python3

import numpy as np
import pandas as pd
import logging
import autoprop, toml

from Bio.pairwise2 import align, format_alignment
from Bio.SubsMat.MatrixInfo import blosum80
from pyrosetta import init as init_rosetta, pose_from_file
from pyrosetta.rosetta.protocols import loophash
from pyrosetta.rosetta.utility import fixedsizearray1_double_6_t as Real6

class LoophashDatabases:

    @classmethod
    def load(cls, db_dir):
        # Only load the small databases, because we're looking for turns.
        lhms = {}
        for i in range(3, 14+1):
            path = db_dir / f'loopdb.{i}.db'
            lhms[i] = loophash.LoopHashMap(i)
            lhms[i].read_db(str(path))

        bbs = loophash.BackboneDB()
        bbs.read_db(str(db_dir / 'backbone.db'), True)

        return cls(lhms, bbs)

    def __init__(self, lhms, bbs):
        self.lhms = lhms
        self.bbs = bbs

class AlignedPose:

    @classmethod
    def from_file(cls, msa, path, chain):
        pose = pose_from_file(str(path))
        return cls.from_pose(msa, pose, chain)

    @classmethod
    def from_pose(cls, msa, pose, chain):
        mapped_indices = map_pose_indices(
                msa.ref_ungapped,
                pose,
                chain,
        )
        return cls(msa.ref_ungapped, pose, mapped_indices)


    def __init__(self, ref_ungapped, pose, mapped_indices):
        self.pose = pose
        self.resis = mapped_indices
        self.ref_ungapped = ref_ungapped

    def has_resis(self, *indices):
        return all(i in self.resis for i in indices)

@autoprop
class GapFilters:

    def to_dict(self):
        return dict(
                n_residues=self.n_residues,
                n_missing=self.n_missing,
                filters=self._counts,
                filter_order=self._order
        )

    @classmethod
    def from_dict(cls, d):
        filters = cls(d['n_residues'], d['n_missing'])
        filters._counts = d['filters']
        filters._order = d['filter_order']
        return filters


    def to_toml(self, path):
        with open(path, 'w') as f:
            toml.dump(self.to_dict(), f)

    @classmethod
    def from_toml(cls, path):
        return cls.from_dict(toml.load(path))


    def __init__(self, n_residues, n_missing):
        self.n_residues = n_residues
        self.n_missing = n_missing
        self._counts = {}
        self._order = {}

    def __str__(self):
        s = ""
        extra_info = {
                "All possible residue pairs":
                    f"({self.n_residues} residues)",
                "Residue missing from pose":
                    f"({self.n_missing} missing residues; expect −{self.n_missing_combos})",
        }

        def format(name, remaining, count=None, indent=0):
            cols = [
                    '{0:{1}}'.format(x, fmt)
                    for x, fmt in [
                        (f'{"→ " if indent else ""}{name}: ',       '<35'),
                        (remaining,                                 '>6'),
                        (f'−{count}' if count else '',              '>7'),
                        (extra_info.get(name, ''),                  's'),
                    ]
            ]
            return '  '.join(cols).strip() + '\n'

        s += format("All possible residue pairs", self.n_combos)
        n = self.n_combos
        for filter, count in self:
            n -= count
            s += format(filter, n, count, indent=2)
        return s

    def __len__(self):
        return len(self._counts)

    def __iter__(self):
        from operator import itemgetter
        yield from sorted(
                self._counts.items(),
                key=lambda x: self._order[x[0]],
        )

    def __getitem__(self, key):
        return self._counts[key]

    def __setitem__(self, key, value):
        self._counts[key] = value
        if key not in self._order:
            self._order[key] = len(self._order)

    def get_n_combos(self):
        return self.n_residues * (self.n_residues - 1) // 2

    def get_n_missing_combos(self):
        return self.n_missing * (self.n_residues - self.n_missing) + \
               self.n_missing * (self.n_missing - 1) // 2

    def get_n_accepted(self):
        return self.n_combos - sum(self._counts.values())

def choose_gaps_to_delete(scores, gaps, filters):
    hits = []
    score_threshold = np.mean(scores)
    n_res = len(scores)

    max_loop = 6
    max_deletion = 0.1

    filters[f'Spanning loop >{max_loop} residues'] = 0
    filters[f'Deleted >{100*max_deletion:.0f}% of the protein'] = 0
    filters['Below-average deletion score'] = 0

    for i, row in gaps.iterrows():
        if row['len_spanning_loop'] > max_loop:
            filters[f'Spanning loop >{max_loop} residues'] += 1
            continue

        if row['len_deletion'] / n_res > max_deletion:
            filters[f'Deleted >{100*max_deletion:.0f}% of the protein'] += 1
            continue

        del_start, del_end, del_score = pick_deletion_window(
                scores,
                row['gap_start'],
                row['gap_end'],
                row['len_deletion'],
        )
        if del_score <= score_threshold:
            filters['Below-average deletion score'] += 1
            continue

        hits += [{
            'del_start'         : del_start,
            'del_end'           : del_end,
            'del_score'         : del_score,
            **row.to_dict(),
        }]

    df = pd.DataFrame(hits)
    df_uniq = df.\
            groupby(['del_start', 'del_end'], as_index=False).\
            agg(np.mean)

    filters['Duplicates'] = len(df) - len(df_uniq)

    return df_uniq

def load_spannable_gaps(work_dels, msa):
    if work_dels.spannable_hdf5.exists():
        gaps = pd.read_hdf(work_dels.spannable_hdf5)
        filters = GapFilters.from_toml(work_dels.filters_toml)

    else:
        logging.basicConfig(
                filename=work_dels.rosetta_log,
                filemode='w',
        )

        init_rosetta(set_logging_handler='logging')
        db = LoophashDatabases.load(work_dels.shared.loophash_db)
        aligned_pose = AlignedPose.from_file(
                msa,
                work_dels.shared.target_pdb,
                work_dels.shared.target_pdb_chain,
        )

        gaps, filters = find_spannable_gaps(db, aligned_pose)

        with work_dels.touch(work_dels.spannable_hdf5) as p:
            gaps.to_hdf(p, 'spannable_gaps')
        with work_dels.touch(work_dels.filters_toml) as p:
            filters.to_toml(p)

    return gaps, filters

def find_spannable_gaps(db, aligned_pose):
    from itertools import combinations
    from tqdm import tqdm

    hits = []
    n_res = len(aligned_pose.ref_ungapped)
    n_missing = n_res - len(aligned_pose.resis)

    filters = GapFilters(n_res, n_missing)
    filters['Residue missing from pose'] = 0
    filters['No spanning loop found'] = 0
    filters['Gap smaller than spanning loop'] = 0

    all_pairs = tqdm(
            combinations(range(n_res), 2),
            total=filters.n_combos,
            unit=' pair'
    )
    for i, j in all_pairs:

        if not aligned_pose.has_resis(i, j):
            filters['Residue missing from pose'] += 1
            continue

        # Leaps involving the first residue of the pose will trigger a segfault 
        # in `get_rt_over_leap_without_foldtree_bs()`, due to it attempting to 
        # access an atom from the previous residue to complete the stub.  The 
        # other `get_rt_over_leap*()` functions handle this case more 
        # gracefully (returning False and reporting an error), but still fail 
        # to calculate the leap.  
        #
        # The workaround to skip the first residue, when it comes up.  At first 
        # I tried to modify the iterator to exclude the first position in the 
        # reference sequence, but that position may or may not correspond to 
        # the first position in the pose (e.g. if the pose in missing some 
        # N-terminal residues).  So the check has to be done here.

        if aligned_pose.resis[i] == 1:
            filters['No spanning loop found'] += 1
            continue

        n_loop_res = find_smallest_spanning_loop(
                db, 
                aligned_pose.pose,
                aligned_pose.resis[i],
                aligned_pose.resis[j],
        )
        if n_loop_res == 0:
            filters['No spanning loop found'] += 1
            continue

        n_gap = j - i
        n_del = n_gap - n_loop_res
        if n_del < 1:
            filters['Gap smaller than spanning loop'] += 1
            continue

        hits += [{
            'gap_start': i,
            'gap_end': j,
            'len_gap': n_gap,
            'len_spanning_loop': n_loop_res,
            'len_deletion': n_del,
        }]

    return pd.DataFrame(hits), filters

def find_smallest_spanning_loop(db, pose, resi_n, resi_c, radius=1):
    """
    Return the length of the smallest loop (in the given loophash databases) 
    that can span the gap between the given residues in the given pose.  Return 
    0 if no such loops could be found.
    """
    # See lab notebook and `check_loophash_recovery`:
    # radius=1 seems to be the best for recovering known loop lengths.  Since 
    # I'm trying to use loophash as a ruler here, that is a good metric for 
    # what I want.

    # Find the 6D transformation between the ends of the gap.
    rt = Real6();
    success = loophash.get_rt_over_leap_without_foldtree_bs(
            pose, resi_n, resi_c, rt
    )
    if not success:
        raise ValueError(f"Cannot determine 6D transformation between residues {resi_n} and {resi_c} (pose-numbering).")

    # Find the shortest loop that can bridge the gap.
    for k, lhm in sorted(db.lhms.items()):
        if lhm.radial_count(radius, rt) > 0:
            return k

    return 0

def pick_deletion_window(scores, i, j, n_del):
    # Note that `i` and `j` refer to individual residues, while `del_start` and 
    # `del_end` refer to pseudo-indices between residues, e.g. for use with 
    # slicing.
    del_scores = np.convolve(scores[i:j+1], np.ones(n_del), 'valid') / n_del
    i_max = np.argmax(del_scores)
    del_start = i + i_max
    del_end = del_start + n_del
    del_score = del_scores[i_max]
    return del_start, del_end, del_score

def map_pose_indices(ref_ungapped, pose, chain):
    pose_seq, pose_indices = get_chain_info(pose, chain)

    # Align the pose sequence with the (ungapped) reference sequence:
    #
    # I decided to use BLOSUM80 because the two sequences in this case should 
    # be very similar (not that it should matter much).  I used the 
    # corresponding gap penalties from BLAST 2.2.27, which I found in the 
    # reference below:
    #
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3848038/

    alignments = align.globalds(
            pose_seq,
            ref_ungapped,
            blosum80,
            -10, -1,
    )
    aligned_pose, aligned_ref, score, start, end = alignments[0]

    i_ref = i_pose = 0
    mapped_indices = {}

    for aa_ref, aa_pose in zip(aligned_ref, aligned_pose):
        hit_ref = aa_ref not in '-.'
        hit_pose = aa_pose not in '-.'

        if hit_ref and hit_pose:
            mapped_indices[i_ref] = pose_indices[i_pose]

        i_ref += hit_ref
        i_pose += hit_pose

    return mapped_indices

def get_chain_info(pose, chain):
    """
    Return the sequence and the pose-numbered residue indices for the given 
    chain, including only residues that are actually present in the pose.
    """
    seq = ""
    resis = []

    for i in range(pose.conformation().size()):
        i += 1  # Rosetta is 1-indexed.
        aa = pose.residue_type(i).name1()
        aa_chain = pose.pdb_info().chain(i)

        if aa_chain == chain and aa not in 'JOBZUX':
            seq += aa
            resis.append(i)

    return seq, resis


