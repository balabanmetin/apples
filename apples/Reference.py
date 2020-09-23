import itertools
import subprocess
import tempfile
from abc import ABC, abstractmethod
from apples.distance import *
import numpy as np
from apples.fasta2dic import fasta2dic
import heapq


class Reference(ABC):
    """prot flag true if protein sequences, false for nucleotide sequences"""
    def __init__(self, ref_fp, prot_flag):
        self.refs = fasta2dic(ref_fp, prot_flag, False)  # never mask low confidence bases in the reference
        self.prot_flag = prot_flag

        if prot_flag:
            self.dist_function = scoredist
        else:
            self.dist_function = jc69
    @abstractmethod
    def get_obs_dist(self, query_seq, query_name):
        pass


class Full_reference(Reference):
    def __init__(self, ref_fp, prot_flag):
        Reference.__init__(self, ref_fp, prot_flag)

    def get_obs_dist(self, query_seq, query_tag):
        obs_dist = {query_tag: 0}
        for tagr, seqr in self.refs:
            obs_dist[tagr] = self.dist_function(query_seq, seqr)
        return obs_dist


class Reduced_reference(Reference):
    def __init__(self, ref_fp, prot_flag, tree_file, threshold, baseobs):
        Reference.__init__(self, ref_fp, prot_flag)
        self.threshold = threshold
        self.baseobs = baseobs

        tc_output_file = tempfile.NamedTemporaryFile(delete=True, mode='w+t').name
        nldef = tempfile.NamedTemporaryFile(delete=True, mode='w+t')
        s = ["TreeCluster.py", "-t", str(self.threshold * 1.2), "-i", tree_file, "-m", "max", "-o", tc_output_file]
        subprocess.call(s, stdout=nldef, stderr=nldef)

        self.representatives = []
        with open(tc_output_file, "r") as tc_output:
            tc_output.readline()
            lines = map(lambda x: x.strip().split('\t'), tc_output.readlines())
            lines_sorted = sorted(lines, key=lambda x: x[1])
            for key, group in itertools.groupby(lines_sorted, lambda x: x[1]):
                # print(key, list(group))
                if key == "-1":
                    for thing in group:
                        self.representatives.append((self.refs[thing[0]], [thing[0]]))
                else:
                    things = [g[0] for g in group]
                    self.representatives.append((self._find_representative(things), things))

    def _find_representative(self, group):

        if self.prot_flag:
            alphabet = np.array(['A', 'C', 'D', 'E', 'F',
                                 'G', 'H', 'I', 'K', 'L',
                                 'M', 'N', 'P', 'Q', 'R',
                                 'S', 'T', 'V', 'W', 'Y', '-'], dtype="S1")

        else:
            alphabet = np.array([b'A', b'C', b'G', b'T', b'-'], dtype="S1")

        lookup = {n: i for i, n in enumerate(alphabet)}

        def get_consensus(arr):
            n = arr.shape[1]
            frequency_matrix = np.zeros((len(alphabet), n))
            for dna in arr:
                for base in alphabet:
                    frequency_matrix[lookup[base]] += dna == base
            return alphabet[np.argmax(frequency_matrix, axis=0)]

        group_seqs = [self.refs[g] for g in group]
        group_mat = np.array(np.vstack(group_seqs))
        return get_consensus(group_mat)

    def get_obs_dist(self, query_seq, query_tag):
        obs_dist = {}
        obs_num = 0
        representative_dists = []
        i = 0
        for consensus_seq, group in self.representatives:
            dist = self.dist_function(query_seq, consensus_seq)
            representative_dists.append((dist, i))
            i += 1
        heapq.heapify(representative_dists)
        while representative_dists:
            head = heapq.heappop(representative_dists)
            if head[0] <= self.threshold or (head[0] > self.threshold and obs_num < self.baseobs):
                _, group = self.representatives[head[1]]
                for thing in group:
                    obs_dist[thing] = self.dist_function(query_seq, self.refs[thing])
                    obs_num += 1
            else:
                break

        return obs_dist