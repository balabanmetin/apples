import itertools
import subprocess
import tempfile
from abc import ABC, abstractmethod

from apples.PoolRepresentativeWorker import PoolRepresentativeWorker
from apples.distance import *
from apples.fasta2dic import fasta2dic
import heapq
import time
import logging
import multiprocessing as mp


class Reference(ABC):
    """prot flag true if protein sequences, false for nucleotide sequences"""
    def __init__(self, ref_fp: str, prot_flag: bool):
        self.refs = fasta2dic(ref_fp, prot_flag, False)  # never mask low confidence bases in the reference
        self.prot_flag = prot_flag

        if prot_flag:
            self.dist_function = scoredist
        else:
            self.dist_function = jc69

    @abstractmethod
    def get_obs_dist(self, query_seq, query_name):
        pass


class FullReference(Reference):
    def __init__(self, ref_fp: str, prot_flag: bool):
        Reference.__init__(self, ref_fp, prot_flag)

    def get_obs_dist(self, query_seq, query_tag):
        obs_dist = {query_tag: 0}
        for tagr, seqr in self.refs:
            obs_dist[tagr] = self.dist_function(query_seq, seqr)
        return obs_dist


class ReducedReference(Reference):
    def __init__(self, ref_fp, prot_flag, tree_file, threshold, num_thread):
        Reference.__init__(self, ref_fp, prot_flag)
        self.threshold = threshold

        start = time.time()
        tc_output_file = tempfile.NamedTemporaryFile(delete=True, mode='w+t').name
        nldef = tempfile.NamedTemporaryFile(delete=True, mode='w+t')
        s = ["TreeCluster.py", "-t", str(self.threshold * 1.2), "-i", tree_file, "-m", "max", "-o", tc_output_file]
        subprocess.call(s, stdout=nldef, stderr=nldef)
        logging.info(
            "[%s] TreeCluster is completed in %.3f seconds." % (time.strftime("%H:%M:%S"),(time.time() - start)))

        start = time.time()
        with open(tc_output_file, "r") as tc_output:
            tc_output.readline()
            lines = map(lambda x: x.strip().split('\t'), tc_output.readlines())
            lines_sorted = sorted(lines, key=lambda x: x[1])
            clusters = [(key, [i[0] for i in list(group)]) for key, group in itertools.groupby(lines_sorted, lambda x: x[1])]
            partition_worker = PoolRepresentativeWorker()
            partition_worker.set_class_attributes(self.refs, self.prot_flag)
            pool = mp.Pool(num_thread)
            results = pool.map(partition_worker.worker, clusters)
            pool.close()
            pool.join()
            self.representatives = [item for sublist in results for item in sublist]

        logging.info(
            "[%s] Representative sequences are computed in %.3f seconds." %
            (time.strftime("%H:%M:%S"), (time.time() - start)))


    def set_baseobs(self, baseobs):
        self.baseobs = baseobs


    def get_obs_dist(self, query_seq, query_tag):
        obs_dist = {}
        obs_num = 0
        representative_dists = []
        i = 0
        for consensus_seq, group in self.representatives:
            dist = self.dist_function(query_seq, consensus_seq)
            if dist < 0: # invalid distance
                continue
            representative_dists.append((dist, i))
            i += 1
        heapq.heapify(representative_dists)
        while representative_dists:
            head = heapq.heappop(representative_dists)
            if head[0] <= self.threshold or (head[0] > self.threshold and obs_num < self.baseobs):
                _, group = self.representatives[head[1]]
                for thing in group:
                    thing_d = self.dist_function(query_seq, self.refs[thing])
                    if not thing_d < 0:
                        obs_dist[thing] = thing_d
                        obs_num += 1
            else:
                break

        return obs_dist
