import itertools
import subprocess
import tempfile
from abc import ABC, abstractmethod

from apples.PoolRepresentativeWorker import PoolRepresentativeWorker
from apples.distance import scoredist, jc69
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
    def get_obs_dist(self, query_seq, query_name, overlap_frac):
        pass


class FullReference(Reference):
    def __init__(self, ref_fp: str, prot_flag: bool):
        """
        Initializes the object with the provided reference file path and protein flag.

        Parameters:
            ref_fp (str): The file path to the reference.
            prot_flag (bool): The protein flag indicating whether the sequences are protein of nucleotide.

        Returns:
            None
        """
        Reference.__init__(self, ref_fp, prot_flag)

    def get_obs_dist(self, query_seq, query_tag, overlap_frac):
        """
        Calculate the observed distances for a query sequence and taxon tag, based on the provided overlap fraction.

        Parameters:
            query_seq (str): The query sequence for which the observed distances is calculated.
            query_tag (str): The tag associated with the query sequence.
            overlap_frac (float): The fraction of overlap to be considered in the calculation.

        Returns:
            dict: Dictionary representing the observed distances, with the query tag as the key and the distance value.
        """
        obs_dist = {query_tag: 0}
        for tagr, seqr in self.refs:
            obs_dist[tagr] = self.dist_function(query_seq, seqr, overlap_frac)
        return obs_dist


class ReducedReference(Reference):
    def __init__(self, ref_fp, prot_flag, tree_file, threshold, num_thread):
        """
        Constructor __init__ that takes in several parameters and initializes the object.
        It then proceeds to run an external command using the subprocess module to call a program called TreeCluster,
        which clusters the input tree into subsets of trees.
        After that, it reads and processes the output of TreeCluster,
        and then uses multiprocessing to compute representative sequences for each subset.
        Finally, it logs the time taken for these operations.

        Parameters:
            ref_fp: str - file path for the reference
            prot_flag: bool - flag for protocol
            tree_file: str - file path for the tree
            threshold: float - threshold value for Treecluster
            num_thread: int - number of threads
        """
        Reference.__init__(self, ref_fp, prot_flag)
        self.threshold = threshold

        start = time.time()
        tc_output_file = tempfile.NamedTemporaryFile(delete=True, mode='w+t').name
        nldef = tempfile.NamedTemporaryFile(delete=True, mode='w+t')
        s = ['TreeCluster.py', '-t', str(self.threshold * 1.2), '-i', tree_file, '-m', 'max', '-o', tc_output_file]
        subprocess.call(s, stdout=nldef, stderr=nldef)
        logging.info(
            '[%s] TreeCluster is completed in %.3f seconds.' % (time.strftime('%H:%M:%S'), (time.time() - start))
        )

        start = time.time()
        with open(tc_output_file, 'r') as tc_output:
            tc_output.readline()
            lines = map(lambda x: x.strip().split('\t'), tc_output.readlines())
            lines_sorted = sorted(lines, key=lambda x: x[1])
            clusters = [
                (key, [i[0] for i in list(group)]) for key, group in itertools.groupby(lines_sorted, lambda x: x[1])
            ]
            partition_worker = PoolRepresentativeWorker()
            partition_worker.set_class_attributes(self.refs, self.prot_flag)
            pool = mp.Pool(num_thread)
            results = pool.map(partition_worker.worker, clusters)
            pool.close()
            pool.join()
            self.representatives = [item for sublist in results for item in sublist]

        logging.info(
            '[%s] Representative sequences are computed in %.3f seconds.'
            % (time.strftime('%H:%M:%S'), (time.time() - start))
        )

    def set_baseobs(self, baseobs):
        self.baseobs = baseobs

    def get_obs_dist(self, query_seq, query_tag, overlap_frac):
        """
        Calculate the distance between the query sequence and the representatives, and return the observed distances.
        It iterates through the representatives, calculates distances using a specified function,
        and populates the obs_dist dictionary with the observed distances.
        The method also uses a heap data structure to efficiently process representative distances.
        The code computes the distance between the query and the references and prioritizes the clusters
        whose representatives have a lower distance to the query. It uses heaps to efficiently prioritize and process
        the representative distances, ensuring that the clusters with lower distances are processed first.

        Parameters:
        - query_seq: The query sequence to compare with the representatives.
        - query_tag: The tag associated with the query sequence.
        - overlap_frac: The overlap fraction for distance calculation.

        Returns:
        - obs_dist: A dictionary containing the observed distances between the query sequence and the representatives.
        """
        obs_dist = {}
        obs_num = 0
        representative_dists = []
        for i, content in enumerate(self.representatives):
            consensus_seq, group = content
            dist = self.dist_function(query_seq, consensus_seq, overlap_frac)
            if dist >= 0:  # valid distance
                representative_dists.append((dist, i))
        heapq.heapify(representative_dists)
        while representative_dists:
            head = heapq.heappop(representative_dists)
            if head[0] <= self.threshold or obs_num < self.baseobs:
                _, group = self.representatives[head[1]]
                for thing in group:
                    thing_d = self.dist_function(query_seq, self.refs[thing], overlap_frac)
                    if not thing_d < 0:
                        obs_dist[thing] = thing_d
                        obs_num += 1
            else:
                break
        # for k,v in obs_dist.items():
        #    print(k+"\t"+str(v))
        return obs_dist
