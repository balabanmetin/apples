import numpy as np


class PoolRepresentativeWorker:
    refs = None
    prot_flag = None

    @classmethod
    def set_class_attributes(cls, refs, prot_flag):
        cls.refs = refs
        cls.prot_flag = prot_flag

    @staticmethod
    def _find_representative(group, refs, prot_flag):

        if prot_flag:
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

        group_seqs = [refs[g] for g in group]
        group_mat = np.array(np.vstack(group_seqs))
        return get_consensus(group_mat)

    @classmethod
    def worker(cls, cluster):
        key, group = cluster
        if key == "-1":
            return [(cls.refs[thing], [thing]) for thing in group]
        else:
            return [(cls._find_representative(group, cls.refs, cls.prot_flag), group)]