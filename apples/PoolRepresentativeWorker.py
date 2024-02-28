import numpy as np


class PoolRepresentativeWorker:
    refs = None
    prot_flag = None

    @classmethod
    def set_class_attributes(cls, refs, prot_flag):
        """
        Set class attributes with the given references and protection flag.
        """
        cls.refs = refs
        cls.prot_flag = prot_flag

    @staticmethod
    def _find_representative(group, refs, prot_flag):
        """
        Finds the representative sequence for a given group based on a set of reference sequences
        and a flag indicating whether the sequences are protein or not.
        It calculates the consensus sequence for the input group of sequences using a specific alphabet
        and returns the consensus sequence.

        Args:
            group: List of indices representing the group of sequences
            refs: Dictionary mapping sequence indices to actual sequences
            prot_flag: Boolean indicating whether the sequences are protein sequences or not

        Returns:
            Numpy array representing the representative sequence for the given group
        """

        if prot_flag:
            alphabet = np.array(
                [
                    'A',
                    'C',
                    'D',
                    'E',
                    'F',
                    'G',
                    'H',
                    'I',
                    'K',
                    'L',
                    'M',
                    'N',
                    'P',
                    'Q',
                    'R',
                    'S',
                    'T',
                    'V',
                    'W',
                    'Y',
                    '-',
                ],
                dtype='S1',
            )

        else:
            alphabet = np.array([b'A', b'C', b'G', b'T', b'-'], dtype='S1')

        lookup = {n: i for i, n in enumerate(alphabet)}

        def get_consensus(arr):
            """
            Compute the consensus sequence from a 2D array of DNA or protein sequences.

            Parameters:
                arr (numpy.ndarray): The 2D array of DNA or protein sequences.

            Returns:
                numpy.ndarray: The consensus sequence.
            """
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
        """
        Class method for processing a worker cluster.

        Args:
            cls: The class itself.
            cluster: A tuple representing the cluster to be processed.

        Returns:
            A list of tuples containing processed data based on the cluster key.
        """
        key, group = cluster
        if key == '-1':
            return [(cls.refs[thing], [thing]) for thing in group]
        else:
            return [(cls._find_representative(group, cls.refs, cls.prot_flag), group)]
