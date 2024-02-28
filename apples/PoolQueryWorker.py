import time
import logging
from apples.OLS import OLS
from apples.FM import FM
from apples.BE import BE
from apples.BME import BME
import sys
from apples.Subtree import Subtree


class PoolQueryWorker:
    reference = None
    options = None
    name_to_node_map = None
    treecore = None

    @classmethod
    def set_class_attributes(cls, reference, options, name_to_node_map):
        """
        Set class attributes for the given reference, options, and name-to-node map.
        """
        cls.reference = reference
        cls.options = options
        cls.name_to_node_map = name_to_node_map
        # cls.tree = tree

    @classmethod
    def runquery(cls, query_name, query_seq, obs_dist):
        """
        A class method that runs a query and returns a data structure with placement information.
        It takes in the query name, query sequence, and observed distances as parameters.
        The return type is a dictionary containing placement information.
        It computes placements for a query based on observed distances, and it handles various cases
        such as insufficient distances and potential misplacements. It also logs computation and completion times.
        """
        jplace = dict()
        jplace['placements'] = [{'p': [[0, 0, 1, 0, 0]], 'n': [query_name]}]

        start_dist = time.time()
        if not obs_dist:
            obs_dist = cls.reference.get_obs_dist(query_seq, query_name, cls.options.minimum_alignment_overlap)
        else:

            def valid_dists(obs_dist, name_to_node_map):
                """
                This function takes in two parameters: obs_dist, a dictionary, and name_to_node_map, a map.
                It iterates through the items in observed distances, sorted by their distance,
                and yields the key-value pairs that do not exceed the filter threshold.
                """
                tx = 0
                for k, v in sorted(obs_dist.items(), key=lambda kv: kv[1]):
                    if v < 0 or k not in name_to_node_map:
                        continue
                    tx += 1
                    if tx > cls.options.base_observation_threshold and v > cls.options.filt_threshold:
                        break
                    yield (k, v)

            obs_dist = dict(valid_dists(obs_dist, cls.name_to_node_map))
        end_dist = time.time() - start_dist

        start_dp = time.time()
        if query_name in cls.name_to_node_map:
            if query_name in obs_dist:
                del obs_dist[query_name]
            logging.warning(
                'The query named %s exists in the backbone. Changing its name to %s-query.' % (query_name, query_name)
            )
            query_name = query_name + '-query'
            jplace['placements'][0]['n'] = [query_name]

        for k, v in obs_dist.items():
            if v == 0:
                jplace['placements'][0]['p'][0][0] = cls.name_to_node_map[k].edge_index
                return jplace

        def not_sufficient_distances():
            """
            A function that handles cases where there are not sufficient distances for placing a taxon.
            It writes an error message to stderr and logs the computation and completion times.
            It then returns the updated jplace data structure.
            """
            sys.stderr.write(
                'Taxon {} cannot be placed. At least three non-infinity distances '
                'should be observed to place a taxon. '
                'Consequently, this taxon is ignored (no output).\n'.format(query_name)
            )
            jplace['placements'][0]['p'][0][0] = -1
            end_dp = time.time() - start_dp
            logging.info(
                '[%s] Distances are computed for query %s in %.3f seconds.\n'
                '[%s] Dynamic programming is completed for query %s in %.3f seconds.'
                % (time.strftime('%H:%M:%S'), query_name, end_dist, time.strftime('%H:%M:%S'), query_name, end_dp)
            )
            return jplace

        if len(obs_dist) <= 2:
            return not_sufficient_distances()

        # extract the subtree
        subtree = Subtree(obs_dist, cls.name_to_node_map)

        # select the least squares weighting scheme
        if cls.options.method_name == 'BE':
            alg = BE(subtree)
        elif cls.options.method_name == 'FM':
            alg = FM(subtree)
        elif cls.options.method_name == 'BME':
            alg = BME(subtree)
        else:
            alg = OLS(subtree)

        # call dynamic programming
        alg.dp_frag()

        # call best placement for all subtree edges
        alg.placement_per_edge(cls.options.negative_branch)
        # find the edge with the least squares error (weighed)
        presult, potential_misplacement_flag = alg.placement(cls.options.criterion_name)
        jplace['placements'][0]['p'] = [presult]
        if potential_misplacement_flag == 1:
            if cls.options.exclude_intplace:
                jplace['placements'][0]['p'][0][0] = -1
                ignoredprompt = ' Consequently, this sequence is ignored (no output).'
            else:
                ignoredprompt = ''
            logging.warning(
                'Best placement for query sequence %s has zero pendant edge length and placed at an internal node '
                'with a non-zero least squares error. This is a potential misplacement.%s' % (query_name, ignoredprompt)
            )

        # reset the subtree edges valid flags
        subtree.unroll_changes()

        end_dp = time.time() - start_dp
        logging.info(
            '[%s] Distances are computed for query %s in %.3f seconds.\n'
            '[%s] Dynamic programming is completed for query %s in %.3f seconds.'
            % (time.strftime('%H:%M:%S'), query_name, end_dist, time.strftime('%H:%M:%S'), query_name, end_dp)
        )
        return jplace
