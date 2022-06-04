import time
import logging
from apples.OLS import OLS
from apples.FM import FM
from apples.BE import BE
from apples.BME import BME
import sys
from apples.Subtree import Subtree
import math


class PoolQueryWorker:
    reference = None
    options = None
    name_to_node_map = None
    treecore = None


    @classmethod
    def set_class_attributes(cls, reference, options, name_to_node_map):
        cls.reference = reference
        cls.options = options
        cls.name_to_node_map = name_to_node_map
        # cls.tree = tree

    @classmethod
    def runquery(cls, query_name, query_seq, obs_dist):
        jplace = dict()
        jplace["placements"] = [{"p": [[0, 0, 1, 0, 0]], "n": [query_name]}]

        start_dist = time.time()
        if not obs_dist:
            obs_dist = cls.reference.get_obs_dist(query_seq, query_name, cls.options.minimum_alignment_overlap)
        else:
            def valid_dists(obs_dist, name_to_node_map):
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

        if len(obs_dist) <= 3:
            start_dp = time.time()
            sys.stderr.write('Taxon {} cannot be placed. At least three non-infinity distances '
                             'should be observed to place a taxon. '
                             'Consequently, this taxon is ignored (no output).\n'.format(query_name))
            jplace["placements"][0]["p"][0][0] = -1
            end_dp = time.time() - start_dp
            logging.info(
                "[%s] Distances are computed for query %s in %.3f seconds.\n"
                "[%s] Dynamic programming is completed for query %s in %.3f seconds." %
                (time.strftime("%H:%M:%S"), query_name, end_dist,
                 time.strftime("%H:%M:%S"), query_name, end_dp))
            return jplace

        start_dp = time.time()
        subtree = Subtree(obs_dist, cls.name_to_node_map, query_name)
        distr = subtree.edge_contribution_counter()
        largest_distr = max(list(distr.values()))
        if largest_distr >= math.floor(math.log2(len(cls.name_to_node_map))):
            for k,v in reversed(obs_dist.items()):
                if distr[k] == largest_distr:
                    #print(query_name, k, largest_distr)
                    discard_index = k
                    break
        else:
            largest_dist = max(obs_dist.values())
            for k,v in reversed(obs_dist.items()):
                if v == largest_dist:
                    #print(query_name, k, largest_dist)
                    discard_index = k
                    break
        subtree.unroll_changes()
        del obs_dist[discard_index]
        subtree = Subtree(obs_dist, cls.name_to_node_map, query_name)

        #print(subtree.query_name, subtree.edge_contribution_counter())
        for k, v in obs_dist.items():
            if v == 0 and k != query_name:
                jplace["placements"][0]["p"][0][0] = cls.name_to_node_map[k].edge_index
                subtree.unroll_changes()
                return jplace


        if cls.options.method_name == "BE":
            alg = BE(subtree)
        elif cls.options.method_name == "FM":
            alg = FM(subtree)
        elif cls.options.method_name == "BME":
            alg = BME(subtree)
        else:
            alg = OLS(subtree)
        alg.dp_frag()

        alg.placement_per_edge(cls.options.negative_branch)
        presult, potential_misplacement_flag = alg.placement(cls.options.criterion_name)
        jplace["placements"][0]["p"] = [presult]
        if potential_misplacement_flag == 1:
            if cls.options.exclude_intplace:
                jplace["placements"][0]["p"][0][0] = -1
                ignoredprompt = " Consequently, this sequence is ignored (no output)."
            else:
                ignoredprompt = ""
            logging.warning(
                "Best placement for query sequence %s has zero pendant edge length and placed at an internal node "
                "with a non-zero least squares error. This is a potential misplacement.%s" % (query_name, ignoredprompt))

        subtree.unroll_changes()

        end_dp = time.time() - start_dp
        logging.info(
            "[%s] Distances are computed for query %s in %.3f seconds.\n"
            "[%s] Dynamic programming is completed for query %s in %.3f seconds." %
            (time.strftime("%H:%M:%S"), query_name, end_dist,
             time.strftime("%H:%M:%S"), query_name, end_dp))
        return jplace
