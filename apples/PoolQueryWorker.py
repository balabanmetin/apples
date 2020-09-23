from apples.criteria import OLS, FM, BE
import sys
from apples.Subtree import Subtree


class PoolQueryWorker:
    reference = None
    options = None
    name_to_node_map = None
    treecore = None


    @classmethod
    def set_class_attributes(cls, reference, options, name_to_node_map, treecore):
        cls.reference = reference
        cls.options = options
        cls.name_to_node_map = name_to_node_map
        cls.treecore = treecore

    @classmethod
    def runquery(cls, query_name, query_seq, obs_dist):
        jplace = dict()
        jplace["placements"] = [{"p": [[0, 0, 1, 0, 0]], "n": [query_name]}]

        if not obs_dist:
            obs_dist = cls.reference.get_obs_dist(query_seq, query_name)
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

        if len(obs_dist) <= 2:
            sys.stderr.write('Taxon {} cannot be placed. At least three non-infinity distances '
                             'should be observed to place a taxon. '
                             'Consequently, this taxon is ignored (no output).\n'.format(query_name))
            jplace["placements"][0]["p"][0][0] = -1
            return jplace

        for k, v in obs_dist.items():
            if v == 0 and k != query_name:
                jplace["placements"][0]["p"][0][0] = cls.name_to_node_map[k].edge_index
                return jplace

        tc = cls.treecore
        subtree = Subtree(obs_dist, cls.name_to_node_map)
        tc.dp_frag(subtree)

        if cls.options.method_name == "BE":
            alg = BE(tc.tree)
        elif cls.options.method_name == "FM":
            alg = FM(tc.tree)
        else:
            alg = OLS(tc.tree)
        alg.placement_per_edge(cls.options.negative_branch)
        jplace["placements"][0]["p"] = [alg.placement(cls.options.criterion_name)]
        subtree.unroll_changes()
        return jplace
