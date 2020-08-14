from apples.criteria import OLS, FM, BE
import sys

class Query_context():
    def __init__(self, treecore, treecore_frag, options, reference):
        self.treecore = treecore
        self.treecore_frag = treecore_frag
        self.options = options
        self.reference = reference

    def runquery(self, query_name, query_seq, obs_dist):
        jplace = dict()
        jplace["placements"] = [{"p": [[0, 0, 1, 0, 0]], "n": [query_name]}]

        if not obs_dist:
            obs_dist = self.reference.get_obs_dist(query_seq, query_name)
        else:
            tx=0
            for k, v in sorted(obs_dist.items(), key=lambda kv: kv[1]):
                if v == -1:
                    continue
                tx += 1
                if tx > self.options.base_observation_threshold and v > self.options.filt_threshold:
                    obs_dist[k] = -1

        for l in self.treecore.tree.traverse_postorder(internal=False):
            if l.label not in obs_dist:
                raise ValueError('Taxon {} should be in distances table.'.format(l.label))
            if obs_dist[l.label] == 0:
                jplace["placements"][0]["p"][0][0] = l.edge_index
                return jplace

        # number of non-infinity distances
        tx = len([ v for k, v in obs_dist.items() if v > -1 ])

        if tx <= 2:
            sys.stderr.write('Taxon {} cannot be placed. At least three non-infinity distances '
                             'should be observed to place a taxon. '
                             'Consequently, this taxon is ignored (no output).\n'.format(query_name))
            jplace["placements"][0]["p"][0][0] = -1
            return jplace

        if -1 not in obs_dist.values():
            tc = self.treecore
            tc.dp(obs_dist)
        else:
            tc = self.treecore_frag
            tc.validate_edges(obs_dist)
            tc.dp_frag(obs_dist)

        if self.options.method_name == "BE":
            alg = BE(tc.tree)
        elif self.options.method_name == "FM":
            alg = FM(tc.tree)
        else:
            alg = OLS(tc.tree)
        alg.placement_per_edge(self.options.negative_branch)
        jplace["placements"][0]["p"] = [alg.placement(self.options.criterion_name, query_name)]
        return jplace
