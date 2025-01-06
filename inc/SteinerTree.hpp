#pragma once

#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <string>
#include <vector>
#include <tuple>

using std::unordered_map;
using std::unordered_set;
using std::to_string;
using std::vector;
using std::string;
using std::pair;

#include "TransitiveMap.hpp"


namespace sv_merge {


class SteinerTree {
    /**
     * The tripartite graph `(S,H,{r})` on which Steiner tree approximations are computed. `S` contains one node per
     * sample and `H` contains one node per possible solution (a solution is either one or two haplotypes that cover all
     * the reads of at least one sample). An `(s,h)` edge means that all the reads from sample `s` can be covered by a
     * haplotype in `h`; the weight is the smallest cost of such a cover. The weight of every `(h,r)` edge is the cost
     * of `h`. The graph is directed outwards from `r`, which is not represented explicitly.
     */
    vector<pair<int64_t,int64_t>> solutions;  // The list of all distinct solutions
    vector<vector<int64_t>> solution_samples;  // The sample neighbors of each solution
    vector<vector<float>> solution_sample_weights;  // The weight of each solution-sample edge

    float hap_cost, edge_cost_multiplier;

public:
    /**
     * Builds the directed graph on which Steiner tree approximations are computed.
     *
     * @param transmap `transmap.sort_adjacency_lists()` must have already been called;
     * @param hap_cost cost of a haplotype;
     * @param edge_cost_multiplier multiplier of the cost of each solution-sample edge.
     */
    SteinerTree(TransMap& transmap, float hap_cost, float edge_cost_multiplier);

    /**
     *
     */
    void dense_bunch();




};


}
