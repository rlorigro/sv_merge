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


namespace std {

template <class T>
inline void hash_combine(size_t& seed, T const& v) {
    seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
}

template <typename A, typename B>
struct hash<tuple<A,B>> {
    std::size_t operator()(const tuple<A,B>& k) const {
        size_t hash_val = std::hash<A>()(std::get<0>(k));
        hash_combine(hash_val,std::get<1>(k));
        return hash_val;
    }
};

template <typename A, typename B, typename C>
struct hash<tuple<A,B,C>> {
    std::size_t operator()(const tuple<A,B,C>& k) const {
        size_t hash_val = std::hash<A>()(std::get<0>(k));
        hash_combine(hash_val,std::get<1>(k));
        hash_combine(hash_val,std::get<2>(k));
        return hash_val;
    }
};

}


namespace sv_merge {

class SteinerTree {
    /**
     * The tripartite graph `(S,H,{r})` on which Steiner tree approximations are computed. `S` contains one node per
     * sample and `H` contains one node per possible solution (a solution is either one or two haplotypes that cover all
     * the reads of at least one sample). An `(s,h)` edge means that all the reads from sample `s` can be covered by a
     * haplotype in `h`; the weight is the smallest cost of such a cover. The weight of every `(h,r)` edge is the cost
     * of `h`. The graph is directed outwards from `r`, which is not represented explicitly.
     */
    size_t n_samples, n_solutions;
    vector<pair<int64_t,int64_t>> solutions;  // List of all distinct solutions, not necessarily sorted.
    vector<vector<int64_t>> solution_samples;  // solutionID -> sampleID
    vector<vector<double>> solution_sample_weights;  // solutionID -> solutionSampleWeight (sum of from read-hap weights in the transmap)
    unordered_map<int64_t,vector<pair<int64_t,double>>> sample_solutions;  // sampleID -> (solutionID,sampleSolutionWeight)
    unordered_map<int64_t,vector<int64_t>> hap_solutions;  // hapID -> solutionID

    /**
     * Data structures needed by `parse_approximation()`.
     */
    unordered_map<tuple<int64_t,int64_t,int64_t>, set<pair<int64_t,int64_t>>> solution_sample_reads;

    /**
     * A Steiner tree approximation
     */
    double objective;
    unordered_set<int64_t> selected_haps;  // Haplotypes in selected solutions.
    vector<bool> selected_solutions;  // Marks `solutions`.
    vector<vector<bool>> selected_edges;  // Marks `solution_samples` cells.

    /**
     * Clears the state of the approximation
     */
    void approximate_clear();

    /**
     * Builds object variable `sample_solutions`.
     *
     * Remark: the procedure does nothing if the variables are already non-empty.
     */
    void build_sample_solutions();

    /**
     * Builds object variable `hap_solutions`, assuming that `solutions` has already been built.
     *
     * Remark: the procedure does nothing if the variables are already non-empty.
     */
    void build_hap_solutions();


public:
    /**
     * Builds the directed graph on which Steiner tree approximations are computed.
     *
     * @param transmap `transmap.sort_adjacency_lists()` must have already been called.
     */
    SteinerTree(TransMap& transmap);

    /**
     * @param minimize computes the smallest (TRUE) or largest (FALSE) value of `d`;
     * @param build_solution TRUE: builds an optimal solution, which can be retrieved with `parse_approximation()`;
     * @return the optimal value of `d`.
     */
    double optimize_d(bool minimize, bool build_solution);

    /**
     * Uses the single-source shortest path tree as an approximation of the directed Steiner tree. This gives an
     * approximation factor of `n_samples`.
     *
     * @param hap_cost cost of a haplotype;
     * @param edge_cost_multiplier multiplier of the cost of each solution-sample edge;
     * @return the objective value of the approximation.
     */
    double approximate_shortest_paths(float hap_cost, float edge_cost_multiplier);

    /**
     * Priority queue comparison criterion
     */
    class Compare {
    public:
        bool operator() (pair<int64_t,double> x, pair<int64_t,double> y) {
            return x.second<y.second;
        }
    };

    /**
     * Greedily adds a solution with largest density, where the density is the ratio between the number of new covered
     * samples and the cost of assigning them to the solution. This implements the algorithm in Figure 1 of:
     *
     * Charikar et al. "Approximation algorithms for directed Steiner problems." Journal of Algorithms 33.1 (1999).
     *
     * This gives an approximation factor of `\sqrt{n_samples}`.
     *
     * @param hap_cost cost of a haplotype;
     * @param edge_cost_multiplier multiplier of the cost of each solution-sample edge;
     * @return the objective value of the approximation.
     */
    double greedy_dense_solution(float hap_cost, float edge_cost_multiplier);

    /**
     * @param solutions_to_update temporary space.
     */
    void greedy_dense_solution_impl(vector<double>& density, priority_queue<pair<int64_t,double>, vector<pair<int64_t,double>>, Compare>& queue, unordered_set<int64_t>& covered_samples, unordered_set<int64_t>& solutions_to_update, float hap_cost, float edge_cost_multiplier);

    /**
     * Initializes `greedy_dense_solution()` from a specific solution.
     *
     * @param solutions_to_update temporary space.
     */
    void greedy_dense_solution_set_initial_conditions(int64_t solution_id, vector<double>& density, unordered_set<int64_t>& covered_samples, unordered_set<int64_t>& solutions_to_update, float hap_cost, float edge_cost_multiplier);

    /**
     * One iteration of procedure `greedy_dense_solution()`.
     *
     * @param covered_samples temporary space;
     * @param covered_samples_new temporary space;
     * @param solutions_to_update  temporary space;
     * @param update_queue TRUE=updates `queue`;
     * @param queue priority queue used in every iteration.
     */
    void greedy_dense_solution_iteration(int64_t solution_id, vector<double>& density, unordered_set<int64_t>& covered_samples, unordered_set<int64_t>& covered_samples_new, unordered_set<int64_t>& solutions_to_update, bool update_queue, priority_queue<pair<int64_t,double>, vector<pair<int64_t,double>>, Compare> queue, float hap_cost, float edge_cost_multiplier);

    /**
     * Writes the current approximation to `output_dir`, and deletes every read-hap edge of `transmap` that does not
     * belong to the approximation. This is similar to `path_optimizer_mathopt.parse_read_model_solution()`.
     */
    void parse_approximation(TransMap& transmap, path output_dir);
};


}
