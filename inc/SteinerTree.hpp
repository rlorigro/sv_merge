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
    size_t n_samples, n_solutions, n_haps;
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
     * Adds to `out[i]` the number of samples that can be assigned to `i` solutions.
     */
    void get_sample_solutions_histogram(vector<int64_t> out);

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
     * Writes the current approximation to `output_dir`, and deletes every read-hap edge of `transmap` that does not
     * belong to the approximation. This is similar to `path_optimizer_mathopt.parse_read_model_solution()`.
     */
    void parse_approximation(TransMap& transmap, path output_dir);




    // ------------------------------------------- GREEDY DENSE SOLUTION -----------------------------------------------
    /**
     * @param covered_samples temporary space;
     */
    double get_density(int64_t solution_id, unordered_set<int64_t>& covered_samples, float hap_cost, float edge_cost_multiplier);

    /**
     * A simpler version of `get_density()` to be used when no solution has been selected yet.
     */
    double get_density_init(int64_t solution_id, float hap_cost, float edge_cost_multiplier);

    /**
     * Priority of the max-queue used by `greedy_dense_solution()`.
     */
    class Compare {
    public:
        bool operator() (pair<int64_t,double> x, pair<int64_t,double> y) { return x.second<y.second; }
    };

    /**
     * Iteratively adds a solution with largest density (or "cost-effectiveness"), where the density is the ratio
     * between the number of new samples that can be covered, and the cost of assigning such samples to the solution.
     * This is inspired by the usual greedy approximation of min-weight set cover and of Steiner tree: see e.g.
     *
     * Charikar et al. "Approximation algorithms for directed Steiner problems." Journal of Algorithms 33.1 (1999).
     *
     * Remark: to minimize/maximize only `d` exactly, use `approximate_shortest_paths()`. To minimize/maximize just `n`,
     * set `edge_cost_multiplier=0`: this gives just an approximation.
     *
     * @param hap_cost cost of a haplotype;
     * @param edge_cost_multiplier multiplier of the cost of each solution-sample edge;
     * @param mode 0=computes a single sequence of greedy steps; 1=starts a sequence of greedy steps from every solution
     * and selects the min;
     * @param n_samples_with_solutions adds to position `i` the number of samples with `i` possible solutions;
     * @return the objective value of the approximation.
     */
    double greedy_dense_solution(float hap_cost, float edge_cost_multiplier, int64_t mode, vector<int64_t>& n_samples_with_solutions);

    /**
     * A sequence of greedy steps driven by solution density.
     *
     * @param density temporary space;
     * @param queue temporary space;
     * @param covered_samples temporary space;
     * @param covered_samples_new temporary space;
     * @param solutions_to_update temporary space.
     */
    void greedy_dense_solution_impl(vector<double>& density, priority_queue<pair<int64_t,double>, vector<pair<int64_t,double>>, Compare>& queue, unordered_set<int64_t>& covered_samples, unordered_set<int64_t>& covered_samples_new, unordered_set<int64_t>& solutions_to_update, float hap_cost, float edge_cost_multiplier);

    /**
     * The main step of `greedy_dense_solution_impl()`: `solution_id` is added to the current solution and some
     * densities are updated.
     *
     * @param density temporary space;
     * @param covered_samples temporary space;
     * @param covered_samples_new temporary space;
     * @param solutions_to_update  temporary space;
     * @param update_queue FALSE=updates only `density`; TRUE=updates also `queue`;
     * @param queue the priority queue used in every iteration.
     */
    void greedy_dense_solution_step(int64_t solution_id, vector<double>& density, unordered_set<int64_t>& covered_samples, unordered_set<int64_t>& covered_samples_new, unordered_set<int64_t>& solutions_to_update, bool update_queue, priority_queue<pair<int64_t,double>, vector<pair<int64_t,double>>, Compare> queue, float hap_cost, float edge_cost_multiplier);




    // ----------------------------------------------- SHORTEST PATHS --------------------------------------------------
    /**
     * Uses a single-source shortest path tree of the tripartite graph (sample,solution,root) as an approximation. This
     * does not model the fact that different solutions can share haplotypes.
     *
     * Remark: to minimize/maximize only `d`, set `hap_cost=0`: this solves the problem exactly. To minimize/maximize
     * just `n`, set `edge_cost_multiplier=0`: this gives just an approximation.
     *
     * @param hap_cost cost of a haplotype;
     * @param edge_cost_multiplier multiplier of the cost of each solution-sample edge;
     * @param n_samples_with_solutions adds to position `i` the number of samples with `i` possible solutions;
     * @return the objective value of the approximation.
     */
    double approximate_shortest_paths(bool minimize, float hap_cost, float edge_cost_multiplier, bool build_solution, vector<int64_t>& n_samples_with_solutions);

};


}
