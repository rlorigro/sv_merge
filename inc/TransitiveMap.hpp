#pragma once

#include "VectorHeteroGraph.hpp"
#include "BinarySequence.hpp"
#include "VariantGraph.hpp"
#include "misc.hpp"

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
using std::byte;
using std::tuple;


namespace sv_merge {


class TransMap {
    /// Attributes
    HeteroGraph<HeteroNode> graph;

    // Pains me to add yet another map but here it is
    unordered_map <int64_t, BinarySequence<uint64_t> > sequences;

    // Pains me to add yet another map but here it is
    unordered_map<int64_t,coord_t> sequence_flanks;

    // Pains me to add yet another map but here it is
    unordered_map<int64_t,bool> sequence_reversals;

    /**
     * Data structures used by compression/decompression procedures.
     */
    vector<int64_t> read_ids, cluster_ids;
    unordered_map<string,tuple<int64_t,int64_t,int64_t,int64_t>> sample_to_sample;  // from_sample_name -> from_sample_first, from_sample_last, to_sample_first, to_sample_last
    unordered_set<int64_t> solved_samples;
    vector<string> partitioned_samples;
    bool is_compressed = false;

public:
    static const string sample_node_name;
    static const string read_node_name;
    static const string path_node_name;
    static const string variant_node_name;

    /**
     * Haps that must be set to one in the ILP
     */
    unordered_set<int64_t> present_haps;

    /**
     * Edges that must be set to one in the ILP. Format: (read_id,hap_id).
     */
    unordered_set<pair<int64_t,int64_t>> present_edges;

    TransMap();

    /// Building
    void reserve_nodes(size_t n);
    void reserve_edges(size_t n);
    void reserve_sequences(size_t n);
    void add_sample(const string& name);
    void add_read(const string& name);
    void add_flank_coord(const string& name, int32_t start, int32_t stop);
    void add_read(const string& name, const string& sequence);
    void add_read_with_move(string& name, BinarySequence<uint64_t>& sequence);
    void add_read_with_move(string& name, BinarySequence<uint64_t>& sequence, bool reversal);
    void add_path(const string& name);
    void add_variant(const string& name);
    void add_edge(const string& a, const string& b);
    void add_edge(int64_t a, int64_t b, float weight);
    void add_edge(const string& a, const string& b, float weight);

    void remove_edge(int64_t a, int64_t b);
    void remove_node(int64_t id);

    /// Accessing
    bool empty() const;
    bool get_is_compressed() const;

    size_t get_edge_count() const;
    size_t get_read_count() const;
    size_t get_path_count() const;
    size_t get_sample_count() const;

    int64_t get_id(const string& name) const;
    pair<bool,int64_t> try_get_id(const string& name) const;
    pair<bool,float> try_get_edge_weight(int64_t id_a, int64_t id_b) const;
    bool has_edge(int64_t a, int64_t b) const;
    bool has_node(const string& name) const;
    bool has_node(int64_t id) const;
    int64_t get_node_count() const;
    int64_t get_edge_count(int64_t id) const;
    const unordered_map<int64_t,coord_t>& get_flank_map() const;
    coord_t get_flank_coord(const string& name) const;
    coord_t get_flank_coord(int64_t id) const;
    void construct_named_flank_map(unordered_map<string,interval_t>& flank_map) const;

    const HeteroNode& get_node(int64_t id) const;
    const HeteroNode& get_node(const string& name) const;
    void get_sequence(const string& name, string& result) const;
    void get_sequence(int64_t id, string& result) const;
    size_t get_sequence_size(int64_t id) const;
    size_t get_sequence_size(const string& name) const;

    void for_each_neighbor_of_type(int64_t id, char type, const function<void(int64_t id)>& f) const;
    void for_each_neighbor_of_type(int64_t id, char type, const function<void(int64_t id, float w)>& f) const;

    void for_each_sample(const function<void(const string& name, int64_t id)>& f) const;
    bool has_sample(const string& sample_name) const;
    void for_each_read(const function<void(const string& name, int64_t id)>& f) const;
    void for_each_read(const function<void(const string& name, const BinarySequence<uint64_t>& sequence)>& f) const;
    void for_each_read_id(const function<void(int64_t id)>& f) const;
    void for_each_path(const function<void(const string& name, int64_t id)>& f) const;


    void for_each_read_of_sample(const string& sample_name, const function<void(const string& name, int64_t id)>& f) const;
    void for_each_read_of_sample(int64_t sample_id, const function<void(int64_t read_id)>& f) const;
    void for_each_read_of_path(int64_t path_id, const function<void(int64_t id)>& f) const;
    void for_each_variant_of_path(int64_t path_id, const function<void(int64_t id)>& f) const;

    // TODO: refactor and handle read->sample get_ functions for when Transmap was compressed
    string get_sample_of_read(const string& read_name) const;
    void get_sample_of_read(const string& read_name, string& result) const;
    void get_sample_of_read(int64_t read_id, string& result) const;
    int64_t get_sample_of_read(int64_t read_id) const;
    void for_each_sample_of_read(const string& read_name, const function<void(const string& name, int64_t id)>& f) const;
    void for_each_sample_of_read(const int64_t& read_id, const function<void(int64_t id)>& f) const;
    void for_each_sample_of_path(const string& path_name, const function<void(const string& name, int64_t id)>& f) const;
    void for_each_sample_of_path(int64_t id, const function<void(const string& name, int64_t id)>& f) const;

    void for_each_path_of_sample(const string& sample_name, const function<void(const string& name, int64_t id)>& f) const;
    void for_each_path_of_sample(int64_t sample_id, const function<void(const string& name, int64_t id)>& f) const;
    void for_each_path_of_read(int64_t read_id, const function<void(int64_t path_id)>& f) const;
    void for_each_path_of_read(int64_t read_id, const function<void(const string& path_name, int64_t path_id)>& f) const;

    void for_each_phased_variant_of_sample(const string& sample_name, const function<void(const string& name, int64_t id, bool phase)>& f) const;
    void for_each_phased_variant_of_sample(int64_t sample_id, const function<void(const string& name, int64_t id, bool phase)>& f) const;

    void for_each_read_to_path_edge(const function<void(int64_t read_id, int64_t path_id, float weight)>& f) const;

    void for_node_in_bfs(
            const string& start_name,
            float min_edge_weight,
            const function<bool(const HeteroNode& node)>& criteria,
            const function<void(const HeteroNode& node, int64_t id)>& f) const;

    void for_edge_in_bfs(
            const string& start_name,
            float min_edge_weight,
            const function<bool(const HeteroNode& node)>& criteria,
            const function<void(const HeteroNode& a, int64_t a_id, const HeteroNode& b, int64_t b_id)>& f) const;

    /**
     * @return (n_paths,last_path_id)
     */
    pair<int64_t,int64_t> get_n_paths_of_read(int64_t read_id) const;

    /// Writing
    void write_edge_info_to_csv(path output_path, const VariantGraph& variant_graph, bool use_sample_id = false) const;

    /// Clearing
    void clear_non_samples();

    /// Copying/subsetting
    void extract_sample_as_transmap(const string& sample_name, TransMap& result);

    /// Reversible modification
    // Restructure the transmap by duplicating any paths/haps that are used by multiple samples, so that the resulting
    // transmap is essentially a collection of independent sample->read->hap mappings. For this process to be reversible
    // it requires that we keep track of which parent path corresponded to which duplicated child path. We don't
    // store this in the Transmap members because it is rarely used and would be a waste.
    // !! Currently UNTESTED w.r.t. any compression methods !!
    void detangle_sample_paths(unordered_map<string, string> &hapmap);

    void retangle_sample_paths(const unordered_map<string, string> &hapmap);

    /// Compressing

    /**
     * Used just for debugging.
     */
    bool are_edges_distinct() const;

    /**
     * Sorts the neighbors of every node by `(type,id)`. This order is the same for every set of neighbors.
     */
    void sort_adjacency_lists();

    /**
     * Updates the `first_of_type` index of `graph`, assuming that `sort_adjacency_lists()` has already been called.
     * See `VectorHeteroGraph.update_first_of_type()` for details.
     */
    void update_first_of_type();

    void clear_present_haps_edges();

    /**
     * Splits the read-path graph into its connected components (which are at least as many as the connected components
     * of the sample-path graph). Stores in object variable `partitioned_samples` the samples whose reads were assigned
     * to 2 transmaps.
     *
     * Remark: the procedure assumes that the adjacencies of every node are already sorted in an order that is the same
     * for every node.
     *
     * @param maps output array, with one transmap per connected component; every transmap contains the corresponding
     * samples.
     */
    void partition(vector<TransMap>& maps);

    /**
     * @return a transmap with the following connected components (id, n_reads, n_paths, n_samples, n_edges):
     * 0, 2, 1, 2, 9
     * 1, 2, 1, 2, 9
     * 2, 3, 1, 2, 12
     * 3, 2, 1, 1, 8
     * One sample is partitioned into 3 components.
     */
    static TransMap partition_get_test_transmap();

    /**
     * Two reads are considered identical iff they belong to the same sample and they connect to the same haplotypes
     * with the same weights (possibly after quantization). The procedure collapses all identical reads onto a single
     * node and sums edge weights.
     *
     * Remark: samples in `solved_samples` are skipped.
     *
     * Remark: the procedure assumes that the adjacencies of every node are already sorted in an order that is the same
     * for every node.
     *
     * @param weight_quantum if nonzero, read-haplotype weights are divided by this and floored before being compared
     * exactly.
     */
    void compress_reads(float weight_quantum, bool verbose = false);

    /**
     * Two samples are considered identical iff there is a bijection between equivalent reads: only one sample per
     * equivalence class is kept, and edge weights are summed. A sample is contained in another sample iff there is an
     * injection between equivalent reads: every contained sample such that all its reads have only one assignment is
     * collapsed into a container sample, and edge weights are summed.
     *
     * Remark: there could be multiple haplotypes per sample, even though every read in the sample is assigned to
     * exactly one haplotype.
     *
     * Remark: samples in `solved_samples` are skipped.
     *
     * Remark: the procedure assumes that the adjacencies of every node are already sorted in an order that is the same
     * for every node.
     *
     * Remark: the procedure sets object variables `sample_to_sample, read_ids, cluster_ids`.
     *
     * @param weight_quantum if nonzero, read-haplotype weights are divided by this and floored before being compared
     * exactly.
     */
    void compress_samples(float weight_quantum = 0);

    /**
     * Adds every weight of every record in `weights[from_first..from_last]` to a distinct record in
     * `weights[to_first..to_last]` with the same cluster ID.
     *
     * @param cluster_ids one cluster ID per read;
     * @param weights one array of neighbor weights per read;
     * @param used temporary space.
     */
    void compress_samples_update_weights(int64_t from_first, int64_t from_last, int64_t to_first, int64_t to_last, vector<vector<float>>& weights, vector<bool>& used);

    /**
     * Reintroduces a node for every compressed sample using object variables `sample_to_sample, read_ids, cluster_ids`.
     *
     * Remark: decompressed samples are connected to reads of the remaining samples, so the transmap after decompression
     * breaks the assumption that every read is connected to exactly one sample.
     *
     * Remark: this procedure can be called only once, since all object variables it uses are cleared by it.
     *
     * @param used temporary space.
     */
    void decompress_samples();

    /**
     * Adds to object variable `present_haps` all the mandatory haplotypes, and to `present_edges` the corresponding
     * read-hap edges.
     *
     * @return the number of mandatory haplotypes.
     */
    int64_t get_mandatory_haplotypes();

    /**
     * Removes globally-equivalent and globally-contained haplotypes.
     *
     * Remark: the procedure assumes that the adjacencies of every node are already sorted in an order that is the same
     * for every node.
     *
     * @param weight_quantum if nonzero, haplotype-read weights are divided by this and floored before being compared
     * exactly.
     */
    void compress_haplotypes_global(float weight_quantum = 0);

    /**
     * @param neighbors every row is assumed to be sorted;
     * @return true iff every element in `neighbors[from]` occurs in `neighbors[to]` with smaller or equal weight.
     */
    static bool is_haplotype_contained(int64_t from, int64_t to, const vector<vector<int64_t>>& neighbors, const vector<vector<float>>& weights);

    /**
     * @return true iff there is a read-hap weight that is not smaller than `n_weight/d_weight`.
     */
    bool has_large_weight(float n_weight, float d_weight);

    /**
     * Removes edges to dominated haplotypes per-sample. The procedure assumes that the objective function has the form
     * $nWeight \cdot \sum_i h_i + dWeight \cdot \sum_i d_i$ where every $d_i$ is non-negative and $nWeight,dWeight$
     * are non-negative.
     *
     * Remark: the procedure assumes that the adjacencies of every node are already sorted in an order that is the same
     * for every node.
     *
     * @param weight_quantum if nonzero, haplotype-read weights are divided by this and floored before being compared
     * exactly.
     */
    void compress_haplotypes_local(float n_weight, float d_weight, float weight_quantum = 0);

    /**
     * @param neighbors every row is assumed to be sorted;
     * @return true iff every element in `neighbors[from]` occurs in `neighbors[to]` with a weight that is smaller by at
     * least `delta`.
     */
    static bool is_haplotype_dominated(float delta, int64_t from, int64_t to, const vector<vector<int64_t>>& neighbors, const vector<vector<float>>& weights);

    /**
     * Assigns some samples to one or two haplotypes and simplifies the transmap accordingly. The procedure assumes that
     * the objective function has the form $nWeight \cdot \sum_i h_i + dWeight \cdot \sum_i d_i$ where every $d_i$ is
     * non-negative and $nWeight,dWeight$ are non-negative.
     *
     * Remark: the procedure adds to object variable `present_haps` all the haplotypes that were used to solve a sample,
     * and it adds to `present_edges` all the edges that were used to solve a sample (all the other edges of a solved
     * sample are removed from the transmap). The procedure also sets object variable `solved_samples`.
     *
     * Remark: the procedure assumes that the adjacencies of every node are already sorted in an order that is the same
     * for every node.
     *
     * @param weight_quantum if nonzero, haplotype-read weights are divided by this and floored before being compared
     * exactly.
     */
    void solve_easy_samples(float n_weight, float d_weight, float weight_quantum);

    /**
     * @return a transmap to test `solve_easy_samples()`: 2 out of 3 samples are easy, and `solve_easy_samples()`
     * should solve 1 one-hap sample and 1 two-hap sample; it should remove 12 edges and it should set to one 3
     * haplotypes and 6 edges.
     */
    static TransMap solve_easy_samples_get_test_transmap(float n_weight, float d_weight);

    /**
     * @return `weight` if `weight_quantum=0`; otherwise, `weight` rounded to the nearest multiple of `weight_quantum`.
     */
    static float get_edge_weight(float weight, float weight_quantum);
};

/// WARNING: DOES NOT COMPARE EDGE WEIGHTS, ONLY COMPARES EDGE PRESENCE/ABSENCE
bool operator==(const TransMap& a, const TransMap& b);

}
