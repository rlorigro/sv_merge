#pragma once

#include "VectorHeteroGraph.hpp"
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
    unordered_map<int64_t,string> sequences;

    // Pains me to add yet another map but here it is
    unordered_map<int64_t,coord_t> sequence_flanks;

    // Pains me to add yet another map but here it is
    unordered_map<int64_t,bool> sequence_reversals;

    string sample_node_name;
    string read_node_name;
    string path_node_name;
    string variant_node_name;

    /**
     * Data structures used by compression/decompression procedures.
     */
    unordered_map<string,pair<string,int64_t>> sample_to_sample;

public:
    TransMap();

    /// Building
    void reserve_nodes(size_t n);
    void reserve_edges(size_t n);
    void reserve_sequences(size_t n);
    void add_sample(const string& name);
    void add_read(const string& name);
    void add_flank_coord(const string& name, int32_t start, int32_t stop);
    void add_read(const string& name, const string& sequence);
    void add_read_with_move(string& name, string& sequence);
    void add_read_with_move(string& name, string& sequence, bool reversal);
    void add_path(const string& name);
    void add_variant(const string& name);
    void add_edge(const string& a, const string& b);
    void add_edge(int64_t a, int64_t b, float weight);
    void add_edge(const string& a, const string& b, float weight);

    void remove_edge(int64_t a, int64_t b);
    void remove_node(int64_t id);

    /// Accessing
    bool empty() const;
    int64_t get_id(const string& name) const;
    pair<bool,int64_t> try_get_id(const string& name) const;
    pair<bool,float> try_get_edge_weight(int64_t id_a, int64_t id_b) const;
    bool has_edge(int64_t a, int64_t b) const;
    bool has_node(const string& name) const;
    int64_t get_node_count() const;
    int64_t get_edge_count(int64_t id) const;
    const unordered_map<int64_t,coord_t>& get_flank_map() const;
    coord_t get_flank_coord(const string& name) const;
    coord_t get_flank_coord(int64_t id) const;
    void construct_named_flank_map(unordered_map<string,interval_t>& flank_map) const;

    const HeteroNode& get_node(int64_t id) const;
    const HeteroNode& get_node(const string& name) const;
    const string& get_sequence(const string& name) const;
    const string& get_sequence(int64_t id) const;

    void for_each_neighbor_of_type(int64_t id, char type, const function<void(int64_t id)>& f) const;

    void for_each_sample(const function<void(const string& name, int64_t id)>& f) const;
    bool contains_sample(const string& sample_name) const;
    void for_each_read(const function<void(const string& name, int64_t id)>& f) const;
    void for_each_read(const function<void(const string& name, const string& sequence)>& f) const;
    void for_each_read_id(const function<void(int64_t id)>& f) const;
    void for_each_path(const function<void(const string& name, int64_t id)>& f) const;

    void get_read_sample(const string& read_name, string& result) const;
    void get_read_sample(int64_t read_id, string& result) const;

    void for_each_read_of_sample(const string& sample_name, const function<void(const string& name, int64_t id)>& f) const;
    void for_each_read_of_sample(int64_t sample_id, const function<void(int64_t read_id)>& f) const;
    void for_each_read_of_path(int64_t path_id, const function<void(int64_t id)>& f) const;
    void for_each_variant_of_path(int64_t path_id, const function<void(int64_t id)>& f) const;

    string get_sample_of_read(const string& read_name) const;
    void for_each_sample_of_read(const string& read_name, const function<void(const string& name, int64_t id)>& f) const;
    void for_each_sample_of_read(const int64_t& read_id, const function<void(int64_t id)>& f) const;
    void for_each_sample_of_path(const string& path_name, const function<void(const string& name, int64_t id)>& f) const;

    void for_each_path_of_sample(const string& sample_name, const function<void(const string& name, int64_t id)>& f) const;
    void for_each_path_of_sample(int64_t sample_id, const function<void(const string& name, int64_t id)>& f) const;
    void for_each_path_of_read(int64_t read_id, const function<void(int64_t path_id)>& f) const;

    void for_each_phased_variant_of_sample(const string& sample_name, const function<void(const string& name, int64_t id, bool phase)>& f) const;
    void for_each_phased_variant_of_sample(int64_t sample_id, const function<void(const string& name, int64_t id, bool phase)>& f) const;

    void for_each_read_to_path_edge(const function<void(int64_t read_id, int64_t path_id, float weight)>& f) const;

    void for_node_in_bfs(
            const string& start_name,
            float min_edge_weight,
            const function<bool(const HeteroNode& node)>& criteria,
            const function<void(const HeteroNode& node, int64_t id)>& f) const;

    int64_t get_n_reads() const;
    int64_t get_n_samples() const;
    int64_t get_n_paths() const;
    int64_t get_n_edges() const;


    /// Writing
    void write_edge_info_to_csv(path output_path, const VariantGraph& variant_graph) const;

    /// Clearing
    void clear_non_samples();

    /// Compressing

    /**
     * Splits the read-path graph into its connected components (which are at least as many as the connected components
     * of the sample-path graph).
     *
     * @param maps output array, with one transmap per connected component; every transmap contains the corresponding
     * samples;
     * @param partitioned_samples output array, lists the samples whose reads were assigned to 2 transmaps.
     */
    void partition(vector<TransMap>& maps, vector<string>& partitioned_samples) const;

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
     * node.
     *
     * @param weight_quantum if nonzero, read-haplotype weights are divided by this and floored before being compared
     * exactly;
     * @param mode the weight of every read-hap edge is set to the max (mode=0), min (mode=1), sum (mode=2) or avg
     * (mode=3) of all the edges that were collapsed onto it;
     * @param sort_edges if false, the procedure assumes that the adjacencies of every node are already sorted in an
     * order that is the same for every node.
     */
    void compress_reads(float weight_quantum, uint64_t mode, bool sort_edges = true, bool verbose = false);

    /**
     * Two samples are considered identical iff there is a bijection between equivalent reads: only one sample per
     * equivalence class is kept, and edge weights are summed. A sample is contained in another sample iff there is an
     * injection between equivalent reads: every contained sample such that all its reads have only one assignment is
     * collapsed into a container sample.
     *
     * Remark: the procedure sets object variables `sample_to_identical_sample, sample_to_container_sample,
     * node_to_cluster`.
     *
     * @param sort_edges if false, the procedure assumes that the adjacencies of every node are already sorted in an
     * order that is the same for every node.
     */
    void compress_samples(float weight_quantum, bool sort_edges = true);

    /**
     * Adds every weight of every record in `weights[from_first..from_last]` to a distinct record in
     * `weights[to_first..to_last]` with the same cluster ID.
     */
    void compress_samples_update_weights(int64_t from_first, int64_t from_last, int64_t to_first, int64_t to_last, vector<int64_t>& cluster_ids, vector<vector<float>>& weights, vector<bool>& used);

    /**
     * Reintroduces a node for every compressed sample using object variable `sample_to_compressed_sample`. The rest of
     * the graph remains compressed.
     *
     * Remark: this procedure can be called only once, since `sample_to_compressed_sample` is cleared by it.
     */
    void decompress_samples();

    bool are_edges_distinct() const;
};


}
