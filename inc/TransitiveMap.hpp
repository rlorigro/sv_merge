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

using std::unordered_map;
using std::unordered_set;
using std::to_string;
using std::vector;
using std::string;
using std::pair;
using std::byte;


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
     * For every original sample, the sample it was compressed into. Used only by procedures `compress()` and
     * `decompress()`.
     */
    unordered_map<string,string> sample_to_compressed_sample;

public:
    static const string sample_node_name;
    static const string read_node_name;
    static const string path_node_name;
    static const string variant_node_name;

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
    void for_each_read(const function<void(const string& name, int64_t id)>& f) const;
    void for_each_read(const function<void(const string& name, const BinarySequence<uint64_t>& sequence)>& f) const;
    void for_each_read_id(const function<void(int64_t id)>& f) const;
    void for_each_path(const function<void(const string& name, int64_t id)>& f) const;

    void get_read_sample(const string& read_name, string& result) const;
    void get_read_sample(int64_t read_id, string& result) const;
    int64_t get_read_sample(int64_t read_id) const;

    void for_each_read_of_sample(const string& sample_name, const function<void(const string& name, int64_t id)>& f) const;
    void for_each_read_of_sample(int64_t sample_id, const function<void(int64_t read_id)>& f) const;
    void for_each_read_of_path(int64_t path_id, const function<void(int64_t id)>& f) const;
    void for_each_variant_of_path(int64_t path_id, const function<void(int64_t id)>& f) const;

    string get_sample_of_read(const string& read_name) const;
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

    /// Writing
    void write_edge_info_to_csv(path output_path, const VariantGraph& variant_graph) const;

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
     * Two reads are considered identical iff they connect to the same haplotypes with the same weights (possibly
     * after quantization). The procedure collapses all identical reads onto a single node, which becomes connected to
     * all the samples of the reads in its equivalence class (breaking the assumption that a read is connected to just
     * one sample).
     *
     * Two samples are considered identical iff their reads belong to the same set of read clusters (regardless of
     * how many reads in each sample belong to each cluster). Only one sample per equivalence class is kept, and the
     * mapping is stored in object variable `sample_to_compressed_sample`.
     *
     * Remark: the procedure assumes that haplotypes are already distinct from previous steps, and it does not try to
     * compress them.
     *
     * Remark: the procedure sets object variables `n_reads, n_read_clusters, n_samples, n_sample_clusters`.
     *
     * @param weight_quantum if nonzero, read-haplotype weights are divided by this and floored before being compared
     * exactly;
     * @param mode the weight of every read-hap edge is set to the max (mode=0), min (mode=1), sum (mode=2) or avg
     * (mode=3) of all the edges that were collapsed onto it.
     */
    void compress(float weight_quantum, uint64_t mode);

    /**
     * Reintroduces a node for every compressed sample using object variable `sample_to_compressed_sample` (which is
     * cleared when this procedure returns). The rest of the graph remains compressed.
     */
    void decompress_samples();
};


/// WARNING: DOES NOT COMPARE EDGE WEIGHTS, ONLY COMPARES EDGE PRESENCE/ABSENCE
bool operator==(const TransMap& a, const TransMap& b);

}
