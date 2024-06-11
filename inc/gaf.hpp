#pragma once

#include "TransitiveMap.hpp"
#include "VariantGraph.hpp"
#include "Alignment.hpp"

#include <unordered_map>
#include <functional>
#include <stdexcept>
#include <cstdlib>
#include <utility>
#include <fstream>
#include <string>
#include <array>

using std::unordered_map;
using std::runtime_error;
using std::function;
using std::ifstream;
using std::string;
using std::array;
using std::pair;
using std::pair;


namespace sv_merge{


class GafAlignment: public Alignment{
    vector <pair <string,bool> > path;
    vector <CigarTuple> cigar;
    string query_name;
    string ref_name;
    int32_t query_length;
    int32_t query_start;
    int32_t query_stop;
    int32_t path_length;
    int32_t path_start;
    int32_t path_stop;
    int32_t n_match;
    int32_t alignment_length;
    int32_t map_quality;
    bool reversal;
    bool primary;

public:
    /// Modifying
    void set_path(const vector<pair<string,bool> >& p);
    void set_path(const string& path);
    void set_query_name(const string& query_name);
    void set_ref_name(const string& ref_name);
    void set_query_length(int32_t query_length);
    void set_query_start(int32_t query_start);
    void set_query_stop(int32_t query_stop);
    void set_path_length(int32_t path_length);
    void set_path_start(int32_t path_start);
    void set_path_stop(int32_t path_stop);
    void set_reversal(bool reversal);
    void set_n_match(int32_t n);
    void set_alignment_length(int32_t length);
    void set_map_quality(int32_t q);
    void set_is_primary(bool p);
    void add_tag(const string& tag);
    void clear_cigar();

    /// Accessing
    [[nodiscard]] string get_path_string() const;
    [[nodiscard]] string get_query_name() const;
    [[nodiscard]] string get_ref_name() const;
    [[nodiscard]] int32_t get_query_length() const override;
    [[nodiscard]] int32_t get_query_start() const override;
    [[nodiscard]] int32_t get_query_stop() const;
    [[nodiscard]] int32_t get_path_length() const;
    [[nodiscard]] int32_t get_path_start() const;
    [[nodiscard]] int32_t get_path_stop() const;
    [[nodiscard]] const pair<string,bool>& get_path_step(int32_t index) const;
    [[nodiscard]] int32_t get_n_match() const;
    [[nodiscard]] int32_t get_alignment_length() const;
    [[nodiscard]] int32_t get_map_quality() const;
    [[nodiscard]] bool is_primary() const override;
    [[nodiscard]] bool is_supplementary() const override;
    [[nodiscard]] bool is_reverse() const override;
    void for_step_in_path(const function<void(const string& step_name, bool is_reverse)>& f) const;
    [[nodiscard]] const pair<string,bool>& get_step_of_path(size_t index) const;
    [[nodiscard]] const vector<pair<string,bool> >& get_path() const;

    /// Helper
    static bool parse_path_reversal_token(char c);
    static void parse_string_as_path(const string& p, vector<pair<string,bool> >& result);
    void load_cigar(const string& cigar_string);

    /// Iterators and Alignment API fulfillment
    void for_each_cigar_interval(bool unclip_coords, const function<void(const CigarInterval& cigar)>& f) override;
    void for_each_cigar_tuple(const function<void(const CigarTuple& cigar)>& f) override;
    void get_query_sequence(string& result) override;
    void get_query_sequence(string& result, int32_t start, int32_t stop) override;
    void get_qualities(vector<uint8_t>& result) override;
    void get_tag_as_string(const string& name, string& result, bool allow_missing=false) const override;
    void get_query_name(string& result) const override;
    [[nodiscard]] int32_t get_ref_start() const override;
    [[nodiscard]] int32_t get_ref_stop() const override;
    [[nodiscard]] bool is_unmapped() const override;
};


bool parse_reversal_token(const string& token);

void for_alignment_in_gaf(const path& gaf_path, const function<void(GafAlignment& alignment)>& f);

void for_alignment_in_gaf(const path& gaf_path, const function<void(Alignment& alignment)>& f);

class AlignmentSummary{
public:
    int32_t start;
    int32_t stop;
    float n_match;
    float n_mismatch;
    float n_insert;
    float n_delete;

    AlignmentSummary();
    AlignmentSummary(int32_t start, int32_t stop);
    void update(const CigarInterval& cigar_interval, bool is_ref);
    float compute_identity() const;
};


class GafSummary{
public:
    unordered_map<string,int32_t> node_lengths;
//    unordered_map<string,int32_t> query_lengths;

    unordered_map <string, vector<AlignmentSummary> > ref_summaries;
    unordered_map <string, vector<AlignmentSummary> > query_summaries;

    // Paths observed in alignments
    unordered_map <string, vector <vector <pair<string,bool> > > > query_paths;

    // Transmap holds all sequence-related data (lengths, flanks)
    const TransMap& transmap;
    bool apply_flanks;

    /**
     * For testing purposes only, does not properly initialize
     */
    GafSummary(const TransMap& transmap);

    /**
     * This constructor needs some objects that can inform the GafSummary of all the lengths of the graph nodes/queries.
     * TODO: For future applications, consider passing a GFA path, or simply a map of node/query lengths.
     * @param variant_graph a fully built variant graph that corresponds to the graph aligned as target in the GAF
     * @param trans_map a transitive map that contains only the reads that will be aligned to the graph.
     * @param apply_flanks if true, use the info stored in TransMap to subset the portion of the query that is evaluated
     */
    GafSummary(const VariantGraph& variant_graph, const TransMap& trans_map, bool apply_flanks);

    /**
     * Update alignment stats, for a given ref node
     * @param node_name ref node to be assigned this cigar operation
     * @param c cigar interval obtained from iterating an alignment
     * @param insert a new alignment should start with this operation (separate alignments will be overlap-resolved)
     */
    void update_node(const string& node_name, const CigarInterval& c, bool insert);

    /**
     * Update alignment stats, for a given query sequence
     * @param query_name query to be assigned this cigar operation
     * @param c cigar interval obtained from iterating an alignment
     * @param insert a new alignment should start with this operation (separate alignments will be overlap-resolved)
     */
    void update_query(const string& query_name, const CigarInterval& c, bool insert);

    void for_each_ref_summary(const function<void(const string& name, int32_t length, float identity, float coverage)>& f) const;
    void for_each_query_summary(const function<void(const string& name, int32_t length, float identity, float coverage)>& f) const;

    void resolve_all_overlaps();
    static void resolve_overlaps(vector<AlignmentSummary>& alignments);

    /**
     * Accumulate data from a GAF that can summarize the alignment in terms of the query sequences and the graph. If
     * `apply_flanks` was previously specified, then this accumulation of stats on the queries will only pertain to
     * regions inside the inner flank boundaries of the queries, as provided by the transmap given at construction time.
     * @param gaf_path GAF to be parsed
     */
    void compute(const path& gaf_path);

private:
    /**
     * Accumulate data from a GAF that can summarize the alignment in terms of the query sequences and the graph
     * @param gaf_path GAF to be parsed
     */
    void compute_without_flanks(const path& gaf_path);

    /**
     * Accumulate data from a GAF that can summarize the alignment in terms of the query sequences and the graph. In this
     * variant of the compute() fn, we exclude flanking regions from the accumulation of cigar operations on the query
     * @param gaf_path GAF to be parsed
     */
    void compute_with_flanks(const path& gaf_path);
};


class HalfInterval{
public:
    size_t id;
    int32_t position;
    bool is_start;

    HalfInterval(size_t id, int32_t position, bool is_start);
};


}
