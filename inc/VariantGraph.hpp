#pragma once

#include "Filesystem.hpp"
#include "Region.hpp"
#include "VcfReader.hpp"
#include "bdsg/hash_graph.hpp"
#include "misc.hpp"

using sv_merge::Region;
using sv_merge::VcfRecord;
using sv_merge::VcfReader;
using ghc::filesystem::path;
using bdsg::nid_t;
using bdsg::HashGraph;
using bdsg::handle_t;
using bdsg::edge_t;
using bdsg::path_handle_t;

#include <string>
#include <iterator>
#include <unordered_map>
#include <unordered_set>
#include <ostream>

using std::string;
using std::unordered_map;
using std::unordered_set;
using std::ofstream;

namespace sv_merge {

/**
 * Converts a set of VCF records into a bidirected graph, converts the graph back to a set of VCF records, and
 * serializes the graph. The class supports replacement VCF records, BND VCF records, BND VCF records with inserted
 * sequence, and "single BND" records. The bidirected graph might contain more than two nodes with edges only on one
 * side (e.g. with BND records).
 */
class VariantGraph {
public:
    /**
     * VCF records given in input to graph construction.
     * Some VCF records might not have been used for building the graph.
     */
    vector<VcfRecord> vcf_records;

    /**
     * Bidirected graph. Users are assumed to interact with it directly.
     */
    HashGraph graph;

    /**
     * @param chromosomes map id -> sequence;
     * @param tandem_track map chromosome -> sorted list of zero-based, non-overlapping intervals, with format `[x..y)`.
     * A chromosome might not appear in the map. A chromosome might appear in the map and its list of intervals might be
     * empty. The tandem track is used for adding flanking regions with enough non-tandem sequence to the graph, which
     * is useful for seeding graph alignments: see `build()` for more details.
     */
    VariantGraph(const unordered_map<string,string>& chromosomes, const unordered_map<string,vector<interval_t>>& tandem_track = {});

    /**
     * Given a list of VCF records from the same chromosome, the procedure builds a corresponding bidirected graph and
     * keeps track of which records support which non-reference edge.
     *
     * Remark: CNVs are modeled as DUPs even though they are abundance claims, not adjacency claims (i.e. the copies can
     * be anywhere in the genome). Copy number estimates in the GT fields are ignored.
     *
     * Remark: if a complex SV consists of multiple BND records in multiple chromosomes, only the BNDs with one end in
     * the chromosome of the CHROM field are represented in the graph.
     *
     * Remark: every BND is loaded, even though it might be part of an event that affects fewer bases than those
     * conventionally assigned to SVs.
     *
     * Remark: the procedure can be called multiple times on different instances of `records`.
     *
     * @param records a list of VCF records, not necessarily sorted. We assume that every POS belongs to the same
     * chromosome, and that every record is biallelic. If an event starts at the first position of a chromosome, POS
     * must be zero and REF must start with N (this is the usual VCF convention; the VCF spec seems to deviate from this
     * convention at the first position of a chromosome). Duplicate records, if any, do not create duplicate edges in
     * the graph; instead, they are all mapped to the same edge. There can be records that are identical in everything
     * except their ID, e.g. if they were created by different callers (although such records would have been merged by
     * `bcftools merge` and their IDs would have been concatenated).
     * @param flank_length ensures that an interval of this length, with no overlap to the tandem track, is present
     * before the leftmost breakpoint (after the rightmost breakpoint) of every chromosome in the graph;
     * @param interior_flank_length consider two consecutive breakpoints A, B on the same chromosome; the procedure
     * computes the shortest window to the right of A (resp. to the left of B) that contains an interval of this length
     * with no overlap to the tandem track; if the windows of A and of B do not overlap, two disconnected nodes are
     * created; otherwise, a single node [A..B) is created; this is useful to avoid creating long nodes when the
     * breakpoints of a chromosome form clusters that are far away from one another; this parameter should be set to a
     * small multiple of average read length;
     * @param x,y (zero-based) if these are specified, and if `x` is smaller (resp. `y` is greater) than all the
     * breakpoints in `records` on the chromosome in the CHROM field: compute flanks based on `x` and `y` rather than on
     * the leftmost/rightmost breakpoints in `records`;
     * @param deallocate_ref_alt TRUE: the procedure deallocates every REF and ALT field in `records`; this can be
     * useful for reducing space, since these fields might contain explicit DNA sequences;
     * @param callers caller names (lowercase), used just for printing statistics; caller names must occur in the ID
     * field of a VCF record in order to be considered.
     */
    void build(vector<VcfRecord>& records, int32_t flank_length, int32_t interior_flank_length = INT32_MAX, int32_t x = INT32_MAX, int32_t y = INT32_MAX, bool deallocate_ref_alt = false, const vector<string>& callers = {});

    /**
     * If `p<q` (zero-based), the procedure builds a trivial graph that contains one node for string `chromosome[p..q)`
     * and one node for each of its flanking sequences (if they exist). If `p=q`, the central node is replaced by a
     * single edge.
     *
     * @param flank_length ensures that an interval of this length, with no overlap to the tandem track, is present
     * before `p` and at or after `q`.
     */
    void build(const string& chromosome, int32_t p, int32_t q, int32_t flank_length);

    /**
     * @return TRUE iff a graph built from `records` would contain at least one non-reference edge.
     * This can be useful for deciding which `build()` function to call.
     */
    bool would_graph_be_nontrivial(vector<VcfRecord>& records);

    /**
     * Serializes the bidirected graph, including paths, if any.
     */
    void to_gfa(const path& gfa_path) const;

    /**
     * The procedure sets only variable `graph`. GFA paths are loaded, if any.
     *
     * @param gfa_path all node IDs are assumed to be distinct;
     * @return lexicographically sorted list of node IDs; the integer IDs in `graph` are positions in this list.
     */
    vector<string> load_gfa(const path& gfa_path);

    /**
     * Consider the set of VCF records that were used for building `graph`. Every such VCF record R corresponds to a
     * sequence of non-reference edges. The procedure writes to a file every R such that there is a (possibly circular)
     * path P in `graph` that traverses the entire sequence of non-reference edges of R, or its reverse. The procedure
     * writes to another file every R for which this does not hold.
     *
     * Remark: in the current implementation, DELs that remove a prefix or suffix of a chromosome create no edge in the
     * graph, so they cannot be supported by a path and they are not printed in output. This could be solved by creating
     * an artificial source (resp. sink) node before the beginning (resp. after the end) of every chromosome.
     *
     * Remark: the output VCF is in the same order as `vcf_records`, it does not have a header, and it contains every
     * field in `vcf_records`.
     *
     * Remark: for every pair of mated BND records, the procedure outputs just one of them.
     *
     * @param print_all_records TRUE=disregard paths and print every VCF record to `supported` (including those that
     * were not used for building the graph), and no VCF record to `unsupported`;
     * @param callers caller names (lowercase), used just for printing statistics; caller names must occur in the ID
     * field of a VCF record to be counted.
     */
    void print_supported_vcf_records(ofstream& supported, ofstream& unsupported, bool print_all_records, const vector<string>& callers = {});

    /**
     * Writes every path handle in `graph` as a distinct VCF record:
     * 1. a replacement record if the first and last node of the path are on the same chromosome, in the same
     * orientation, and with compatible positions;
     * 2. a BND record with inserted sequence if the first and last node are not on the same chromosome, or have
     * different orientations, or have incompatible positions;
     * 3. a "single BND" record with inserted sequence if either the first or the last node (but not both) do not belong
     * to a chromosome (i.e. they are new sequence);
     * 4. paths where neither the first nor the last node belong to a chromosome are ignored.
     *
     * Remark: circular paths are simply linearized according to `graph`. This could be improved by trying to find a
     * linearization that is compatible with categories 1-3 above.
     *
     * Remark: the output file is not sorted and it does not have a header.
     */
    void paths_to_vcf_records(ofstream& outstream);

    /**
     * If `mode=TRUE`, makes `out` the set of all IDs of nodes that have neighbors only on one side. Otherwise, prints
     * information about dangling nodes to STDERR.
     */
    void get_dangling_nodes(bool mode, unordered_set<nid_t>& out) const;

    size_t get_n_dangling_nodes() const;

    /**
     * @return TRUE iff the node has neighbors only on one side.
     */
    bool is_dangling_node(const handle_t& node_handle) const;

    /**
     * @return TRUE iff the node has neighbors only on one side.
     */
    bool is_dangling_node(const nid_t& node_id) const;

    /**
     * @return TRUE iff the node has neighbors only on one side, and belongs to the chromosome where all the POS
     * fields of `vcf_records` are located.
     */
    bool is_flanking_node(const handle_t& node_handle) const;

    /**
     * @return TRUE iff the node has neighbors only on one side, and belongs to the chromosome where all the POS
     * fields of `vcf_records` are located.
     */
    bool is_flanking_node(const nid_t& node_id) const;

    /**
     * @return the number of edges between two reference nodes that belong to different chromosomes. This might be
     * smaller than the number of inter-chromosomal BND records, if some BNDs have an inserted sequence.
     */
    size_t get_n_interchromosomal_edges() const;

    /**
     * @return the number of edges between a reference node and a new node not in the reference. This might be greater
     * than the number of INS and replacement records, if some BNDs have an inserted sequence.
     */
    size_t get_n_ins_edges() const;

    bool is_reference_node(const nid_t& node_id) const;

    /**
     * @param node_handle in any orientation.
     */
    bool is_reference_node(const handle_t& node_handle) const;

    /**
     * @param edge in any orientation.
     */
    bool is_reference_edge(const edge_t& edge);

    /**
     * Every VCF record corresponds to a sequence of non-reference edges. The procedure makes `out` the set of all the
     * VCF records whose sequence of non-reference edges coincides with the given sequence of edges or its reverse.
     *
     * @param edges each edge can be represented in any orientation;
     * @param out VCF records are in the order in which they appear in `vcf_records`.
     */
    void get_vcf_records_with_edges(const vector<edge_t>& edges, vector<VcfRecord>& out);

    /**
     * Iterates over every VCF record and its sequence of non-reference edges. VCF records that do not contribute any
     * non-reference edge to `graph` are not iterated.
     *
     * @param id unique integer identifier of `record`.
     */
    void for_each_vcf_record(const function<void(size_t id, const vector<edge_t>& edges_of_the_record, const VcfRecord& record)>& callback);

    /**
     * Given a path P, the procedure iterates over every VCF record R (and its non-reference edges) that is supported by
     * P, i.e. such that P traverses the sequence of non-reference edges of R, or its reverse.
     *
     * @param path a sequence of pairs `(node_id, is_reverse)`; node IDs are assumed to come from the set of node IDs in
     * `graph`;
     * @param id a unique integer identifier of `record`.
     */
    void for_each_vcf_record(const vector<pair<string,bool>>& path, const function<void(size_t id, const vector<edge_t>& edges_of_the_record, const VcfRecord& record)>& callback);

    /**
     * Prints to STDERR the number of non-reference edges supported by X VCF records, and the number of VCF records
     * creating X non-reference edges.
     */
    void print_edge_histograms() const;

    /**
     * Prints to STDERR the sorted sequence of all distinct INS, DUP, DEL, and INV lengths (binned) in `vcf_records`,
     * and the number of their occurrences.
     */
    void print_vcf_records_stats(uint8_t bin_size, const vector<string>& callers={}) const;

    /**
     * Prints to `output_path` one line per node `v` of `graph`, with format `seq_1,seq_2,...,seq_k`, where `seq_*` are
     * all and only the strings that can be built with a search tree of height at most `max_steps` from `v`. The list is
     * sorted lexicographically. The lines of the file are sorted lexicographically as well.
     *
     * Remark: strings in `seq_*` are not canonized.
     */
    void print_graph_signature(size_t max_steps, const path& output_path) const;

    /**
     * Sets `edge_to_vcf_record` and `vcf_record_to_edge`.
     *
     * @param map a list of `(edge,vcfID)` pairs, where `edge` can be in any orientation.
     */
    void load_edge_record_map(const vector<pair<edge_t,size_t>>& map, size_t n_vcf_records);

    /**
     * Erases `node_to_chromosome`.
     */
    void node_to_chromosome_clear();

    /**
     * Inserts a mapping into `node_to_chromosome`.
     */
    void node_to_chromosome_insert(const string& node_label, const vector<string>& node_labels, const string& chromosome, int32_t position);

    /**
     * Erases `insertion_handles_set`.
     */
    void insertion_handles_set_clear();

    /**
     * Inserts elements into `insertion_handles_set`.
     */
    void insertion_handles_set_insert(const string& node_label, const vector<string>& node_labels);

private:
    /**
     * GFA constants
     */
    static const char GFA_SEPARATOR;
    static const char GFA_NODE_CHAR;
    static const char GFA_LINK_CHAR;
    static const char GFA_PATH_CHAR;
    static const string GFA_OVERLAP_FIELD;
    static const char GFA_PATH_SEPARATOR;
    static const char GFA_FWD_CHAR;
    static const char GFA_REV_CHAR;
    static const char GFA_LINE_END;
    static const char GAF_FWD_CHAR;
    static const char GAF_REV_CHAR;
    static const uint64_t STREAMSIZE_MAX;

    const unordered_map<string,string>& chromosomes;
    const uint8_t n_chromosomes;
    const unordered_map<string,vector<interval_t>>& tandem_track;
    size_t n_vcf_records;
    string main_chromosome;  // The CHROM field of every VCF record
    int32_t main_chromosome_length;
    unordered_set<string> bnd_ids;  // ID field of every BND record

    /**
     * For every chromosome: every first position of a chunk induced by breakpoints (zero-based).
     */
    unordered_map<string,vector<int32_t>> chunk_first_raw;
    unordered_map<string,vector<int32_t>> chunk_first;
    unordered_map<string,vector<handle_t>> node_handles;

    /**
     * For every node of the graph that belongs to a chromosome: its chromosome and first position (zero-based).
     */
    unordered_map<nid_t,pair<string,int32_t>> node_to_chromosome;

    /**
     * All the new nodes coming from INS, REPLACEMENT, and BND records.
     *
     * Remark: for simplicity the ALT sequence of every VCF record is treated as distinct, even though it may be
     * identical to the ALT sequence of other VCF records.
     */
    vector<handle_t> insertion_handles;
    unordered_set<handle_t> insertion_handles_set;  // Same content as above

    /**
     * Maps every non-reference edge of the graph (in canonical form) to the VCF records that support it.
     */
    unordered_map<edge_t,vector<size_t>> edge_to_vcf_record;

    /**
     * For every VCF record: its sequence of non-reference edges in `graph` (in canonical form).
     */
    vector<vector<edge_t>> vcf_record_to_edge;

    /**
     * Reused temporary space
     */
    vector<bool> printed, initialized;
    vector<vector<size_t>> flags;
    interval_t tmp_interval;
    vector<edge_t> tmp_edges;

    /**
     * @return the index of a closest element to `position` in `chunk_first`.
     */
    inline size_t find_closest(const string& chromosome, int32_t position) const;

    /**
     * Adds `(from,to)` to `graph` if not already present, and connects `record_id` to the edge in `edge_to_vcf_record`
     * and `vcf_record_to_edge`.
     */
    void add_nonreference_edge(const handle_t& from, const handle_t& to, size_t record_id);

    /**
     * Removes duplicated positions from a list.
     */
    static inline void sort_and_compact_positions(vector<int32_t>& positions);

    /**
     * Stores `path_encoding` (assumed to be a valid path in GFA format) in `graph`.
     *
     * @param node_ids string IDs used in the GFA file;
     * @param buffer temporary space.
     */
    path_handle_t load_gfa_path(const string& path_encoding, const vector<string>& node_ids, const string& path_name, string& buffer);

    /**
     * Stores `path_encoding` (assumed to be a valid path in GAF format) in `graph`.
     *
     * @param path_encoding node IDs are assumed to come from the set of node IDs in `graph`.
     */
    path_handle_t load_gaf_path(const string& path_encoding, const string& path_name, string& buffer);

    /**
     * Stores `path_encoding` (assumed to be a valid path) in `graph`.
     *
     * @param path a sequence of pairs `(node_id, is_reverse)`; node IDs are assumed to come from the set of node IDs in
     * `graph`.
     */
    path_handle_t load_gaf_path(vector<pair<string,bool>>& path, const string& path_name);

    /**
     * Destroys all paths in `graph` by iterating over them explicitly.
     */
    void destroy_paths();

    /**
     * Looks up `query` in `edges`, and sets `flags[i]=rank` if the canonized `query` equals `edges[i]`.
     *
     * Remark: for simplicity this is implemented as a linear scan, since `edges` is assumed to be short.
     *
     * @param query in any orientation.
     */
    void mark_edge(const edge_t& query, size_t rank, const vector<edge_t>& edges, vector<size_t>& flags) const;

    /**
     * Every VCF record corresponds to a sequence of non-reference edges. The procedure makes `out` the set of all
     * VCF records whose sequence of non-reference edges (or its reverse) is identical to (if `identical=true`) or
     * contained in (if `identical=false`) the given sequence of edges.
     *
     * @param edges each edge can be represented in any orientation;
     * @param out positions in `vcf_records`, sorted.
     */
    void get_vcf_records_with_edges_impl(const vector<edge_t>& edges, bool identical, vector<size_t>& out);

    /**
     * Remark: for simplicity this is implemented as a linear scan. It could be made faster.
     *
     * @param pos zero-based;
     * @param tandem_track a sorted list of non-overlapping intervals per chromosome, with format `[x..y)`;
     * @return the smallest `y` such that `[x..y]` has no intersection with `tandem_track`, `y-x+1=flank_length`, and
     * `x>=pos`; returns `INT32_MAX` if no such interval exists.
     */
    int32_t get_flank_boundary_right(const string& chromosome_id, int32_t pos, int32_t flank_length);

    /**
     * Remark: for simplicity this is implemented as a linear scan. It could be made faster.
     *
     * @param pos zero-based; can be equal to `chromosome_length`;
     * @param tandem_track a sorted list of non-overlapping intervals per chromosome, with format `[x..y)`;
     * @return the largest `x` such that `[x..y]` has no intersection with `tandem_track`, `y-x+1=flank_length`, and
     * `y<pos`; returns `INT32_MAX` if no such interval exists.
     */
    int32_t get_flank_boundary_left(const string& chromosome_id, int32_t pos, int32_t flank_length);

    /**
     * Remark: this is a recursive procedure that creates new string objects. It should be implemented more efficiently
     * with a stack of characters and without recursion.
     */
    void print_signature_impl(const handle_t& handle, const string& path, size_t steps_performed, size_t max_steps, vector<string>& strings) const;

    /**
     * @param interval1,interval2 format [x..y).
     */
    static inline bool intersect(const interval_t& interval1, const interval_t& interval2);

    static void print_sv_lengths(vector<size_t>& lengths, uint8_t bin_size, const string& prefix) ;

    static void print_bnds(vector<pair<string,string>>& bnds) ;

    static void increment_caller_count(const VcfRecord& record, const vector<string>& callers, vector<vector<size_t>>& caller_count, uint8_t column);
};

}