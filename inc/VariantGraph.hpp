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
     * Bidirected graph. Users are assumed to interact with it directly.
     */
    HashGraph graph;

    /**
     * @param chromosomes map id -> sequence;
     * @param tandem_track a sorted list of zero-based, non-overlapping intervals per chromosome, with format `[x..y)`.
     */
    VariantGraph(const unordered_map<string,string>& chromosomes, const unordered_map<string,vector<interval_t>>& tandem_track);

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
     * the graph; instead, they are all mapped to the same edge.
     * @param flank_length ensures that an interval of this length, with no overlap to the tandem track, is present
     * before the leftmost breakpoint (after the rightmost breakpoint) of every chromosome in the graph; if two
     * consecutive breakpoints are farther away than `flank_length`, two disconnected nodes are created;
     * @param deallocate_ref_alt TRUE: the procedure deallocates every REF and ALT field in `records`; this can be
     * useful for reducing space, since these fields might contain explicit DNA sequences;
     * @param callers caller names (lowercase), used just for printing statistics; caller names must occur in the ID
     * field of a VCF record to be counted.
     */
    void build(vector<VcfRecord>& records, uint16_t flank_length, bool deallocate_ref_alt, const vector<string>& callers={});

    /**
     * Serializes the bidirected graph, including paths, if any.
     */
    void to_gfa(const path& gfa_path) const;

    /**
     * The procedure sets only variable `graph`. GFA paths are loaded, if any.
     *
     * @param gfa_path all node IDs are assumed to be distinct;
     * @return sorted list of node IDs; the integer IDs in `graph` are positions in this list.
     */
    vector<string> load_gfa(const path& gfa_path);

    /**
     * Consider the set of VCF records that were used for building `graph`. Every such VCF record corresponds to a
     * sequence of oriented edges. The procedure writes to a file every VCF record R such that there is a (possibly
     * circular) path in `graph` that supports all the edges of R in their original orientation, or all of them in their
     * reverse orientation.
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
     * @param print_all_records TRUE=disregard paths and print every VCF record in input (including those that were not
     * used for building the graph);
     * @param callers caller names (lowercase), used just for printing statistics; caller names must occur in the ID
     * field of a VCF record to be counted.
     */
    void print_supported_vcf_records(ofstream& outstream, bool print_all_records, const vector<string>& callers= {});

    /**
     * Writes every path in `graph` as a distinct VCF record:
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
     * If `mode=TRUE`, makes `out` the set of all IDs of nodes that have edges only on one side. Otherwise, prints
     * information about dangling nodes to STDERR.
     */
    void get_dangling_nodes(bool mode, unordered_set<nid_t>& out) const;

    uint32_t get_n_dangling_nodes() const;

    /**
     * @return the number of edges between two reference nodes that belong to different chromosomes. This might be
     * smaller than the number of inter-chromosomal BND records, if some BNDs have an inserted sequence.
     */
    uint32_t get_n_interchromosomal_edges() const;

    /**
     * @return the number of edges between a reference node and a new node not in the reference. This might be greater
     * than the number of INS and replacement records, if some BNDs have an inserted sequence.
     */
    uint32_t get_n_ins_edges() const;

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
    void print_graph_signature(uint16_t max_steps, const path& output_path) const;

    /**
     * Sets `edge_to_vcf_record` and `vcf_record_to_edge`.
     *
     * @param map a list of (edge,vcfID) pairs.
     */
    void load_edge_record_map(const vector<pair<edge_t,uint32_t>>& map, uint32_t n_vcf_records);

    /**
     * Erases `node_to_chromosome`.
     */
    void node_to_chromosome_clear();

    /**
     * Inserts a mapping into `node_to_chromosome`.
     */
    void node_to_chromosome_insert(const string& node_id, const vector<string> node_ids, const string& chromosome, uint32_t position);

private:
    /**
     * GFA constants
     */
    static const char GFA_SEPARATOR;
    static const char GFA_NODE_CHAR;
    static const char GFA_LINK_CHAR;
    static const char GFA_PATH_CHAR;
    static const char GFA_OVERLAP_FIELD;
    static const char GFA_PATH_SEPARATOR;
    static const char GFA_FWD_CHAR;
    static const char GFA_REV_CHAR;
    static const char GFA_LINE_END;
    static const uint64_t STREAMSIZE_MAX;

    const unordered_map<string,string>& chromosomes;
    const uint8_t n_chromosomes;
    const unordered_map<string,vector<interval_t>>& tandem_track;
    vector<VcfRecord> vcf_records;
    uint32_t n_vcf_records;
    unordered_set<string> bnd_ids;  // ID field of every BND record

    /**
     * For every chromosome: every first position of a chunk induced by breakpoints (zero-based).
     */
    unordered_map<string,vector<uint32_t>> chunk_first_raw;
    unordered_map<string,vector<uint32_t>> chunk_first;
    unordered_map<string,vector<handle_t>> node_handles;

    /**
     * For every node of the graph that belongs to a chromosome: its chromosome and first position (zero-based).
     */
    unordered_map<nid_t,pair<string,uint32_t>> node_to_chromosome;

    /**
     * All the new nodes coming from INS, REPLACEMENT, and BND records.
     *
     * Remark: for simplicity the ALT sequence of every VCF record is treated as distinct, even though it may be
     * identical to the ALT sequence of other VCF records.
     */
    vector<handle_t> insertion_handles;

    /**
     * Maps every non-reference edge of the graph to the VCF records that support it.
     */
    unordered_map<edge_t,vector<uint32_t>> edge_to_vcf_record;

    /**
     * For every VCF record: its sequence of non-reference edges in `graph`.
     */
    vector<vector<edge_t>> vcf_record_to_edge;

    /**
     * Reused temporary space
     */
    vector<bool> printed;
    vector<vector<uint8_t>> flags;
    interval_t tmp_interval;

    /**
     * @return the index of a closest element to `position` in `chunk_first`.
     */
    inline uint16_t find_closest(const string& chromosome, uint32_t position) const;

    /**
     * Adds `(from,to)` to `graph` if not already present, and adds `record` to the edge in `edge_to_vcf_record`.
     */
    void add_nonreference_edge(const handle_t& from, const handle_t& to, uint32_t record_id);

    /**
     * Removes duplicated positions from a list.
     */
    static inline void sort_and_compact_positions(vector<uint32_t>& positions);

    /**
     * Transforms `path_encoding` (a path in GFA format) to a path in `graph`.
     *
     * @param node_ids string IDs used in the GFA file;
     * @param buffer temporary space.
     */
    void load_gfa_path(const string& path_encoding, const vector<string>& node_ids, const string& path_name, string& buffer);

    /**
     * Destroys all paths in `graph` by iterating over them explicitly.
     */
    void destroy_paths();

    /**
     * Looks for `query` in `edges`, and sets `flags[i]=1` (resp. 2; resp. 3) if `query` equals `edges[i]` in the
     * forward (resp. reverse; resp. both forward and reverse) orientation.
     *
     * Remark: for simplicity this is implemented as a linear scan, since `edges` is assumed to be short.
     */
    void mark_edge(const edge_t& query, const vector<edge_t>& edges, vector<uint8_t>& flags) const;

    /**
     * Remark: for simplicity this is implemented as a linear scan. It could be made faster.
     *
     * @param pos zero-based;
     * @param tandem_track a sorted list of non-overlapping intervals per chromosome, with format `[x..y)`;
     * @return the smallest `y` such that `[x..y]` has no intersection with `tandem_track`, `y-x+1=flank_length`, and
     * `x>pos`; returns `UINT32_MAX` if no such interval exists.
     */
    uint32_t get_flank_boundary_right(const string& chromosome_id, uint32_t pos, uint16_t flank_length);

    /**
     * Remark: for simplicity this is implemented as a linear scan. It could be made faster.
     *
     * @param pos zero-based; can be equal to `chromosome_length`;
     * @param tandem_track a sorted list of non-overlapping intervals per chromosome, with format `[x..y)`;
     * @return the largest `x` such that `[x..y]` has no intersection with `tandem_track`, `y-x+1=flank_length`, and
     * `y<pos`; returns `UINT32_MAX` if no such interval exists.
     */
    uint32_t get_flank_boundary_left(const string& chromosome_id, uint32_t pos, uint16_t flank_length);

    /**
     * Remark: this is a recursive procedure that creates new string objects. It should be implemented more efficiently
     * with a stack of characters and without recursion.
     */
    void print_signature_impl(const handle_t& handle, const string& path, uint16_t steps_performed, uint16_t max_steps, vector<string>& strings) const;

    /**
     * @param interval1,interval2 format [x..y).
     */
    static inline bool intersect(const interval_t& interval1, const interval_t& interval2);

    static void print_sv_lengths(vector<uint32_t>& lengths, uint8_t bin_size, const string& prefix) ;

    static void print_bnds(vector<pair<string,string>>& bnds) ;

    static void increment_caller_count(const VcfRecord& record, const vector<string>& callers, vector<vector<uint32_t>>& caller_count, uint8_t column);
};

}