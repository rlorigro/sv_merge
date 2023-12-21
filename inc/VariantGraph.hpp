#pragma once

#include "Filesystem.hpp"
#include "Region.hpp"
#include "VcfReader.hpp"
#include "bdsg/hash_graph.hpp"
#include "misc.hpp"

using ghc::filesystem::path;
using sv_merge::Region;
using sv_merge::VcfRecord;
using sv_merge::VcfReader;

#include <string>
#include <iterator>
#include <unordered_map>
#include <unordered_set>

using std::string;
using std::unordered_map;
using std::unordered_set;
using std::sort;
using std::unique;
using std::min;
using std::max;
using bdsg::nid_t;
using bdsg::HashGraph;
using bdsg::handle_t;
using bdsg::edge_t;
using bdsg::path_handle_t;
using bdsg::step_handle_t;


namespace sv_merge {

/**
 * Converts a set of VCF records into a bidirected graph, converts the graph back to a set of VCF records, and
 * serializes the graph. Replacement VCF records, BND VCF records, and BND VCF records with inserted sequence are
 * supported. The bidirected graph might have several nodes with edges only on one side (e.g. for BND records).
 */
class VariantGraph {
public:
    /**
     * Bidirected graph. Users are assumed to interact with it directly.
     */
    HashGraph graph;

    /**
     * @param chromosomes map id -> sequence;
     * @param tandem_track a sorted list of non-overlapping intervals per chromosome, with format `[x..y)`.
     */
    VariantGraph(const unordered_map<string,string>& chromosomes, const unordered_map<string,vector<interval_t>>& tandem_track);

    /**
     * Given a list of VCF records from the same chromosome, the procedure builds a corresponding bidirected graph and
     * keeps track of which records support each node and edge.
     *
     * Remark: if a complex SV consists of multiple BND records in multiple chromosomes, only the BNDs in one chromosome
     * are represented in the graph.
     *
     * Remark: the procedure can be called multiple times on different instances of `records`.
     *
     * @param records a list of VCF records, not necessarily sorted; we assume that every POS belongs to the same
     * chromosome, and that every record is biallelic;
     * @param flank_length ensure that an interval of this length, with no overlap to the tandem track, is present
     * before the leftmost breakpoint (after the rightmost breakpoint) of every chromosome in the graph;
     * @param deallocate_ins TRUE: the procedure deallocates the ALT sequence of every non-symbolic INS in `records`.
     */
    void build(vector<VcfRecord>& records, uint16_t flank_length, bool deallocate_ins);

    /**
     * Serializes the bidirected graph. Paths, if any, are not saved.
     */
    void to_gfa(const path& gfa_path) const;

    /**
     * Consider the set of VCF records that were used to build the graph. Every such VCF record corresponds to a
     * sequence of oriented edges. The procedure writes to a file every VCF record R such that there is a (possibly
     * circular) path in `graph` that supports all the edges of R in their original orientation, or all of them in their
     * reverse orientation.
     *
     * Remark: the output file is in the same order as `vcf_records` and it does not have a header.
     * Remark: for every two mated BND records, the procedure outputs just one of them.
     */
    void to_vcf(const path& vcf_path);

    /**
     * Writes every path in the graph as a distinct VCF record, of replacement type or of BND type with inserted
     * sequence.
     *
     * Remark: the procedure only outputs paths in the graph that begin and end at a node that belongs to a chromosome
     * (i.e. not an INS node). Circular paths are linearized according to `graph`, and this choice might not satisfy
     * the constraint.
     *
     * Remark: the output file is not sorted and it does not have a header.
     */
    void to_vcf_paths(const path& vcf_path);

    /**
     * Resets `out` to the set of all IDs of nodes that have edges only on one side.
     */
    void get_ids_of_dangling_nodes(unordered_set<nid_t> out) const;


private:
    const unordered_map<string,string>& chromosomes;
    const uint8_t n_chromosomes;
    const unordered_map<string,vector<interval_t>>& tandem_track;
    vector<VcfRecord> vcf_records;
    uint32_t n_vcf_records;
    unordered_set<string> bnd_ids;  // ID field of every BND record

    /**
     * For every chromosome, every first position of a chunk induced by breakpoints (zero-based).
     */
    unordered_map<string,vector<uint32_t>> chunk_first;
    unordered_map<string,vector<handle_t>> node_handles;

    /**
     * For every node of the graph that belongs to a chromosome, its chromosome and first position (zero-based).
     */
    unordered_map<nid_t,pair<string,uint32_t>> node_to_chromosome;

    /**
     * All the new nodes coming from INS, REPLACEMENT, and BND.
     *
     * Remark: for simplicity the ALT sequence of every VCF record is treated as distinct, even though it may be
     * identical to the ALT sequence of other VCF records.
     */
    vector<handle_t> insertion_handles;

    /**
     * Maps every edge in the graph to the VCF records that support it.
     */
    unordered_map<edge_t,vector<uint32_t>> edge_to_vcf_record;

    /**
     * For every VCF record, its sequence of edges in `graph`.
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
    uint16_t find_closest(const string& chromosome, uint32_t position) const;

    /**
     * Adds `(from,to)` to `graph` if not already present, and adds `record` to the edge in `edge_to_vcf_record`.
     */
    void add_edge(const handle_t& from, const handle_t& to, uint32_t record_id);

    /**
     * Removes duplicates
     */
    static void sort_and_compact_positions(vector<uint32_t>& positions);

    /**
     * Looks for `query` in `edges`, and sets `flags[i]=1` (resp. =2; =3) if `query` equals `edges[i]` in the forward
     * (resp. reverse; both forward and reverse) orientation.
     *
     * Remark: this is implemented as a linear scan, since `edges` is assumed to be short.
     */
    void mark_edge(const edge_t& query, const vector<edge_t>& edges, vector<uint8_t>& flags) const;

    /**
     * Remark: for simplicity this is implemented as a linear scan. It could be made faster.
     *
     * @param pos zero-based;
     * @param tandem_track a sorted list of non-overlapping intervals per chromosome, with format `[x..y)`;
     * @return the smallest `y` such that `[x..y]` has no intersection with `tandem_track`, `y-x+1=flank_length`, and
     * `x>pos`; if no such interval exists, returns UINT32_MAX.
     */
    uint32_t get_flank_boundary_right(const string& chromosome_id, uint32_t pos, uint16_t flank_length);

    /**
     * Remark: for simplicity this is implemented as a linear scan. It could be made faster.
     *
     * @param pos zero-based;
     * @param tandem_track a sorted list of non-overlapping intervals per chromosome, with format `[x..y)`;
     * @return the largest `x` such that `[x..y]` has no intersection with `tandem_track`, `y-x+1=flank_length`, and
     * `y<pos`; if no such interval exists, returns UINT32_MAX.
     */
    uint32_t get_flank_boundary_left(const string& chromosome_id, uint32_t pos, uint16_t flank_length);

    /**
     * @param interval1,interval2 with format [x..y).
     */
    static bool intersect(const interval_t& interval1, const interval_t& interval2);
};

}