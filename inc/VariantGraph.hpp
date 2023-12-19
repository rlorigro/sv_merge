#pragma once

#include "Filesystem.hpp"
#include "Region.hpp"
#include "VcfReader.hpp"
#include "bdsg/hash_graph.hpp"

using ghc::filesystem::path;
using sv_merge::Region;
using sv_merge::VcfRecord;
using sv_merge::VcfReader;

#include <algorithm>
#include <string>
#include <iterator>
#include <unordered_map>

using std::string;
using std::unordered_map;
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
 *
 *
 */
class VariantGraph {
public:
    /**
     * Variant graph. Users are assumed to interact with it directly.
     */
    HashGraph graph;

    /**
     * @param chromosomes map id->sequence
     */
    VariantGraph(const unordered_map<string,string>& chromosomes);

    /**
     * Given a list of VCF records, the procedure builds a corresponding bidirected graph and keeps track of which
     * records support each node and edge.
     *
     * Remark: the procedure supports replacement VCF records, and breakend VCF records with inserted sequence.
     * Remark: the procedure can be called multiple times on different instances of `records`.
     *
     * @param records a list of VCF records, not necessarily sorted; we assume that every POS belongs to the same
     * chromosome, and that every record is biallelic;
     * @param flank_length maximum amount of sequence before the leftmost breakpoint (after the rightmost breakpoint)
     * of a chromosome;
     * @param deallocate_ins TRUE: the procedure deallocates the ALT sequence of every non-symbolic INS in `records`.
     */
    void build(vector<VcfRecord>& records, unordered_map<string,vector<interval_t>>& tandem_track, uint16_t flank_length, bool deallocate_ins);

    /**
     * Serializes `graph`. Paths are not saved.
     */
    void to_gfa(const path& gfa_path) const;

    /**
     * Consider the set of VCF records that were used to build `graph`. Every such VCF record corresponds to a
     * sequence of oriented edges. The procedure writes to a file every VCF record R such that: there is a (possibly
     * circular) path in `graph` that supports all the edges of R in their original orientation, or all of them in their
     * reverse orientation.
     *
     * Remark: the output file is in the same order as `vcf_records` and it does not have a header.
     * Remark: for every two mated breakend records, the procedure outputs just one of them.
     */
    void to_vcf(const path& vcf_path);

    /**
     * Writes every path in `graph` as a distinct VCF record, of replacement type or of breakend type with inserted
     * sequence. Circular paths are linearized according to `graph`.
     *
     * Remark: the output file is not sorted and it does not have a header.
     */
    void to_vcf_paths(const path& vcf_path);

private:
    const unordered_map<string,string>& chromosomes;
    const uint8_t n_chromosomes;
    vector<VcfRecord> vcf_records;
    uint32_t n_vcf_records;
    unordered_set<string> bnd_ids;  // ID field of every BND record

    /**
     * For every chromosome, every first position of a chunk induced by breakpoints (zero-based).
     */
    unordered_map<string,vector<uint32_t>> chunk_first;
    unordered_map<string,vector<handle_t>> node_handles;

    /**
     * For every node of the graph, its chromosome and its first position (zero-based).
     */
    unordered_map<nid_t,pair<string,uint32_t>> node_to_chromosome;    //   <-------------- !!!!!!!!!!!!

    /**
     * All the new nodes coming from INS and REPLACEMENT, and their corresponding VCF records.
     *
     * Remark: the ALT sequence of every VCF record is treated as distinct for simplicity, even though it may be
     * identical to the ALT sequence of other VCF records.
     */
    vector<handle_t> insertion_handles;
    vector<uint32_t> insertion_to_vcf_record;

    /**
     * Maps every edge of the graph to the VCF records that support it, if any.
     */
    unordered_map<edge_t,vector<uint32_t>> edge_to_vcf_record;

    /**
     * For every VCF record, its corresponding insertion or replacement nodes, and its edges, in `graph`.
     */
    vector<vector<handle_t>> vcf_record_to_insertion;
    vector<vector<edge_t>> vcf_record_to_edge;

    /**
     * Temporary space for printing VCFs.
     */
    vector<bool> printed;
    vector<vector<uint8_t>> flags;

    /**
     * @return the index of a closest element to `position` in `chunk_first`.
     */
    uint16_t find_closest(const string& chromosome, const uint32_t position) const;

    /**
     * Adds `(from,to)` to `graph` if not already present, and adds `record` to the edge in `edge_to_vcf_record`.
     */
    void add_edge(const handle_t& from, const handle_t& to, const uint32_t record_id);

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
};

}