#pragma once

#include "Filesystem.hpp"
#include "Region.hpp"
#include "VcfReader.hpp"
#include "bdsg/hash_graph.hpp"

using ghc::filesystem::path;
using sv_merge::Region;
using sv_merge::VcfRecord;
using sv_merge::VcfReader;

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


namespace sv_merge {

/**
 *
 *
 */
class VariantGraph {
public:
    /**
     * Basic chromosome constants
     */
    static const uint8_t N_CHROMOSOMES;
    static const string CHR_STR;
    static const uint8_t CHR_STR_LENGTH;
    static const string X_STR;
    static const string X_STR_PRIME;
    static const string Y_STR;
    static const string Y_STR_PRIME;
    static const string M_STR;
    static const string M_STR_PRIME;
    static const string MT_STR;
    static const string MT_STR_PRIME;

    /**
     * @param chromosomes chromosome sequences in canonical order.
     */
    VariantGraph(const vector<string>& chromosomes);

    /**
     * Given an array of VCF records, the procedure builds a corresponding bidirected graph and keeps track of which
     * records support each node and edge.
     *
     * @param records an array of VCF records, not necessarily sorted; we assume that every POS belongs to the same
     * chromosome, and that every record is biallelic;
     * @param flank_length maximum amount of sequence before the leftmost breakpoint (after the rightmost breakpoint)
     * of a chromosome;
     * @param deallocate_ins TRUE: the procedure deallocates the ALT sequence of every non-symbolic INS in `records`.
     */
    void build(vector<VcfRecord>& records, uint16_t flank_length, bool deallocate_ins);  // <-- Add TRF track and compute flanks adaptively for each BND

    /**
     * Removes everything that does not belong to a path
     */
    void clean();


    // create also a VCF replacement allele for every path as a VCF substitution.

    // ^ support VCF substitutions during construction as well, since we will need to repeat this during evaluation with our own substitution paths.


    /**
     * Serializes `graph`.
     */
    void to_gfa(const path& gfa_path);

    /**
     * Expresses as a VCF the set of nodes and edges that are supported by at least one path in `paths`.
     *
     * Remark: this needs to keep reference coordinates.
     */
    void to_vcf(const path& vcf_path);





private:
    /**
     * Chromosome sequences in canonical order.
     */
    const unordered_map<string,string>& chromosomes;

    /**
     * All the new nodes coming from INS, and their corresponding VCF records.
     *
     * Remark: the INS sequence of every VCF record is treated as distinct for simplicity, even though it may be
     * identical to the INS sequence of other VCF records.
     */
    vector<handle_t> insertion_handles;
    unordered_map<handle_t,VcfRecord> insertion_to_vcf_record;

    /**
     * For every chromosome, every first position of a chunk induced by breakpoints (zero-based).
     */
    unordered_map<string,vector<uint32_t>> chunk_first;
    unordered_map<string,vector<handle_t>> node_handles;

    /**
     * Maps every edge of the graph to the VCF records that support it
     */
    unordered_map<edge_t,vector<VcfRecord>> edge_to_vcf_record;

    /**
     * Variant graph
     */
    HashGraph graph;

    /**
     * @return UINT16_MAX for non-standard chromosomes.
     */
    static uint16_t chrom2id(const string& chrom);

    /**
     * @return the index of a closest element to `position` in `chunk_first`.
     */
    uint16_t find_closest(uint8_t chromosome_id, uint32_t position);

    /**
     * Adds `(from,to)` to `graph` if not already present, and associates a copy of `record` to the edge in
     * `edge_to_vcf_record`.
     */
    void add_edge(handle_t from, handle_t to, const VcfRecord& record);

    /**
     * Removes duplicates
     */
    static void sort_and_compact_positions(vector<uint32_t>& positions);
};

}