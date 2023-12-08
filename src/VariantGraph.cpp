#include "VariantGraph.hpp"

#include <string_view>

using std::min;
using std::max;
using std::find;
using std::distance;

namespace sv_merge {

const uint8_t VariantGraph::N_CHROMOSOMES = 25;  // 22+X+Y+M
const string VariantGraph::CHR_STR = "chr";
const string VariantGraph::X_STR_PRIME = "X";
const string VariantGraph::Y_STR_PRIME = "Y";
const string VariantGraph::M_STR_PRIME = "M";
const string VariantGraph::MT_STR_PRIME = "MT";
const uint8_t VariantGraph::CHR_STR_LENGTH = CHR_STR.length();
const string VariantGraph::X_STR = CHR_STR+X_STR_PRIME;
const string VariantGraph::Y_STR = CHR_STR+Y_STR_PRIME;
const string VariantGraph::M_STR = CHR_STR+M_STR_PRIME;
const string VariantGraph::MT_STR = CHR_STR+MT_STR_PRIME;


VariantGraph::VariantGraph(const vector<string>& chromosomes):
        chromosomes(chromosomes)
{
    this->chunk_first.reserve(N_CHROMOSOMES);
    this->node_handles.reserve(N_CHROMOSOMES);
}


uint16_t VariantGraph::chrom2id(const string& chrom) {
    const uint8_t length = chrom.length();
    if (length>=CHR_STR_LENGTH && chrom.starts_with(CHR_STR)) {
        if (chrom==X_STR) return 23;
        else if (chrom==Y_STR) return 24;
        else if (chrom==M_STR || chrom==MT_STR) return 25;
        else if (length<=CHR_STR_LENGTH+2) return (uint16_t)stoul(chrom.substr(CHR_STR_LENGTH,length-CHR_STR_LENGTH));
        else return UINT16_MAX;
    }
    else {
        if (chrom==X_STR_PRIME) return 23;
        else if (chrom==Y_STR_PRIME) return 24;
        else if (chrom==M_STR_PRIME || chrom==MT_STR_PRIME) return 25;
        else if (length<=2) return (uint16_t)stoul(chrom);
        else return UINT16_MAX;
    }
}


void VariantGraph::sort_and_compact_positions(vector<uint32_t>& positions) {
    sort(positions.begin(),positions.end());
    auto iterator = unique(positions.begin(),positions.end());
    positions.resize(distance(positions.begin(),iterator));
}


uint16_t VariantGraph::find_closest(uint8_t chromosome_id, uint32_t position) {
    const vector<uint32_t> positions = chunk_first.at(chromosome_id);
    auto iterator = lower_bound(positions.begin(),positions.end(),position);
    if (iterator==positions.end()) return positions.size()-1;
    uint16_t p = distance(positions.begin(),iterator);
    if (p==0 || positions.at(p)==position) return p;
    return position-positions.at(p-1)<=positions.at(p)-position?p-1:p;
}


void VariantGraph::add_edge(handle_t from, handle_t to, const VcfRecord& record) {
    if (graph.has_edge(from,to)) edge_to_vcf_record[graph.edge_handle(from,to)].push_back(record.clone());
    else {
        graph.create_edge(from, to);
        edge_to_vcf_record[graph.edge_handle(from, to)] = {record.clone()};
    }
}

/**
 *
 * in theory support a variable flank different for left and right, and moreover for each BND we should have a separate flank length.
 *
 */
void VariantGraph::build(vector<VcfRecord>& records, uint16_t flank_length, bool deallocate_ins) {
    if (records.empty()) return;
    const uint16_t CHROMOSOME_ID = chrom2id(chromosomes[chrom2id(records.at(0).chrom)]);
    uint8_t orientation, orientation_trans;
    uint16_t i, j, n_positions, n_handles, trans_chr_id;
    uint32_t p, q, first_pos, last_pos, trans_pos;
    string tmp_buffer;
    pair<uint32_t, uint32_t> tmp_pair;
    handle_t handle_from, handle_to, reversed_handle;

    for (i=0; i<N_CHROMOSOMES; i++) { chunk_first.at(i).clear(); node_handles.at(i).clear(); }
    insertion_handles.clear(); insertion_to_vcf_record.clear(); edge_to_vcf_record.clear();
    graph.clear();

    // Collecting every distinct first position of a reference node. Building all INS nodes.
    vector<uint32_t>& first_positions = chunk_first.at(CHROMOSOME_ID);
    for (auto& record: records) {
        record.get_reference_coordinates(false,tmp_pair);
        if (tmp_pair.first==UINT32_MAX || tmp_pair.second==UINT32_MAX) continue;
        if (tmp_pair.first!=tmp_pair.second) { first_positions.push_back(tmp_pair.first); first_positions.push_back(tmp_pair.second); }
        else {
            if (record.sv_type==VcfReader::TYPE_INSERTION && !record.is_symbolic) {
                first_positions.push_back(tmp_pair.first);
                handle_t handle_ins = graph.create_handle(record.alt.substr(1));
                insertion_handles.push_back(handle_ins);
                insertion_to_vcf_record[handle_ins]=record;
                if (deallocate_ins) record.alt.clear();
            }
            else if (record.sv_type==VcfReader::TYPE_BREAKEND) {
                orientation=record.get_breakend_orientation_cis();
                if (orientation==1) first_positions.push_back(tmp_pair.first+1);
                else if (orientation==2) first_positions.push_back(tmp_pair.first);
                trans_pos=record.get_breakend_pos();
                if (trans_pos==UINT32_MAX) {
                    // NOP: virtual telomeric breakend.
                }
                else {
                    trans_pos--; // Zero-based
                    record.get_breakend_chromosome(tmp_buffer);
                    trans_chr_id=chrom2id(tmp_buffer);
                    orientation=record.get_breakend_orientation_trans();
                    if (orientation==1) chunk_first.at(trans_chr_id).push_back(trans_pos+1);
                    else if (orientation==2) chunk_first.at(trans_chr_id).push_back(trans_pos);
                }
            }
        }
    }
    for (i=0; i<N_CHROMOSOMES; i++) {
        if (chunk_first.at(i).size()<=1) continue;
        sort_and_compact_positions(chunk_first.at(i));
    }

    // Building all the non-INS nodes and all reference edges
    for (i=0; i<N_CHROMOSOMES; i++) {
        first_positions = chunk_first.at(i);
        n_positions = first_positions.size();
        if (n_positions == 0) continue;
        vector<handle_t>& handles = node_handles.at(i);
        handles.clear();
        string& chrom_sequence = chromosomes.at(i);
        p = first_positions.at(0);
        first_pos = max(p - flank_length, 0);
        if (p != first_pos) handles.push_back(graph.create_handle(chrom_sequence.substr(first_pos, p - first_pos)));
        for (j = 1; j < n_positions; j++) {
            p = first_positions.at(j - 1);
            handles.push_back(graph.create_handle(chrom_sequence.substr(p, first_positions.at(j) - p)));
        }
        p = first_positions.at(n_positions - 1);
        last_pos = min(p + flank_length, chrom_sequence.length());
        if (p != last_pos) handles.push_back(graph.create_handle(chrom_sequence.substr(p, last_pos - p)));
        n_handles=handles.size();
        for (j=1; j<n_handles; j++) graph.create_edge(handles.at(j-1),handles.at(j));
    }

    // Building non-reference edges
    const vector<handle_t>& handles = node_handles.at(CHROMOSOME_ID);
    i=0;
    for (auto& record: records) {
        record.get_reference_coordinates(false,tmp_pair);
        if (tmp_pair.first==UINT32_MAX || tmp_pair.second==UINT32_MAX) continue;
        if (tmp_pair.first!=tmp_pair.second) {
            p=find_closest(CHROMOSOME_ID,tmp_pair.first);
            q=find_closest(CHROMOSOME_ID,tmp_pair.second);
            if (record.sv_type==VcfReader::TYPE_DELETION) add_edge(handles.at(p), handles.at(q+1), record);
            else if (record.sv_type==VcfReader::TYPE_DUPLICATION) add_edge(handles.at(q),handles.at(p+1),record);
            else if (record.sv_type==VcfReader::TYPE_INVERSION) {
                reversed_handle = graph.flip(handles.at(q));
                add_edge(handles.at(p), reversed_handle, record);
                reversed_handle = graph.flip(handles.at(p + 1));
                add_edge(reversed_handle, handles.at(q + 1), record);
            }
        }
        else {
            if (record.sv_type==VcfReader::TYPE_INSERTION && !record.is_symbolic) {
                p = find_closest(CHROMOSOME_ID, tmp_pair.first);
                add_edge(handles.at(p), insertion_handles.at(i), record);
                add_edge(insertion_handles.at(i), handles.at(p + 1), record);
                i++;
            }
            else if (record.sv_type==VcfReader::TYPE_BREAKEND) {
                orientation=record.get_breakend_orientation_cis();
                p=find_closest(CHROMOSOME_ID,tmp_pair.first);
                if (orientation==1) handle_from=handles.at(p);
                else if (orientation==2) handle_from=graph.flip(handles.at(p+1));
                trans_pos=record.get_breakend_pos();
                if (trans_pos==UINT32_MAX) {
                    // NOP: virtual telomeric breakend.
                }
                else {
                    record.get_breakend_chromosome(tmp_buffer);
                    p=find_closest(chrom2id(tmp_buffer),trans_pos-1);
                    orientation=record.get_breakend_orientation_trans();
                    if (orientation==1) handle_to=graph.flip(node_handles.at(trans_chr_id).at(p));
                    else if (orientation==2) handle_to=node_handles.at(trans_chr_id).at(p+1);
                }
                add_edge(handle_from,handle_to,record);
            }
        }
    }
    graph.optimize(true);

    // Clearing temporary space
    for (i=0; i<N_CHROMOSOMES; i++) { chunk_first.at(i).clear(); node_handles.at(i).clear(); }
    insertion_handles.clear();
}



void VariantGraph::clean() {

    graph.for_each_path_handle(const Iteratee &iteratee) {

    }


}


}