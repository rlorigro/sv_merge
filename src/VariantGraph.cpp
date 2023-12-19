#include "VariantGraph.hpp"
#include "Sequence.hpp"

#include <cmath>
#include <string_view>
#include <ostream>


// use coord_t everywhere for chr coords from <misc.h>

using std::min;
using std::max;
using std::find;
using std::distance;
using std::ofstream;
using std::to_string;
using sv_merge::reverse_complement;

namespace sv_merge {


VariantGraph::VariantGraph(const unordered_map<string,string>& chromosomes):
        chromosomes(chromosomes),
        n_chromosomes(chromosomes.size())
{ }


void VariantGraph::sort_and_compact_positions(vector<uint32_t>& positions) {
    sort(positions.begin(),positions.end());
    auto iterator = unique(positions.begin(),positions.end());
    positions.resize(distance(positions.begin(),iterator));
}


uint16_t VariantGraph::find_closest(const string& chromosome, const uint32_t position) const {
    const vector<uint32_t> positions = chunk_first.at(chromosome);
    auto iterator = lower_bound(positions.begin(),positions.end(),position);
    if (iterator==positions.end()) return positions.size()-1;
    uint16_t p = distance(positions.begin(),iterator);
    if (p==0 || positions.at(p)==position) return p;
    return position-positions.at(p-1)<=positions.at(p)-position?p-1:p;
}


void VariantGraph::add_edge(const handle_t& from, const handle_t& to, const uint32_t record_id) {
    if (graph.has_edge(from,to)) {
        const edge_t edge = graph.edge_handle(from,to);
        if (!edge_to_vcf_record.contains(edge)) edge_to_vcf_record.at(edge)={record_id};
        else edge_to_vcf_record.at(edge).push_back(record_id);
        vcf_record_to_edge.at(record_id).push_back(edge);
    }
    else {
        graph.create_edge(from, to);
        edge_t edge = graph.edge_handle(from, to);
        edge_to_vcf_record.at(edge)={record_id};
        vcf_record_to_edge.at(record_id).push_back(edge);
    }
}


uint32_t next_pos_outside_tandem(string chromosome, uint32_t pos, unordered_map<string,vector<interval_t>>& tandem_track, uint16_t flank_length) {

}


/**
 *
 * in theory support a variable flank different for left and right, and moreover for each BND we should have a separate flank length.
 *
 *
 *
 */
void VariantGraph::build(vector<VcfRecord>& records, unordered_map<string,vector<interval_t>>& tandem_track, uint16_t flank_length, bool deallocate_ins) {
    graph.clear();
    this->vcf_records=std::move(records);
    this->n_vcf_records=vcf_records.size();
    chunk_first.clear(); chunk_first.reserve(n_chromosomes);
    node_handles.clear(); node_handles.reserve(n_chromosomes);
    insertion_handles.clear(); insertion_to_vcf_record.clear();
    edge_to_vcf_record.clear();
    vcf_record_to_insertion.clear(); vcf_record_to_edge.clear();
    bnd_ids.clear();
    if (vcf_records.empty()) return;

    const string CHROMOSOME_ID = vcf_records.at(0).chrom;
    uint8_t orientation_cis, orientation_trans;
    uint16_t n_positions, n_handles;
    uint32_t i, j, p, q, first_pos, last_pos, trans_pos;
    string tmp_buffer;
    pair<uint32_t, uint32_t> tmp_pair;

    // Collecting every distinct first position of a reference node. Building all INS, BND-INS, and REPLACEMENT nodes.
    vector<uint32_t>& first_positions = chunk_first.at(CHROMOSOME_ID);
    for (i=0; i<n_vcf_records; i++) {
        VcfRecord& record = vcf_records.at(i);
        record.get_reference_coordinates(false,tmp_pair);
        if (tmp_pair.first==UINT32_MAX || tmp_pair.second==UINT32_MAX) continue;
        if (tmp_pair.first!=tmp_pair.second) {
            first_positions.push_back(tmp_pair.first); first_positions.push_back(tmp_pair.second);
            if (record.sv_type==VcfReader::TYPE_REPLACEMENT) {
                const handle_t handle_alt = graph.create_handle(record.alt.substr(1));
                insertion_handles.emplace_back(handle_alt);
                insertion_to_vcf_record.emplace_back(i);
                vcf_record_to_insertion.at(i).emplace_back(handle_alt);
                if (deallocate_ins) { record.ref.clear(); record.alt.clear(); }
            }
        }
        else {
            if (record.sv_type==VcfReader::TYPE_INSERTION && !record.is_symbolic) {
                first_positions.emplace_back(tmp_pair.first);
                const handle_t handle_ins = graph.create_handle(record.alt.substr(1));
                insertion_handles.emplace_back(handle_ins);
                insertion_to_vcf_record.emplace_back(i);
                vcf_record_to_insertion.at(i).emplace_back(handle_ins);
                if (deallocate_ins) record.alt.clear();
            }
            else if (record.sv_type==VcfReader::TYPE_BREAKEND && !record.is_breakend_single() && !record.is_breakend_virtual(chromosomes)) {
                // Virtual telomeric breakends and single breakends carry no information.
                p=record.get_info_field(VcfReader::MATEID_STR,0,tmp_buffer);
                if (p!=string::npos && bnd_ids.contains(tmp_buffer)) {
                    // NOP: breakend already processed in its mate entry
                }
                else {
                    bnd_ids.emplace(tmp_buffer);
                    orientation_cis=record.get_breakend_orientation_cis();
                    if (orientation_cis==1) first_positions.emplace_back(tmp_pair.first+1);
                    else if (orientation_cis==2) first_positions.emplace_back(tmp_pair.first);
                    trans_pos=record.get_breakend_pos();
                    trans_pos--; // Zero-based
                    orientation_trans=record.get_breakend_orientation_trans();
                    record.get_breakend_chromosome(tmp_buffer);
                    if (orientation_trans==1) chunk_first[tmp_buffer].emplace_back(trans_pos+1);
                    else if (orientation_trans==2) chunk_first[tmp_buffer].emplace_back(trans_pos);
                    record.get_breakend_inserted_sequence(tmp_buffer);
                    if (!tmp_buffer.empty()) {
                        const handle_t handle_ins = graph.create_handle(tmp_buffer);
                        insertion_handles.emplace_back(handle_ins);
                        insertion_to_vcf_record.emplace_back(i);
                        vcf_record_to_insertion.at(i).emplace_back(handle_ins);
                        if (deallocate_ins) record.alt.clear();
                    }
                }
            }
        }
    }
    for (auto& [key,value]: chunk_first) {
        if (value.size()<=1) continue;
        sort_and_compact_positions(value);
    }

    // Building all non-INS nodes and all reference edges
    for (auto& [key,value]: chunk_first) {
        first_positions = value;
        n_positions = first_positions.size();
        if (n_positions == 0) continue;
        vector<handle_t>& handles = node_handles.at(key);
        handles.clear();
        const string& chrom_sequence = chromosomes.at(key);
        p = first_positions.at(0);
        first_pos = max(p - flank_length, 0U);
        if (p != first_pos) handles.emplace(handles.begin(),graph.create_handle(chrom_sequence.substr(first_pos, p - first_pos)));
        for (j = 1; j < n_positions; j++) {
            p = first_positions.at(j - 1);
            handles.emplace_back(graph.create_handle(chrom_sequence.substr(p, first_positions.at(j) - p)));
        }
        p = first_positions.at(n_positions - 1);
        last_pos = min(p + flank_length, (uint32_t)chrom_sequence.length());
        if (p != last_pos) handles.emplace_back(graph.create_handle(chrom_sequence.substr(p, last_pos - p)));
        n_handles=handles.size();
        for (j=1; j<n_handles; j++) graph.create_edge(handles.at(j-1),handles.at(j));
    }

    // Building all non-reference edges
    const vector<handle_t>& handles = node_handles.at(CHROMOSOME_ID);
    bnd_ids.clear();
    for (i=0; i<n_vcf_records; i++) {
        VcfRecord& record = vcf_records.at(i);
        record.get_reference_coordinates(false,tmp_pair);
        if (tmp_pair.first==UINT32_MAX || tmp_pair.second==UINT32_MAX) continue;
        if (tmp_pair.first!=tmp_pair.second) {
            p=find_closest(CHROMOSOME_ID,tmp_pair.first);
            q=find_closest(CHROMOSOME_ID,tmp_pair.second);
            if (record.sv_type==VcfReader::TYPE_DELETION) add_edge(handles.at(p), handles.at(q+1), i);
            else if (record.sv_type==VcfReader::TYPE_DUPLICATION) add_edge(handles.at(q),handles.at(p+1),i);
            else if (record.sv_type==VcfReader::TYPE_INVERSION) {
                handle_t reversed_handle = graph.flip(handles.at(q));
                add_edge(handles.at(p), reversed_handle, i);
                reversed_handle = graph.flip(handles.at(p + 1));
                add_edge(reversed_handle, handles.at(q + 1), i);
            }
            else if (record.sv_type==VcfReader::TYPE_REPLACEMENT) {
                add_edge(handles.at(p), insertion_handles.at(i), i);
                add_edge(insertion_handles.at(i), handles.at(q+1), i);
            }
        }
        else {
            if (record.sv_type==VcfReader::TYPE_INSERTION && !record.is_symbolic) {
                p = find_closest(CHROMOSOME_ID, tmp_pair.first);
                add_edge(handles.at(p), insertion_handles.at(i), i);
                add_edge(insertion_handles.at(i), handles.at(p + 1), i);
            }
            else if (record.sv_type==VcfReader::TYPE_BREAKEND && !record.is_breakend_single() && !record.is_breakend_virtual(chromosomes)) {
                // Virtual telomeric breakends and single breakends carry no information.
                p=record.get_info_field(VcfReader::MATEID_STR,0,tmp_buffer);
                if (p!=string::npos && bnd_ids.contains(tmp_buffer)) {
                    // NOP: breakend already processed in its mate entry
                }
                else {
                    bnd_ids.emplace(tmp_buffer);
                    trans_pos=record.get_breakend_pos();
                    orientation_cis=record.get_breakend_orientation_cis();
                    orientation_trans=record.get_breakend_orientation_trans();
                    if (orientation_cis==0 || orientation_trans==0) {
                        // NOP: orientation could not be determined
                    }
                    else {
                        p=find_closest(CHROMOSOME_ID,tmp_pair.first);
                        handle_t handle_from = orientation_cis==1?handles.at(p):graph.flip(handles.at(p+1));
                        record.get_breakend_chromosome(tmp_buffer);
                        p=find_closest(tmp_buffer,trans_pos-1);
                        handle_t handle_to = orientation_trans==1?graph.flip(node_handles[tmp_buffer].at(p)):node_handles[tmp_buffer].at(p+1);
                        record.get_breakend_inserted_sequence(tmp_buffer);
                        if (tmp_buffer.empty()) add_edge(handle_from,handle_to,i);
                        else {
                            handle_t handle_ins = insertion_handles.at(i);
                            if (orientation_cis==2) graph.flip(handle_ins);
                            add_edge(handle_from,handle_ins,i);
                            add_edge(handle_ins,handle_to,i);
                        }
                    }
                }
            }
        }
    }
    graph.optimize(true);

    // Clearing temporary space
    chunk_first.clear(); node_handles.clear(); insertion_handles.clear(); bnd_ids.clear();
}


void VariantGraph::to_gfa(const path& gfa_path) const {
    const char GFA_SEPARATOR = '\t';
    const char OVERLAP_FIELD = '*';
    ofstream outstream;
    outstream.open(gfa_path);
    graph.for_each_handle([&](handle_t handle) { outstream << "S" << GFA_SEPARATOR << graph.get_id(handle) << GFA_SEPARATOR << graph.get_sequence(handle) << "\n"; });
    graph.for_each_edge([&](edge_t edge) { outstream << "L" << GFA_SEPARATOR << graph.get_id(edge.first) << GFA_SEPARATOR << (graph.get_is_reverse(edge.first)?'-':'+') << GFA_SEPARATOR << graph.get_id(edge.second) << GFA_SEPARATOR << (graph.get_is_reverse(edge.second)?'-':'+') << OVERLAP_FIELD << '\n'; });
    outstream.close();
}


void VariantGraph::mark_edge(const edge_t& query, const vector<edge_t>& edges, vector<uint8_t>& flags) const {
    const uint8_t N_EDGES = edges.size();
    const edge_t reverse(graph.flip(query.second),graph.flip(query.first));

    for (uint8_t i=0; i<N_EDGES; i++) {
        edge_t edge = edges.at(i);
        if (edge==query) {
            if (flags.at(i)==2) flags.at(i)=3;
            else flags.at(i)=1;
            return;
        }
        else if (reverse==query) {
            if (flags.at(i)==1) flags.at(i)=3;
            else flags.at(i)=2;
            return;
        }
    }
}


void VariantGraph::to_vcf(const path& vcf_path) {
    uint32_t i;
    uint8_t j, n_edges, n_forward, n_reverse, n_both;
    ofstream outstream;

    printed.reserve(n_vcf_records); flags.reserve(n_vcf_records);
    printed.clear();
    for (i=0; i<n_vcf_records; i++) {
        printed.emplace_back(false);
        n_edges=vcf_record_to_edge.at(i).size();
        if (i>=flags.size()) flags.emplace_back(vector<uint8_t>(n_edges));
        else { flags.at(i).clear(); flags.at(i).reserve(n_edges); }
        flags.at(i).clear();
        for (j=0; j<n_edges; j++) flags.at(i).emplace_back(0);
    }
    graph.for_each_path_handle([&](path_handle_t path) {
        for (i=0; i<n_vcf_records; i++) {
            if (printed.at(i)) continue;
            n_edges=vcf_record_to_edge.at(i).size();
            for (j=0; j<n_edges; j++) flags.at(i).at(j)=0;
        }
        step_handle_t from = graph.path_begin(path);  // `graph` breaks circular paths
        step_handle_t last = graph.path_back(path);  // `graph` breaks circular paths
        const step_handle_t& source = from;
        while (graph.has_next_step(from)) {
            const step_handle_t to = graph.get_next_step(from);
            if (to==last) break;
            edge_t edge = graph.edge_handle(graph.get_handle_of_step(from), graph.get_handle_of_step(to));
            for (const uint32_t& k: edge_to_vcf_record.at(edge)) mark_edge(edge, vcf_record_to_edge.at(k), flags.at(k));
            from=to;
        }
        for (i=0; i<n_vcf_records; i++) {
            if (printed.at(i)) continue;
            n_edges=vcf_record_to_edge.at(i).size();
            n_forward=0; n_reverse=0; n_both=0;
            for (j=0; j<n_edges; j++) {
                switch (flags.at(i).at(j)) {
                    case 1: n_forward++; break;
                    case 2: n_reverse++; break;
                    case 3: n_both++; break;
                }
            }
            if (n_forward+n_both==n_edges || n_reverse+n_both==n_edges) printed.at(i)=true;
        }
    });
    outstream.open(vcf_path);
    for (i=0; i<n_vcf_records; i++) {
        if (printed.at(i)) { vcf_records.at(i).print(outstream); outstream << '\n'; }
    }
    outstream.close();
}


void VariantGraph::to_vcf_paths(const path& vcf_path) {
    const string ID_PREFIX = "path";  // Arbitrary
    const uint8_t DEFAULT_QUAL = 60;  // Arbitrary
    char pos_base;
    uint32_t pos, id_generator;
    string old_sequence, new_sequence, alt_sequence;
    ofstream outstream;

    outstream.open(vcf_path); id_generator=0;
    graph.for_each_path_handle([&](path_handle_t path) {
        id_generator++;
        const step_handle_t source = graph.path_begin(path);  // `graph` breaks circular paths
        const step_handle_t destination = graph.path_back(path);  // `graph` breaks circular paths
        const handle_t source_handle = graph.get_handle_of_step(source);
        const handle_t destination_handle = graph.get_handle_of_step(destination);
        const pair<string,uint32_t> source_coordinate = node_to_chromosome.at(graph.get_id(source_handle));
        const pair<string,uint32_t> destination_coordinate = node_to_chromosome.at(graph.get_id(destination_handle));
        bool is_source_reversed = graph.get_is_reverse(source_handle);
        bool is_destination_reversed = graph.get_is_reverse(destination_handle);
        new_sequence.clear();
        if ( source_coordinate.first==destination_coordinate.first &&
             ( (!is_source_reversed && !is_destination_reversed && source_coordinate.second<destination_coordinate.second) ||
               (is_source_reversed && is_destination_reversed && source_coordinate.second>destination_coordinate.second)
               )
             ) {
            // Source and destination are on the same chromosome and in the same orientation. We represent this as a VCF
            // replacement.
            if (!is_source_reversed && !is_destination_reversed) {
                pos=source_coordinate.second+graph.get_length(source_handle);  // One-based
                pos_base=chromosomes.at(source_coordinate.first).at(pos-1);
                old_sequence=chromosomes.at(source_coordinate.first).substr(pos,destination_coordinate.second-pos);
                step_handle_t from = source;
                while (graph.has_next_step(from)) {
                    const step_handle_t to = graph.get_next_step(from);
                    if (to==destination) break;
                    // `get_handle_of_step()` guarantees that the sequence is already in the correct orientation.
                    new_sequence.append(graph.get_sequence(graph.get_handle_of_step(to)));
                    from=to;
                }
            }
            else {
                pos=destination_coordinate.second+graph.get_length(destination_handle);  // One-based
                pos_base=chromosomes.at(destination_coordinate.first).at(pos-1);
                old_sequence=chromosomes.at(destination_coordinate.first).substr(pos,source_coordinate.second-pos);
                step_handle_t from = destination;
                while (graph.has_previous_step(from)) {
                    const step_handle_t to = graph.get_previous_step(from);
                    if (to==source) break;
                    new_sequence.append(graph.get_sequence(graph.flip(graph.get_handle_of_step(to))));
                    from=to;
                }
            }
            outstream << source_coordinate.first << VcfReader::VCF_SEPARATOR << pos << VcfReader::VCF_SEPARATOR << ID_PREFIX << id_generator << VcfReader::VCF_SEPARATOR << pos_base << old_sequence << VcfReader::VCF_SEPARATOR << pos_base << new_sequence << VcfReader::VCF_SEPARATOR << DEFAULT_QUAL << VcfReader::VCF_SEPARATOR << VcfReader::PASS_STR << VcfReader::VCF_SEPARATOR << VcfReader::SVLEN_STR << VcfReader::INFO_ASSIGNMENT << (new_sequence.length()>old_sequence.length()?new_sequence.length()-old_sequence.length():old_sequence.length()-new_sequence.length()) << '\n';
        }
        else {
            // Any other configuration is represented as a single VCF breakend with inserted sequence.
            if (!is_source_reversed) {
                pos=source_coordinate.second+graph.get_length(source_handle);  // One-based
                pos_base=chromosomes.at(source_coordinate.first).at(pos-1);
                step_handle_t from = source;
                while (graph.has_next_step(from)) {
                    const step_handle_t to = graph.get_next_step(from);
                    if (to==destination) break;
                    // `get_handle_of_step()` guarantees that the sequence is already in the correct orientation.
                    new_sequence.append(graph.get_sequence(graph.get_handle_of_step(to)));
                    from=to;
                }
                alt_sequence.clear(); alt_sequence.append(pos_base+new_sequence);
                if (!is_destination_reversed) {
                    alt_sequence.append(VcfReader::BREAKEND_CHAR_OPEN+destination_coordinate.first);
                    alt_sequence.append(VcfReader::BREAKEND_SEPARATOR+to_string(destination_coordinate.second+1)+VcfReader::BREAKEND_CHAR_OPEN);
                }
                else {
                    alt_sequence.append(VcfReader::BREAKEND_CHAR_CLOSE+destination_coordinate.first);
                    alt_sequence.append(VcfReader::BREAKEND_SEPARATOR+to_string(destination_coordinate.second+graph.get_length(destination_handle))+VcfReader::BREAKEND_CHAR_CLOSE);
                }
                outstream << source_coordinate.first << VcfReader::VCF_SEPARATOR << pos << VcfReader::VCF_SEPARATOR << ID_PREFIX << id_generator << VcfReader::VCF_SEPARATOR << pos_base << VcfReader::VCF_SEPARATOR << alt_sequence << VcfReader::VCF_SEPARATOR << DEFAULT_QUAL << VcfReader::VCF_SEPARATOR << VcfReader::PASS_STR << VcfReader::VCF_SEPARATOR << VcfReader::SVTYPE_STR << VcfReader::INFO_ASSIGNMENT << VcfReader::BND_STR << "\n";
            }
            else {
                pos=source_coordinate.second+1;  // One-based
                pos_base=chromosomes.at(source_coordinate.first).at(pos-1);
                step_handle_t from = source;
                while (graph.has_next_step(from)) {
                    const step_handle_t to = graph.get_next_step(from);
                    if (to==destination) break;
                    // `get_handle_of_step()` guarantees that the sequence is already in the correct orientation.
                    new_sequence.append(graph.get_sequence(graph.get_handle_of_step(to)));
                    from=to;
                }
                reverse_complement(new_sequence);
                alt_sequence.clear();
                if (!is_destination_reversed) {
                    alt_sequence.append(VcfReader::BREAKEND_CHAR_OPEN+destination_coordinate.first);
                    alt_sequence.append(VcfReader::BREAKEND_SEPARATOR+to_string(destination_coordinate.second+1)+VcfReader::BREAKEND_CHAR_OPEN);
                }
                else {
                    alt_sequence.append(VcfReader::BREAKEND_CHAR_CLOSE+destination_coordinate.first);
                    alt_sequence.append(VcfReader::BREAKEND_SEPARATOR+to_string(destination_coordinate.second+graph.get_length(destination_handle))+VcfReader::BREAKEND_CHAR_CLOSE);
                }
                alt_sequence.append(new_sequence+pos_base);
                outstream << source_coordinate.first << VcfReader::VCF_SEPARATOR << pos << VcfReader::VCF_SEPARATOR << ID_PREFIX << id_generator << VcfReader::VCF_SEPARATOR << pos_base << VcfReader::VCF_SEPARATOR << alt_sequence << VcfReader::VCF_SEPARATOR << DEFAULT_QUAL << VcfReader::VCF_SEPARATOR << VcfReader::PASS_STR << VcfReader::VCF_SEPARATOR << VcfReader::SVTYPE_STR << VcfReader::INFO_ASSIGNMENT << VcfReader::BND_STR << "\n";
            }
        }
    });
    outstream.close();
}



}