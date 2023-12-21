#include "VariantGraph.hpp"
#include "Sequence.hpp"

using sv_merge::reverse_complement;

#include <cmath>
#include <ostream>

using std::min;
using std::max;
using std::find;
using std::distance;
using std::to_string;
using std::ofstream;


// use coord_t everywhere for chr coords from <misc.h>


namespace sv_merge {

VariantGraph::VariantGraph(const unordered_map<string,string>& chromosomes, const unordered_map<string,vector<interval_t>>& tandem_track):
        chromosomes(chromosomes),
        n_chromosomes(chromosomes.size()),
        tandem_track(tandem_track)
{ }


void VariantGraph::sort_and_compact_positions(vector<uint32_t>& positions) {
    sort(positions.begin(),positions.end());
    auto iterator = unique(positions.begin(),positions.end());
    positions.resize(distance(positions.begin(),iterator));
}


uint16_t VariantGraph::find_closest(const string& chromosome, uint32_t position) const {
    const vector<uint32_t>& positions = chunk_first.at(chromosome);
    auto iterator = lower_bound(positions.begin(),positions.end(),position);
    if (iterator==positions.end()) return positions.size()-1;
    uint16_t p = distance(positions.begin(),iterator);
    if (p==0 || positions.at(p)==position) return p;
    return position-positions.at(p-1)<=positions.at(p)-position?p-1:p;
}


void VariantGraph::add_nonreference_edge(const handle_t& from, const handle_t& to, uint32_t record_id) {
    if (graph.has_edge(from,to)) {
        const edge_t edge = graph.edge_handle(from,to);
        if (!edge_to_vcf_record.contains(edge)) edge_to_vcf_record[edge]={record_id};
        else edge_to_vcf_record.at(edge).emplace_back(record_id);
        vcf_record_to_edge.at(record_id).emplace_back(edge);
    }
    else {
        graph.create_edge(from,to);
        edge_t edge = graph.edge_handle(from,to);
        edge_to_vcf_record[edge]={record_id};
        vcf_record_to_edge.at(record_id).emplace_back(edge);
    }
}


bool VariantGraph::intersect(const interval_t& interval1, const interval_t& interval2) {
    const uint32_t from1 = interval1.first;
    const uint32_t to1 = interval1.second-1;
    const uint32_t from2 = interval2.first;
    const uint32_t to2 = interval2.second-1;
    return (from2>=from1 && from2<=to1) || (to2>=from1 && to2<=to1);
}


uint32_t VariantGraph::get_flank_boundary_right(const string& chromosome_id, uint32_t pos, uint16_t flank_length) {
    const uint32_t CHROMOSOME_LENGTH = (uint32_t)chromosomes.at(chromosome_id).length();
    const vector<interval_t>& intervals = tandem_track.at(chromosome_id);
    const uint32_t N_INTERVALS = intervals.size();

    if (N_INTERVALS==0) return min(pos+flank_length,CHROMOSOME_LENGTH-1);
    tmp_interval.first=pos+1;
    tmp_interval.second=pos+flank_length;
    auto iter = std::lower_bound(intervals.begin(),intervals.end(),tmp_interval);
    uint32_t p = iter-intervals.begin();
    if (p>0 && intersect(tmp_interval,intervals.at(p-1))) {
        tmp_interval.first=intervals.at(p-1).second;
        tmp_interval.second=tmp_interval.first+flank_length;
    }
    while (p<N_INTERVALS && intersect(tmp_interval,intervals.at(p))) {
        tmp_interval.first=intervals.at(p).second;
        tmp_interval.second=tmp_interval.first+flank_length;
        p++;
    }
    return tmp_interval.first>=CHROMOSOME_LENGTH?UINT32_MAX:min((uint32_t)tmp_interval.second-1,CHROMOSOME_LENGTH-1);
}


uint32_t VariantGraph::get_flank_boundary_left(const string& chromosome_id, uint32_t pos, uint16_t flank_length) {
    const vector<interval_t>& intervals = tandem_track.at(chromosome_id);
    const uint32_t N_INTERVALS = intervals.size();

    if (N_INTERVALS==0) return pos>=flank_length?pos-flank_length:0;
    tmp_interval.first=pos>=flank_length?pos-flank_length:0;
    tmp_interval.second=pos>0?pos-1:0;
    auto iter = std::lower_bound(intervals.begin(),intervals.end(),tmp_interval);
    uint32_t p = iter-intervals.begin();
    if (p<N_INTERVALS && intersect(tmp_interval,intervals.at(p))) {
        tmp_interval.second=intervals.at(p).first;
        tmp_interval.first=tmp_interval.second>=flank_length?tmp_interval.second-flank_length:0;
    }
    if (p>0) {
        p--;
        while (intersect(tmp_interval,intervals.at(p))) {
            tmp_interval.second=intervals.at(p).first;
            tmp_interval.first=tmp_interval.second>=flank_length?tmp_interval.second-flank_length:0;
            if (p==0) break;
            else p--;
        }
    }
    return tmp_interval.second==0?UINT32_MAX:(uint32_t)tmp_interval.first;
}


void VariantGraph::build(vector<VcfRecord>& records, uint16_t flank_length, bool deallocate_ins) {
    graph.clear();
    this->vcf_records=std::move(records);
    n_vcf_records=vcf_records.size();
    bnd_ids.clear(); bnd_ids.reserve(n_vcf_records);
    chunk_first.clear(); chunk_first.reserve(n_chromosomes);
    node_handles.clear(); node_handles.reserve(n_chromosomes);
    node_to_chromosome.clear(); node_to_chromosome.reserve(n_vcf_records);
    insertion_handles.clear(); insertion_handles.reserve(n_vcf_records);
    edge_to_vcf_record.clear(); edge_to_vcf_record.reserve(n_vcf_records);
    if (vcf_records.empty()) {
        vcf_record_to_edge.clear();
        return;
    }

    const string CHROMOSOME_ID = vcf_records.at(0).chrom;
    const uint8_t BIN_SIZE = 50;  // Arbitrary
    uint8_t orientation_cis, orientation_trans;
    uint16_t n_positions, n_handles;
    uint32_t i, j, p, q, first_pos, last_pos, trans_pos;
    size_t r;
    string tmp_buffer;
    pair<uint32_t, uint32_t> tmp_pair;

    // Collecting every distinct first position of a reference node, and building all non-reference nodes.
    cerr << "VCF records: " << n_vcf_records << '\n';
    cerr << "Histogram: SV type | SV length (binned by " << BIN_SIZE << "bp) | Number of occurrences\n";
    print_sv_lengths(BIN_SIZE);  // Arbitrary
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
                if (deallocate_ins) { record.ref.clear(); record.alt.clear(); }
            }
        }
        else {
            if (record.sv_type==VcfReader::TYPE_INSERTION && !record.is_symbolic) {
                first_positions.emplace_back(tmp_pair.first);
                const handle_t handle_ins = graph.create_handle(record.alt.substr(1));
                insertion_handles.emplace_back(handle_ins);
                if (deallocate_ins) record.alt.clear();
            }
            else if (record.sv_type==VcfReader::TYPE_BREAKEND && !record.is_breakend_single() && !record.is_breakend_virtual(chromosomes)) {
                // Virtual telomeric breakends and single breakends carry no information.
                r=record.get_info_field(VcfReader::MATEID_STR,0,tmp_buffer);
                if (r!=string::npos && bnd_ids.contains(tmp_buffer)) {
                    // NOP: the breakend was already processed at its mate record.
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
                        if (deallocate_ins) record.alt.clear();
                    }
                }
            }
        }
    }
    cerr << "Table: Chromosome involved | N. distinct breakpoints\n";
    for (auto& [key,value]: chunk_first) {
        uint32_t n_breakpoints = value.size();
        if (n_breakpoints==0) continue;
        if (n_breakpoints>1) sort_and_compact_positions(value);
        cerr << key << ',' << value.size() << '\n';
    }
    cerr << "Number of non-reference nodes: " << insertion_handles.size() << '\n';

    // Building all reference nodes and all reference edges
    for (auto& [key,value]: chunk_first) {
        first_positions=value;
        n_positions=first_positions.size();
        if (n_positions==0) continue;
        vector<handle_t>& handles = node_handles.at(key);
        handles.clear();
        const string& chrom_sequence = chromosomes.at(key);
        p=first_positions.at(0);
        first_pos=get_flank_boundary_left(key,p,flank_length);
        if (first_pos==UINT32_MAX) first_pos=0;
        if (p!=first_pos) {
            handle_t reference_handle = graph.create_handle(chrom_sequence.substr(first_pos,p-first_pos));
            handles.emplace_back(reference_handle);
            node_to_chromosome[graph.get_id(reference_handle)]=pair<string,uint32_t>(key,p);
        }
        for (j=1; j<n_positions; j++) {
            p=first_positions.at(j-1);
            handle_t reference_handle = graph.create_handle(chrom_sequence.substr(p,first_positions.at(j)-p));
            handles.emplace_back(reference_handle);
            node_to_chromosome[graph.get_id(reference_handle)]=pair<string,uint32_t>(key,p);
        }
        p=first_positions.at(n_positions-1);
        last_pos=get_flank_boundary_right(key,p,flank_length);
        if (last_pos==UINT32_MAX) last_pos=chrom_sequence.length()-1;
        if (p!=last_pos) {
            handle_t reference_handle = graph.create_handle(chrom_sequence.substr(p,last_pos+1-p));
            handles.emplace_back(reference_handle);
            node_to_chromosome[graph.get_id(reference_handle)]=pair<string,uint32_t>(key,p);
        }
        n_handles=handles.size();
        for (j=1; j<n_handles; j++) graph.create_edge(handles.at(j-1),handles.at(j));
    }
    cerr << "Table: Chromosome involved | N. reference nodes\n";
    for (auto& [key,value]: node_handles) {
        const uint32_t n_chunks = value.size();
        if (n_chunks>0) cerr << key << ',' << n_chunks << '\n';
    }
    const uint32_t n_reference_edges = graph.get_edge_count();
    cerr << "Number of reference edges: " << n_reference_edges << '\n';

    // Building all non-reference edges
    vcf_record_to_edge.reserve(n_vcf_records);
    for (i=0; i<vcf_record_to_edge.size(); i++) vcf_record_to_edge.at(i).clear();
    for (i=vcf_record_to_edge.size(); i<n_vcf_records; i++) vcf_record_to_edge.emplace_back();
    const vector<handle_t>& handles = node_handles.at(CHROMOSOME_ID);
    bnd_ids.clear();
    for (i=0; i<n_vcf_records; i++) {
        VcfRecord& record = vcf_records.at(i);
        record.get_reference_coordinates(false,tmp_pair);
        if (tmp_pair.first==UINT32_MAX || tmp_pair.second==UINT32_MAX) continue;
        if (tmp_pair.first!=tmp_pair.second) {
            p=find_closest(CHROMOSOME_ID,tmp_pair.first);
            q=find_closest(CHROMOSOME_ID,tmp_pair.second);
            if (record.sv_type==VcfReader::TYPE_DELETION) add_nonreference_edge(handles.at(p), handles.at(q+1), i);
            else if (record.sv_type==VcfReader::TYPE_DUPLICATION || record.sv_type==VcfReader::TYPE_CNV) add_nonreference_edge(handles.at(q),handles.at(p+1),i);
            else if (record.sv_type==VcfReader::TYPE_INVERSION) {
                handle_t reversed_handle = graph.flip(handles.at(q));
                add_nonreference_edge(handles.at(p), reversed_handle, i);
                reversed_handle=graph.flip(handles.at(p+1));
                add_nonreference_edge(reversed_handle,handles.at(q+1),i);
            }
            else if (record.sv_type==VcfReader::TYPE_REPLACEMENT) {
                add_nonreference_edge(handles.at(p),insertion_handles.at(i),i);
                add_nonreference_edge(insertion_handles.at(i),handles.at(q+1),i);
            }
        }
        else {
            if (record.sv_type==VcfReader::TYPE_INSERTION && !record.is_symbolic) {
                p=find_closest(CHROMOSOME_ID,tmp_pair.first);
                handle_t& insertion_handle = insertion_handles.at(i);
                add_nonreference_edge(handles.at(p),insertion_handle,i);
                add_nonreference_edge(insertion_handle,handles.at(p+1),i);
            }
            else if (record.sv_type==VcfReader::TYPE_BREAKEND && !record.is_breakend_single() && !record.is_breakend_virtual(chromosomes)) {
                // Virtual telomeric breakends and single breakends carry no information.
                r=record.get_info_field(VcfReader::MATEID_STR,0,tmp_buffer);
                if (r!=string::npos && bnd_ids.contains(tmp_buffer)) {
                    // NOP: the breakend was already processed at its mate record.
                }
                else {
                    bnd_ids.emplace(tmp_buffer);
                    trans_pos=record.get_breakend_pos();
                    orientation_cis=record.get_breakend_orientation_cis();
                    orientation_trans=record.get_breakend_orientation_trans();
                    if (orientation_cis==0 || orientation_trans==0) {
                        // NOP: orientation could not be determined.
                    }
                    else {
                        p=find_closest(CHROMOSOME_ID,tmp_pair.first);
                        handle_t handle_from = orientation_cis==1?handles.at(p):graph.flip(handles.at(p+1));
                        record.get_breakend_chromosome(tmp_buffer);
                        p=find_closest(tmp_buffer,trans_pos-1);
                        handle_t handle_to = orientation_trans==1?graph.flip(node_handles.at(tmp_buffer).at(p)):node_handles.at(tmp_buffer).at(p+1);
                        record.get_breakend_inserted_sequence(tmp_buffer);
                        if (tmp_buffer.empty()) add_nonreference_edge(handle_from,handle_to,i);
                        else {
                            handle_t handle_ins = insertion_handles.at(i);
                            if (orientation_cis==2) graph.flip(handle_ins);
                            add_nonreference_edge(handle_from,handle_ins,i);
                            add_nonreference_edge(handle_ins,handle_to,i);
                        }
                    }
                }
            }
        }
    }
    cerr << "Number of non-reference edges: " << (graph.get_edge_count()-n_reference_edges) << '\n';
    cerr << "Number of inter-chromosomal edges: " << get_n_interchromosomal_edges() << '\n';
    cerr << "Number of nodes with edges on just one side: " << get_n_dangling_nodes() << '\n';
    print_edge_histograms();

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

    printed.reserve(n_vcf_records); printed.clear();
    flags.reserve(n_vcf_records);
    for (i=flags.size(); i<n_vcf_records; i++) flags.emplace_back();
    for (i=0; i<n_vcf_records; i++) {
        printed.emplace_back(false);
        flags.at(i).clear(); flags.at(i).reserve(vcf_record_to_edge.at(i).size());
    }
    graph.for_each_path_handle([&](path_handle_t path) {
        for (i=0; i<n_vcf_records; i++) {
            if (printed.at(i)) continue;
            n_edges=vcf_record_to_edge.at(i).size();
            for (j=0; j<n_edges; j++) flags.at(i).emplace_back(0);
        }
        step_handle_t from = graph.path_begin(path);  // `graph` breaks circular paths
        const step_handle_t last = graph.path_back(path);  // `graph` breaks circular paths
        while (graph.has_next_step(from)) {
            const step_handle_t to = graph.get_next_step(from);
            if (to==last) break;
            edge_t edge = graph.edge_handle(graph.get_handle_of_step(from),graph.get_handle_of_step(to));
            for (const uint32_t& k: edge_to_vcf_record.at(edge)) mark_edge(edge,vcf_record_to_edge.at(k),flags.at(k));
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
        const nid_t id1 = graph.get_id(source_handle);
        const nid_t id2 = graph.get_id(destination_handle);
        if (!node_to_chromosome.contains(id1) || !node_to_chromosome.contains(id2)) return;
        const pair<string,uint32_t> source_coordinate = node_to_chromosome.at(id1);
        const pair<string,uint32_t> destination_coordinate = node_to_chromosome.at(id2);
        const bool is_source_reversed = graph.get_is_reverse(source_handle);
        const bool is_destination_reversed = graph.get_is_reverse(destination_handle);
        new_sequence.clear();
        if ( source_coordinate.first==destination_coordinate.first &&
             ( (!is_source_reversed && !is_destination_reversed && source_coordinate.second<destination_coordinate.second) ||
               (is_source_reversed && is_destination_reversed && source_coordinate.second>destination_coordinate.second)
               )
             ) {
            // Source and destination are on the same chromosome and in the same orientation. We represent this as a VCF
            // replacement.
            if (!is_source_reversed) {
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
            outstream << source_coordinate.first << VcfReader::VCF_SEPARATOR << pos << VcfReader::VCF_SEPARATOR << ID_PREFIX << id_generator << VcfReader::VCF_SEPARATOR << pos_base << old_sequence << VcfReader::VCF_SEPARATOR << pos_base << new_sequence << VcfReader::VCF_SEPARATOR << DEFAULT_QUAL << VcfReader::VCF_SEPARATOR << VcfReader::PASS_STR << VcfReader::VCF_SEPARATOR << VcfReader::VCF_MISSING_CHAR << '\n';
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


void VariantGraph::get_ids_of_dangling_nodes(unordered_set<nid_t>& out) const {
    out.clear();
    graph.for_each_handle([&](handle_t node) {
        if (graph.get_degree(node,true)==0 || graph.get_degree(node,false)==0) out.emplace(graph.get_id(node));
    });
}


uint32_t VariantGraph::get_n_dangling_nodes() const {
    uint32_t out = 0;
    graph.for_each_handle([&](handle_t node) {
        if (graph.get_degree(node,true)==0 || graph.get_degree(node,false)==0) out++;
    });
    return out;
}


uint32_t VariantGraph::get_n_interchromosomal_edges() const {
    uint32_t out = 0;
    graph.for_each_edge([&](edge_t edge) {
        if (node_to_chromosome.at(graph.get_id(edge.first)).first!=node_to_chromosome.at(graph.get_id(edge.second)).first) out++;
    });
    return out;
}


void VariantGraph::print_edge_histograms() const {
    uint32_t n_edges = edge_to_vcf_record.size();
    if (n_edges==0) return;

    uint32_t i, last;
    vector<uint32_t> histogram;

    for (i=0; i<=n_vcf_records; i++) histogram.emplace_back(0);
    for (auto& [key,value]: edge_to_vcf_record) histogram.at(value.size())++;
    cerr << "Histogram: X = n. VCF records | Y = n. non-reference edges supported by X VCF records\n";
    for (last=n_vcf_records; ; last--) { if (histogram.at(i)>0) break; }
    for (i=0; i<=last; i++) cerr << i << ',' << histogram.at(i) << '\n';

    histogram.reserve(n_edges+1);
    for (i=0; i<histogram.size(); i++) histogram.at(i)=0;
    for (i=histogram.size(); i<=n_edges; i++) histogram.emplace_back(0);
    for (i=0; i<n_vcf_records; i++) histogram.at(vcf_record_to_edge.at(i).size())++;
    cerr << "Histogram: X = n. non-reference edges | Y = n. VCF records creating X non-reference edges\n";
    for (last=n_edges; ; last--) { if (histogram.at(i)>0) break; }
    for (i=0; i<=last; i++) cerr << i << ',' << histogram.at(i) << '\n';

    n_edges=graph.get_edge_count();
    histogram.reserve(n_edges+1);
    for (i=0; i<histogram.size(); i++) histogram.at(i)=0;
    for (i=histogram.size(); i<=n_edges; i++) histogram.emplace_back(0);
    graph.for_each_handle([&](handle_t node) { histogram.at(max(graph.get_degree(node,false),graph.get_degree(node,true)))++; });
    cerr << "Histogram: X = max degree of a node | Y = n. nodes with max degree X\n";
    for (last=n_edges; ; last--) { if (histogram.at(i)>0) break; }
    for (i=0; i<=last; i++) cerr << i << ',' << histogram.at(i) << '\n';
}


void VariantGraph::print_vcf_records_stats(uint8_t bin_size) const {
    vector<uint32_t> lengths_ins, lengths_del, lengths_inv;
    vector<pair<string,string>> bnds;
    string buffer;

    for (auto& record: vcf_records) {
        if (record.sv_type==VcfReader::TYPE_INSERTION || record.sv_type==VcfReader::TYPE_DUPLICATION) lengths_ins.emplace_back(record.sv_length/bin_size);
        else if (record.sv_type==VcfReader::TYPE_DELETION) lengths_del.emplace_back(record.sv_length/bin_size);
        else if (record.sv_type==VcfReader::TYPE_INVERSION) lengths_inv.emplace_back(record.sv_length/bin_size);
        else if (record.sv_type==VcfReader::TYPE_BREAKEND) {
            record.get_breakend_chromosome(buffer);
            if (record.chrom<buffer) bnds.emplace_back(record.chrom,buffer);
            else bnds.emplace_back(buffer,record.chrom);
        }
    }
    print_sv_lengths(lengths_ins,bin_size,VcfReader::INS_STR);
    print_sv_lengths(lengths_del,bin_size,VcfReader::DEL_STR);
    print_sv_lengths(lengths_inv,bin_size,VcfReader::INV_STR);
}


void VariantGraph::print_sv_lengths(vector<uint32_t>& lengths, uint8_t bin_size, const string& prefix) const {
    uint32_t i, previous, previous_count;
    const uint32_t SIZE = lengths.size();

    if (lengths.size()>1) sort(lengths.begin(),lengths.end());
    previous=lengths.at(0); previous_count=1;
    for (i=1; i<SIZE; i++) {
        if (lengths.at(i)!=previous) {
            cerr << prefix << ',' << (previous*bin_size) << ',' << previous_count << '\n';
            previous=lengths.at(i); previous_count=1;
        }
    }
    cerr << prefix << ',' << (previous*bin_size) << ',' << previous_count << '\n';
}


-------------->
void VariantGraph::print_bnds(vector<pair<string,string>>& bnds) const {
    uint32_t i, previous, previous_count;
    const uint32_t SIZE = lengths.size();

    if (lengths.size()>1) sort(lengths.begin(),lengths.end());
    previous=lengths.at(0); previous_count=1;
    for (i=1; i<SIZE; i++) {
        if (lengths.at(i)!=previous) {
            cerr << prefix << ',' << (previous*bin_size) << ',' << previous_count << '\n';
            previous=lengths.at(i); previous_count=1;
        }
    }
    cerr << prefix << ',' << (previous*bin_size) << ',' << previous_count << '\n';
}


}