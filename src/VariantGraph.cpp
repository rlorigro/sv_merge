#include "VariantGraph.hpp"
#include "Sequence.hpp"

using sv_merge::reverse_complement;
using bdsg::step_handle_t;
using bdsg::path_handle_t;

#include <cmath>
#include <cstring>

using std::min;
using std::max;
using std::sort;
using std::binary_search;
using std::move;
using std::make_move_iterator;
using std::unique;
using std::find;
using std::distance;
using std::to_string;
using std::streamsize;
using std::toupper;
using std::set;


namespace sv_merge {

const char VariantGraph::GFA_SEPARATOR = '\t';
const char VariantGraph::GFA_NODE_CHAR = 'S';
const char VariantGraph::GFA_LINK_CHAR = 'L';
const char VariantGraph::GFA_PATH_CHAR = 'P';
const string VariantGraph::GFA_OVERLAP_FIELD = "0M";
const char VariantGraph::GFA_PATH_SEPARATOR = ',';
const char VariantGraph::GFA_FWD_CHAR = '+';
const char VariantGraph::GFA_REV_CHAR = '-';
const char VariantGraph::GFA_LINE_END = '\n';
const char VariantGraph::GAF_FWD_CHAR = '>';
const char VariantGraph::GAF_REV_CHAR = '<';
const uint64_t VariantGraph::STREAMSIZE_MAX = numeric_limits<streamsize>::max();


VariantGraph::VariantGraph(const unordered_map<string,string>& chromosomes, const unordered_map<string,vector<interval_t>>& tandem_track, int32_t min_variant_length, bool force_uppercase, bool silent):
        chromosomes(chromosomes),
        n_chromosomes(chromosomes.size()),
        tandem_track(tandem_track),
        min_variant_length(min_variant_length),
        force_uppercase(force_uppercase),
        silent(silent)
{ }


void VariantGraph::sort_and_compact_positions(vector<int32_t>& positions) {
    if (positions.size()<=1) return;
    sort(positions.begin(),positions.end());
    const auto iterator = unique(positions.begin(),positions.end());
    positions.resize(distance(positions.begin(),iterator));
}


size_t VariantGraph::find_closest(const string& chromosome, int32_t position) const {
    const vector<int32_t>& positions = chunk_first.at(chromosome);
    const auto iterator = lower_bound(positions.begin(),positions.end(),position);
    if (iterator==positions.end()) return positions.size()-1;
    const auto p = distance(positions.begin(),iterator);
    if (p==0 || positions.at(p)==position) return p;
    return position-positions.at(p-1)<=positions.at(p)-position?p-1:p;
}


void VariantGraph::add_nonreference_edge(const handle_t& from, const handle_t& to, size_t record_id) {
    if (graph.has_edge(from,to)) {
        const edge_t edge = graph.edge_handle(from,to);
        if (record_id!=SIZE_MAX) {
            if (!edge_to_vcf_record.contains(edge)) edge_to_vcf_record[edge]={record_id};
            else edge_to_vcf_record.at(edge).emplace_back(record_id);
            vcf_record_to_edge.at(record_id).emplace_back(edge);
        }
    }
    else {
        graph.create_edge(from,to);
        if (record_id!=SIZE_MAX) {
            edge_t edge = graph.edge_handle(from,to);
            edge_to_vcf_record[edge]={record_id};
            vcf_record_to_edge.at(record_id).emplace_back(edge);
        }
    }
}


bool VariantGraph::intersect(const interval_t& interval1, const interval_t& interval2) {
    const size_t from1 = interval1.first;
    const size_t to1 = interval1.second-1;
    const size_t from2 = interval2.first;
    const size_t to2 = interval2.second-1;
    return max(from1,from2)<=min(to1,to2);
}


int32_t VariantGraph::get_flank_boundary_right(const string& chromosome_id, int32_t pos, int32_t flank_length) {
    const auto CHROMOSOME_LENGTH = (int32_t)chromosomes.at(chromosome_id).length();
    if (pos==CHROMOSOME_LENGTH-1) return INT32_MAX;
    if (!tandem_track.contains(chromosome_id)) return min(pos+flank_length-1,CHROMOSOME_LENGTH-1);
    const vector<interval_t>& intervals = tandem_track.at(chromosome_id);
    return get_flank_boundary_right_impl(pos,CHROMOSOME_LENGTH,intervals,flank_length);
}


int32_t VariantGraph::get_flank_boundary_right_impl(int32_t pos, int32_t sequence_length, const vector<interval_t>& intervals, int32_t flank_length) {
    const size_t N_INTERVALS = intervals.size();
    if (N_INTERVALS==0) return min(pos+flank_length-1,sequence_length-1);
    tmp_interval.first=pos; tmp_interval.second=pos+flank_length;
    auto iter = std::lower_bound(intervals.begin(),intervals.end(),tmp_interval);
    size_t p = iter-intervals.begin();
    if (p>0 && intersect(tmp_interval,intervals.at(p-1))) {
        tmp_interval.first=intervals.at(p-1).second;
        tmp_interval.second=tmp_interval.first+flank_length;
    }
    while (p<N_INTERVALS && intersect(tmp_interval,intervals.at(p))) {
        tmp_interval.first=intervals.at(p).second;
        tmp_interval.second=tmp_interval.first+flank_length;
        p++;
    }
    return tmp_interval.first>=sequence_length?INT32_MAX:min(tmp_interval.second-1,sequence_length-1);
}


int32_t VariantGraph::get_flank_boundary_left(const string& chromosome_id, int32_t pos, int32_t flank_length) {
    if (pos==0) return INT32_MAX;
    if (!tandem_track.contains(chromosome_id)) return pos>=flank_length?pos-flank_length:0;
    const vector<interval_t>& intervals = tandem_track.at(chromosome_id);
    return get_flank_boundary_left_impl(pos,intervals,flank_length);
}


int32_t VariantGraph::get_flank_boundary_left_impl(int32_t pos, const vector<interval_t>& intervals, int32_t flank_length) {
    const size_t N_INTERVALS = intervals.size();
    if (N_INTERVALS==0) return pos>=flank_length?pos-flank_length:0;
    tmp_interval.first=pos>=flank_length?pos-flank_length:0; tmp_interval.second=pos;
    auto iter = std::lower_bound(intervals.begin(),intervals.end(),tmp_interval);
    size_t p = iter-intervals.begin();
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
    return tmp_interval.second==0?INT32_MAX:tmp_interval.first;
}


void VariantGraph::get_tandem_intervals(nid_t node_id, int32_t node_length, bool is_reverse, int32_t offset, vector<interval_t>& out) {
    size_t i, from;

    const auto OUT_FIRST = out.size();
    const pair<string,int32_t> COORDINATE = node_to_chromosome.find(node_id)->second;
    const string CHROMOSOME_ID = COORDINATE.first;
    if (!tandem_track.contains(CHROMOSOME_ID)) return;
    const vector<interval_t>& INTERVALS = tandem_track.at(CHROMOSOME_ID);
    if (INTERVALS.empty()) return;
    const size_t N_INTERVALS = INTERVALS.size();
    tmp_interval.first=COORDINATE.second;
    tmp_interval.second=tmp_interval.first+node_length;
    auto ITER = std::lower_bound(INTERVALS.begin(),INTERVALS.end(),tmp_interval);
    if (ITER!=INTERVALS.end()) i=ITER-INTERVALS.begin()+1;
    else i=N_INTERVALS;
    from=N_INTERVALS;
    while (i>0) {
        i--;
        if (INTERVALS.at(i).second<=tmp_interval.first) { from=i+1; break; }
    }
    i=from;
    while (i<N_INTERVALS) {
        if (INTERVALS.at(i).first>=tmp_interval.second) break;
        out.emplace_back(max(tmp_interval.first,INTERVALS.at(i).first),min(tmp_interval.second,INTERVALS.at(i).second));
        i++;
    }
    if (!is_reverse) {
        for (i=OUT_FIRST; i<out.size(); i++) {
            out.at(i).first=out.at(i).first-tmp_interval.first+offset;
            out.at(i).second=out.at(i).second-tmp_interval.first+offset;
        }
    }
    else {
        std::reverse(out.begin()+OUT_FIRST,out.end());
        for (i=OUT_FIRST; i<out.size(); i++) {
            std::swap(out.at(i).first,out.at(i).second);
            out.at(i).first=(tmp_interval.second-1)-(out.at(i).first-1)+offset;
            out.at(i).second=(tmp_interval.second-1)-(out.at(i).second-1)+offset;
        }
    }
}


void VariantGraph::get_tandem_intervals(const vector<nid_t>& node_ids, const vector<bool>& is_reverse, const vector<int32_t>& node_lengths, vector<interval_t>& out) {
    const size_t LENGTH = node_ids.size();
    size_t i, j;
    int32_t offset, length, length_prime;
    nid_t node_id;

    // Collecting per-node tandem intervals, including entire non-ref nodes.
    out.clear(); offset=0; i=0;
    while (i<LENGTH) {
        node_id=node_ids.at(i);
        length=node_lengths.at(i);
        if (is_reference_node(node_id)) {
            get_tandem_intervals(node_id,length,is_reverse.at(i),offset,out);
            offset+=length;
            i++;
        }
        else {
            if (!out.empty() && out.at(out.size()-1).second==offset) {
                out.at(out.size()-1).second=offset+length;
                offset+=length;
                i++;
            }
            else if (i<LENGTH-1) {
                j=out.size();
                i++;
                node_id=node_ids.at(i);
                length_prime=node_lengths.at(i);
                if (!is_reference_node(node_id)) throw runtime_error("ERROR: the path contains a non-ref node adjacent to another non-ref node?!");
                get_tandem_intervals(node_id,length_prime,is_reverse.at(i),offset+length,out);
                if (out.size()>j && out.at(j).first==offset+length) out.at(j).first=offset;
                offset+=length+length_prime;
                i++;
            }
            else {
                offset+=length;
                i++;
            }
        }
    }
    if (out.empty()) return;

    // Compacting adjacent intervals
    j=0;
    for (i=1; i<out.size(); i++) {
        if (out.at(i).first==out.at(j).second) out.at(j).second=out.at(i).second;
        else std::swap(out.at(++j),out.at(i));
    }
    out.erase(out.begin()+(j+1),out.end());
}


/**
 * Remark: this is implemented naively in `O(|path|^2)` time.
 *
 * Remark: the procedure uses global buffers `node_ids_buffer, is_reverse_buffer, node_lengths_buffer,
 * intervals_buffer`.
 */
void VariantGraph::vcf_record_to_path_intervals(const vector<pair<string,bool>>& path, const vector<edge_t>& edges_of_the_record, int32_t flank_length, vector<interval_t>& out) {
    const auto PATH_LENGTH = (int32_t)path.size();
    if (PATH_LENGTH<2) { out.clear(); return; }
    const auto N_EDGES = (int32_t)edges_of_the_record.size();
    bool found, rev;
    int32_t i, j, k, p, x, y;
    int32_t length, path_length_bps, first;
    nid_t id;
    handle_t handle1, handle2;

    node_ids_buffer.clear(); is_reverse_buffer.clear(); node_lengths_buffer.clear(); tmp_edges.clear();
    id=stoll(path.at(0).first);
    rev=path.at(0).second;
    node_ids_buffer.emplace_back(id);
    is_reverse_buffer.emplace_back(rev);
    handle1=graph.get_handle(id,rev);
    length=int32_t(graph.get_length(handle1));
    node_lengths_buffer.emplace_back(length);
    path_length_bps=length;
    for (i=1; i<PATH_LENGTH; i++) {
        id=stoll(path.at(i).first);
        rev=path.at(i).second;
        node_ids_buffer.emplace_back(id);
        is_reverse_buffer.emplace_back(rev);
        handle2=graph.get_handle(id,rev);
        length=int32_t(graph.get_length(handle2));
        node_lengths_buffer.emplace_back(length);
        path_length_bps+=length;
        tmp_edges.emplace_back(graph.edge_handle(handle1,handle2));  // Canonized
        handle1=handle2;
    }
    get_tandem_intervals(node_ids_buffer,is_reverse_buffer,node_lengths_buffer,intervals_buffer);
    out.clear(); x=0;
    for (i=0; i<PATH_LENGTH-1; i++) {
        x+=node_lengths_buffer.at(i);
        first=0;
        for (k=0; k<N_EDGES; k++) {
            if (edges_of_the_record.at(k)!=null_edge) continue;
            found=false;
            if (tmp_edges.at(i)==edges_of_the_record.at(first)) {
                // Forward match
                j=i; p=first+1; y=x;
                while (j<(PATH_LENGTH-1)-1 && p<k) {
                    j++;
                    y+=node_lengths_buffer.at(j);
                    if (tmp_edges.at(j)==edges_of_the_record.at(p)) p++;
                }
                if (p==k) found=true;
            }
            if (!found && tmp_edges.at(i)==edges_of_the_record.at(k-1)) {
                // Reverse match
                j=i; p=k-2; y=x;
                while (j<(PATH_LENGTH-1)-1 && p>=first) {
                    j++;
                    y+=node_lengths_buffer.at(j);
                    if (tmp_edges.at(j)==edges_of_the_record.at(p)) p--;
                }
                if (p==first-1) found=true;
            }
            if (found) {
                int32_t left = get_flank_boundary_left_impl(x,intervals_buffer,flank_length);
                if (left==INT32_MAX) left=0;
                int32_t right = get_flank_boundary_right_impl(y,path_length_bps,intervals_buffer,flank_length);
                if (right==INT32_MAX) right=path_length_bps-1;
                out.emplace_back(left,right);
                break;
            }
            first=k+1;
        }
    }
}


handle_t VariantGraph::sequence2handle(const string& sequence, string& tmp_buffer) {
    if (force_uppercase) {
        tmp_buffer.clear();
        tmp_buffer.append(sequence);
        const size_t length = sequence.length();
        for (size_t i=0; i<length; i++) tmp_buffer.at(i)=(char)toupper(tmp_buffer.at(i));
        return graph.create_handle(tmp_buffer);
    }
    else return graph.create_handle(sequence);
}


void VariantGraph::build(vector<VcfRecord>& records, int32_t flank_length, int32_t interior_flank_length, int32_t x, int32_t y, bool deallocate_ref_alt, const vector<string>& callers, bool acyclic, bool graph_closure) {
    graph.clear();
    n_vcf_records=records.size();
    if (n_vcf_records!=0) this->vcf_records=std::move(records);
    else this->vcf_records.clear();
    bnd_ids.clear(); if (n_vcf_records!=0) bnd_ids.reserve(n_vcf_records);
    node_handles.clear(); if (n_vcf_records!=0) node_handles.reserve(n_chromosomes);
    for (const auto& chromosome: chromosomes) node_handles.emplace(chromosome.first,vector<handle_t>());
    node_to_chromosome.clear(); if (n_vcf_records!=0) node_to_chromosome.reserve(n_vcf_records);
    insertion_handles.clear(); if (n_vcf_records!=0) insertion_handles.reserve(n_vcf_records);
    edge_to_vcf_record.clear(); if (n_vcf_records!=0) edge_to_vcf_record.reserve(n_vcf_records);
    if (n_vcf_records==0) { vcf_record_to_edge.clear(); return; }
    if (acyclic) interval_to_insertion_handle.clear();
    null_edge.first={0}; null_edge.second={0};

    main_chromosome=vcf_records.at(0).chrom;
    main_chromosome_length=(int32_t)chromosomes.at(main_chromosome).length();
    const uint8_t BIN_SIZE = 5;  // Arbitrary, just for stats.
    bool previous_handle_exists, split_node;
    uint8_t orientation_cis, orientation_trans, is_bnd_single;
    size_t i, j, r, s, t;
    size_t n_positions, assign_record;
    int32_t p, q, p_prime, q_prime, first_pos, last_pos, trans_pos;
    string tmp_buffer, tmp_buffer_2, str;
    pair<int32_t,int32_t> tmp_pair;
    tuple<int32_t,int32_t,int32_t> tmp_tuple;
    handle_t previous_handle;

    // Collecting every distinct first position of a reference node, and building all non-reference nodes.
    // Remark: it can happen that the list of first positions of a chromosome contains the position right after the last
    // position of the chromosome (e.g. if there is a BND that connects the ends of two chromosomes, or if there is an
    // INS after the end of a chromosome).
    if (!silent) {
        cerr << "Number of VCF records: " << n_vcf_records << '\n';
        print_vcf_records_stats(BIN_SIZE,callers);
    }
    chunk_first_raw.clear(); chunk_first_raw.reserve(n_chromosomes);
    for (const auto& chromosome: chromosomes) chunk_first_raw.emplace(chromosome.first,vector<int32_t>());
    vector<int32_t>& first_positions = chunk_first_raw.at(main_chromosome);
    for (i=0; i<n_vcf_records; i++) {
        VcfRecord& record = vcf_records.at(i);
        const string& chrom_sequence = chromosomes.at(record.chrom);
        record.get_reference_coordinates(false,tmp_pair);
        if (tmp_pair.first==INT32_MAX || tmp_pair.second==INT32_MAX) continue;
        if (tmp_pair.first!=tmp_pair.second) {
            first_positions.push_back(tmp_pair.first);
            if (!acyclic || (record.sv_type!=VcfReader::TYPE_DUPLICATION && record.sv_type!=VcfReader::TYPE_CNV)) {
                if (tmp_pair.second>main_chromosome_length) tmp_pair.second=main_chromosome_length;
                first_positions.push_back(tmp_pair.second);
            }
            if (record.sv_type==VcfReader::TYPE_REPLACEMENT || record.sv_type==VcfReader::TYPE_SNP) {
                const handle_t handle_alt = sequence2handle(record.alt,tmp_buffer_2);
                insertion_handles.emplace_back(handle_alt);
            }
            else if (acyclic) {
                if (record.sv_type==VcfReader::TYPE_DUPLICATION || record.sv_type==VcfReader::TYPE_CNV) {
                    get<0>(tmp_tuple)=VcfReader::TYPE_DUPLICATION; get<1>(tmp_tuple)=tmp_pair.first; get<2>(tmp_tuple)=tmp_pair.second;
                    if (interval_to_insertion_handle.contains(tmp_tuple)) insertion_handles.emplace_back(interval_to_insertion_handle.at(tmp_tuple));
                    else {
                        const handle_t reference_handle = sequence2handle(chrom_sequence.substr(tmp_pair.first,tmp_pair.second-tmp_pair.first),tmp_buffer_2);
                        insertion_handles.emplace_back(reference_handle);
                        interval_to_insertion_handle.emplace(tmp_tuple,reference_handle);
                    }
                }
                else if (record.sv_type==VcfReader::TYPE_INVERSION) {
                    get<0>(tmp_tuple)=record.sv_type; get<1>(tmp_tuple)=tmp_pair.first; get<2>(tmp_tuple)=tmp_pair.second;
                    if (interval_to_insertion_handle.contains(tmp_tuple)) insertion_handles.emplace_back(interval_to_insertion_handle.at(tmp_tuple));
                    else {
                        str=chrom_sequence.substr(tmp_pair.first,tmp_pair.second-tmp_pair.first);
                        reverse_complement(str);
                        const handle_t reference_handle = sequence2handle(str,tmp_buffer_2);
                        insertion_handles.emplace_back(reference_handle);
                        interval_to_insertion_handle.emplace(tmp_tuple,reference_handle);
                    }
                }
            }
        }
        else {
            if (record.sv_type==VcfReader::TYPE_INSERTION && !record.is_symbolic) {
                first_positions.emplace_back(tmp_pair.first);
                const handle_t handle_ins = sequence2handle(record.alt.substr(1),tmp_buffer_2);
                insertion_handles.emplace_back(handle_ins);
            }
            else if (!acyclic) {
                is_bnd_single=record.is_breakend_single();
                if (record.sv_type==VcfReader::TYPE_BREAKEND && is_bnd_single>=1 && !record.is_breakend_virtual(chromosomes)) {
                    // Virtual telomeric breakends, and single breakends with no inserted sequence, carry no information.
                    r=record.get_info_field(VcfReader::MATEID_STR,0,tmp_buffer);
                    if (r!=string::npos && bnd_ids.contains(tmp_buffer)) {
                        // NOP: the breakend was already processed at its mate record.
                        // Remark: If the MATEID field is missing from a mate record, the non-reference nodes created
                        // are the same.
                    }
                    else {
                        bnd_ids.emplace(record.id);
                        orientation_cis=record.get_breakend_orientation_cis();
                        first_positions.emplace_back(orientation_cis==1?tmp_pair.first+1:tmp_pair.first);
                        if (is_bnd_single==2) {
                            trans_pos=record.get_breakend_pos();
                            trans_pos--; // Zero-based
                            orientation_trans=record.get_breakend_orientation_trans();
                            record.get_breakend_chromosome(tmp_buffer);
                            chunk_first_raw[tmp_buffer].emplace_back(orientation_trans==1?trans_pos+1:trans_pos);
                            record.get_breakend_inserted_sequence(tmp_buffer);
                            if (!tmp_buffer.empty()) {
                                const handle_t handle_ins = sequence2handle(tmp_buffer,tmp_buffer_2);
                                insertion_handles.emplace_back(handle_ins);
                            }
                        }
                        else {
                            record.get_breakend_inserted_sequence(tmp_buffer);
                            const handle_t handle_ins = sequence2handle(tmp_buffer,tmp_buffer_2);
                            insertion_handles.emplace_back(handle_ins);
                        }
                    }
                }
            }
        }
    }
    interval_to_insertion_handle.clear();
    for (auto& [chromosome,first_positions]: chunk_first_raw) sort_and_compact_positions(first_positions);
    if (!silent) {
        cerr << "Table: Chromosome involved | N. distinct breakpoints\n";
        for (auto& [chromosome,first_positions]: chunk_first_raw) {
            const uint32_t n_breakpoints = first_positions.size();
            if (n_breakpoints==0) continue;
            cerr << chromosome << ": ";
            for (const auto& position: first_positions) cerr << std::to_string(position) << ',';
            cerr << '\n';
            cerr << chromosome << ',' << first_positions.size() << '\n';
        }
        cerr << "Number of non-reference nodes: " << insertion_handles.size() << '\n';
    }

    // Building all reference nodes and all reference edges, taking into account the tandem repeat track and the
    // flanking requirements.
    // Remark: this procedure works also when the list of first positions of a chromosome contains the position right
    // after the last position of the chromosome.
    chunk_first.clear(); chunk_first.reserve(n_chromosomes);
    for (const auto& [chromosome,first_positions]: chunk_first_raw) {
        n_positions=first_positions.size();
        if (n_positions==0) continue;
        const string& chrom_sequence = chromosomes.at(chromosome);
        const auto chrom_length = (int32_t)chrom_sequence.size();
        vector<int32_t> first_positions_new;
        first_positions_new.reserve(n_positions);
        vector<handle_t>& handles = node_handles.at(chromosome);
        handles.clear();
        p=first_positions.at(0);
        first_pos=get_flank_boundary_left(chromosome,x!=INT32_MAX&&x<p&&chromosome==main_chromosome?x:p,flank_length);
        if (first_pos==INT32_MAX) first_pos=0;
        if (p!=first_pos) {
            const handle_t reference_handle = sequence2handle(chrom_sequence.substr(first_pos,p-first_pos),tmp_buffer_2);
            handles.emplace_back(reference_handle);
            first_positions_new.emplace_back(first_pos);
            node_to_chromosome[graph.get_id(reference_handle)]=pair<string,int32_t>(chromosome,first_pos);
            previous_handle_exists=true; previous_handle=reference_handle;
        }
        else previous_handle_exists=false;
        for (j=1; j<n_positions; j++) {
            p=first_positions.at(j-1);
            q=first_positions.at(j);
            split_node=false;
            if (q-p>=interior_flank_length) {
                p_prime=get_flank_boundary_right(chromosome,p,flank_length);
                q_prime=get_flank_boundary_left(chromosome,q,flank_length);
                if (p_prime!=INT32_MAX && q_prime!=INT32_MAX && p_prime<q_prime) {
                    split_node=true;
                    const handle_t reference_handle_1 = sequence2handle(chrom_sequence.substr(p,p_prime-p+1),tmp_buffer_2);
                    const handle_t reference_handle_2 = sequence2handle(chrom_sequence.substr(q_prime,q-q_prime),tmp_buffer_2);
                    handles.emplace_back(reference_handle_1);
                    handles.emplace_back(reference_handle_2);
                    first_positions_new.emplace_back(p);
                    first_positions_new.emplace_back(q_prime);
                    node_to_chromosome[graph.get_id(reference_handle_1)]=pair<string,int32_t>(chromosome,p);
                    node_to_chromosome[graph.get_id(reference_handle_2)]=pair<string,int32_t>(chromosome,q_prime);
                    if (previous_handle_exists) graph.create_edge(previous_handle,reference_handle_1);
                    previous_handle_exists=true; previous_handle=reference_handle_2;
                }
            }
            if (!split_node) {
                const handle_t reference_handle = sequence2handle(chrom_sequence.substr(p,q-p),tmp_buffer_2);
                handles.emplace_back(reference_handle);
                first_positions_new.emplace_back(p);
                node_to_chromosome[graph.get_id(reference_handle)]=pair<string,int32_t>(chromosome,p);
                if (previous_handle_exists) graph.create_edge(previous_handle,reference_handle);
                previous_handle_exists=true; previous_handle=reference_handle;
            }
        }
        p=first_positions.at(n_positions-1);
        first_positions_new.emplace_back(p);
        if (p<chrom_length) {
            last_pos=get_flank_boundary_right(chromosome,y!=INT32_MAX&&y>p&&chromosome==main_chromosome?y:p,flank_length);
            if (last_pos==INT32_MAX) last_pos=(int32_t)(chrom_length-1);
            const handle_t reference_handle = sequence2handle(chrom_sequence.substr(p,last_pos+1-p),tmp_buffer_2);
            handles.emplace_back(reference_handle);
            node_to_chromosome[graph.get_id(reference_handle)]=pair<string,int32_t>(chromosome,p);
            if (previous_handle_exists) graph.create_edge(previous_handle,reference_handle);
        }
        chunk_first.emplace(chromosome,first_positions_new);
    }
    chunk_first_raw.clear();
    const size_t n_reference_edges = graph.get_edge_count();
    if (!silent) {
        cerr << "Number of reference nodes: " << to_string(graph.get_node_count()-insertion_handles.size()) << '\n';
        cerr << "Total number of nodes: " << to_string(graph.get_node_count()) << '\n';
        cerr << "Table: Chromosome involved | First position of every reference node\n";
        for (const auto& [chromosome,first_positions]: chunk_first) {
            const size_t n_breakpoints = first_positions.size();
            if (n_breakpoints==0) continue;
            cerr << chromosome << ": ";
            for (const auto& position: first_positions) cerr << std::to_string(position) << ',';
            cerr << '\n';
        }
        cerr << "Table: Chromosome involved | N. reference nodes\n";
        for (const auto& [chromosome,chunks]: node_handles) {
            const size_t n_chunks = chunks.size();
            if (n_chunks>0) cerr << chromosome << ',' << n_chunks << '\n';
        }
        cerr << "Number of reference edges: " << n_reference_edges << '\n';
    }

    // Building all non-reference edges
    vcf_record_to_edge.reserve(n_vcf_records);
    for (i=0; i<vcf_record_to_edge.size(); i++) vcf_record_to_edge.at(i).clear();
    for (i=vcf_record_to_edge.size(); i<n_vcf_records; i++) vcf_record_to_edge.emplace_back();
    const vector<handle_t>& handles = node_handles.at(main_chromosome);
    bnd_ids.clear(); j=0;
    for (i=0; i<n_vcf_records; i++) {
        VcfRecord& record = vcf_records.at(i);
        record.get_reference_coordinates(false,tmp_pair);
        if (tmp_pair.first==INT32_MAX || tmp_pair.second==INT32_MAX) continue;
        assign_record=record.sv_length>=min_variant_length?i:SIZE_MAX;
        if (tmp_pair.first!=tmp_pair.second) {
            s=find_closest(main_chromosome,tmp_pair.first);
            t=find_closest(main_chromosome,tmp_pair.second);
            if (record.sv_type==VcfReader::TYPE_DELETION) {
                if (s!=0 && tmp_pair.second!=main_chromosome_length) {
                    add_nonreference_edge(handles.at(s-1),handles.at(t),assign_record);
                    if (assign_record!=SIZE_MAX) vcf_record_to_edge.at(assign_record).emplace_back(null_edge);
                }
            }
            else if (record.sv_type==VcfReader::TYPE_DUPLICATION || record.sv_type==VcfReader::TYPE_CNV) {
                if (acyclic) {
                    const handle_t& insertion_handle = insertion_handles.at(j);
                    if (s!=0) add_nonreference_edge(handles.at(s-1),insertion_handle,assign_record);
                    if (tmp_pair.first!=main_chromosome_length) add_nonreference_edge(insertion_handle,handles.at(s),assign_record);
                    if (assign_record!=SIZE_MAX) vcf_record_to_edge.at(assign_record).emplace_back(null_edge);
                    j++;
                }
                else {
                    add_nonreference_edge(handles.at(t-1),handles.at(s),assign_record);
                    if (assign_record!=SIZE_MAX) vcf_record_to_edge.at(assign_record).emplace_back(null_edge);
                }
            }
            else if (record.sv_type==VcfReader::TYPE_INVERSION) {
                if (acyclic) {
                    const handle_t& insertion_handle = insertion_handles.at(j);
                    if (s!=0) add_nonreference_edge(handles.at(s-1),insertion_handle,assign_record);
                    if (tmp_pair.second!=main_chromosome_length) add_nonreference_edge(insertion_handle,handles.at(t),assign_record);
                    if (assign_record!=SIZE_MAX) vcf_record_to_edge.at(assign_record).emplace_back(null_edge);
                    j++;
                }
                else {
                    if (s!=0) add_nonreference_edge(handles.at(s-1),graph.flip(handles.at(t-1)),assign_record);
                    if (tmp_pair.second!=main_chromosome_length) add_nonreference_edge(graph.flip(handles.at(s)),handles.at(t),assign_record);
                    if (assign_record!=SIZE_MAX) vcf_record_to_edge.at(assign_record).emplace_back(null_edge);
                }
            }
            else if (record.sv_type==VcfReader::TYPE_REPLACEMENT || record.sv_type==VcfReader::TYPE_SNP) {
                const handle_t& insertion_handle = insertion_handles.at(j);
                if (s!=0) add_nonreference_edge(handles.at(s-1),insertion_handle,assign_record);
                if (tmp_pair.second!=main_chromosome_length) add_nonreference_edge(insertion_handle,handles.at(t),assign_record);
                if (assign_record!=SIZE_MAX) vcf_record_to_edge.at(assign_record).emplace_back(null_edge);
                j++;
            }
        }
        else {
            if (record.sv_type==VcfReader::TYPE_INSERTION && !record.is_symbolic) {
                s=find_closest(main_chromosome,tmp_pair.first);
                const handle_t& insertion_handle = insertion_handles.at(j);
                if (s!=0) add_nonreference_edge(handles.at(s-1),insertion_handle,assign_record);
                if (tmp_pair.first!=main_chromosome_length) add_nonreference_edge(insertion_handle,handles.at(s),assign_record);
                if (assign_record!=SIZE_MAX) vcf_record_to_edge.at(assign_record).emplace_back(null_edge);
                j++;
            }
            else if (!acyclic) {
                is_bnd_single=record.is_breakend_single();
                if (record.sv_type==VcfReader::TYPE_BREAKEND && is_bnd_single>=1 && !record.is_breakend_virtual(chromosomes)) {
                    // Virtual telomeric breakends, and single breakends with no inserted sequence, carry no information.
                    r=record.get_info_field(VcfReader::MATEID_STR,0,tmp_buffer);
                    if (r!=string::npos && bnd_ids.contains(tmp_buffer)) {
                        // NOP: the breakend was already processed at its mate record.
                        // Remark: If the MATEID field is missing from a mate record, the non-reference edges created
                        // are the same.
                    }
                    else {
                        bnd_ids.emplace(record.id);
                        orientation_cis=record.get_breakend_orientation_cis();
                        const handle_t& handle_from = orientation_cis==1?handles.at(find_closest(main_chromosome,tmp_pair.first+1)-1):graph.flip(handles.at(find_closest(main_chromosome,tmp_pair.first)));
                        if (is_bnd_single==2) {
                            record.get_breakend_chromosome(tmp_buffer);
                            trans_pos=record.get_breakend_pos()-1;
                            orientation_trans=record.get_breakend_orientation_trans();
                            const handle_t& handle_to = orientation_trans==1?graph.flip(node_handles.at(tmp_buffer).at(find_closest(tmp_buffer,trans_pos+1)-1)):node_handles.at(tmp_buffer).at(find_closest(tmp_buffer,trans_pos));
                            record.get_breakend_inserted_sequence(tmp_buffer);
                            if (tmp_buffer.empty()) {
                                add_nonreference_edge(handle_from,handle_to,assign_record);
                                if (assign_record!=SIZE_MAX) vcf_record_to_edge.at(assign_record).emplace_back(null_edge);
                            }
                            else {
                                const handle_t& insertion_handle = orientation_cis==1?insertion_handles.at(j):graph.flip(insertion_handles.at(j));
                                add_nonreference_edge(handle_from,insertion_handle,assign_record);
                                add_nonreference_edge(insertion_handle,handle_to,assign_record);
                                if (assign_record!=SIZE_MAX) vcf_record_to_edge.at(assign_record).emplace_back(null_edge);
                                j++;
                            }
                        }
                        else {
                            record.get_breakend_inserted_sequence(tmp_buffer);
                            const handle_t& insertion_handle = orientation_cis==1?insertion_handles.at(j):graph.flip(insertion_handles.at(j));
                            add_nonreference_edge(handle_from,insertion_handle,assign_record);
                            if (assign_record!=SIZE_MAX) vcf_record_to_edge.at(assign_record).emplace_back(null_edge);
                            j++;
                        }
                    }
                }
            }
        }
    }
    if (!silent) {
        cerr << "Total number of edges: " << to_string(graph.get_edge_count()) << '\n';
        cerr << "Number of non-reference edges: " << to_string(graph.get_edge_count()-n_reference_edges) << '\n';
        cerr << "Number of chrA-chrB edges: " << to_string(get_n_interchromosomal_edges()) << '\n';
        cerr << "Number of chr-INS edges: " << to_string(get_n_ins_edges()) << '\n';
        cerr << "Number of nodes with edges on just one side: " << to_string(get_n_dangling_nodes()) << '\n';
        unordered_set<nid_t> dangling_nodes;
        get_dangling_nodes(false,dangling_nodes);
        print_edge_histograms();
    }

    // Transforming `insertion_handles` into `insertion_handles_set`.
    insertion_handles_set.clear(); insertion_handles_set.reserve(insertion_handles.size());
    insertion_handles_set.insert(insertion_handles.begin(),insertion_handles.end());
    insertion_handles.clear();

    // Graph closure
    if (graph_closure) build_graph_closure(acyclic);

    // Deallocating temporary space.
    chunk_first.clear(); node_handles.clear(); bnd_ids.clear();
    if (deallocate_ref_alt) {  // Deallocating REF and ALT
        for (auto& record: vcf_records) { record.ref.clear(); record.alt.clear(); }
    }

    // Allocating temporary space: `printed`, `initialized`, `flags`.
    mark_redundant_records();
    printed.clear(); printed.reserve(n_vcf_records);
    for (i=0; i<n_vcf_records; i++) printed.emplace_back(false);
    path_names.clear(); path_names.reserve(n_vcf_records);
    for (i=0; i<n_vcf_records; i++) path_names.emplace_back();
    initialized.clear(); initialized.reserve(n_vcf_records);
    for (i=0; i<n_vcf_records; i++) initialized.emplace_back(false);
    flags.reserve(n_vcf_records);
    for (i=flags.size(); i<n_vcf_records; i++) flags.emplace_back();
    for (i=0; i<n_vcf_records; i++) flags.at(i).reserve(vcf_record_to_edge.at(i).size());
}


/**
 * Closure is applied over multiple iterations, until no new edge is created.
 *
 * Remark: if an edge is processed in one iteration of closure, it is not skipped in the following iterations, since it
 * might get connected to new nodes.
 */
void VariantGraph::build_graph_closure(bool acyclic) {
    size_t i;
    vector<int32_t> tmp_vector;
    vector<tuple<handle_t,handle_t,edge_t,edge_t>> edge_instructions;
    vector<vector<edge_t>> vcf_record_to_edge_next;
    unordered_set<edge_t> new_edges;

    for (auto& pair: edge_to_vcf_record) {
        if (pair.second.size()>1) sort(pair.second.begin(),pair.second.end());
    }
    vcf_record_to_edge_next.reserve(n_vcf_records);
    for (i=0; i<n_vcf_records; i++) vcf_record_to_edge_next.emplace_back();
    while (true) {
        for (auto& pair: edge_to_vcf_record) {  // Every non-reference edge
            const edge_t& e1 = pair.first;
            if (is_reference_node(e1.second)) build_graph_closure_impl(e1,true,acyclic,edge_instructions,tmp_vector);
            if (is_reference_node(e1.first)) build_graph_closure_impl(e1,false,acyclic,edge_instructions,tmp_vector);
        }
        if (edge_instructions.empty()) break;
        for (auto& t: edge_instructions) {
            // `edge_instructions` is guaranteed to contain edges that are not in `graph`.
            graph.create_edge(std::get<0>(t),std::get<1>(t));
            edge_t new_edge = graph.edge_handle(std::get<0>(t),std::get<1>(t));
            new_edges.insert(new_edge);
            edge_t& e1 = std::get<2>(t);
            edge_t& e2 = std::get<3>(t);
            for (auto& record_id: edge_to_vcf_record[e1]) {
                if (!edge_to_vcf_record.contains(new_edge)) edge_to_vcf_record[new_edge]={record_id};
                else edge_to_vcf_record.at(new_edge).emplace_back(record_id);
                build_graph_closure_update_vcf_record_to_edge(record_id,e1,new_edge,vcf_record_to_edge_next);
            }
            for (auto& record_id: edge_to_vcf_record[e2]) {
                if (!edge_to_vcf_record.contains(new_edge)) edge_to_vcf_record[new_edge]={record_id};
                else edge_to_vcf_record.at(new_edge).emplace_back(record_id);
                build_graph_closure_update_vcf_record_to_edge(record_id,e2,new_edge,vcf_record_to_edge_next);
            }
        }
        for (auto& edge: new_edges) {
            if (edge_to_vcf_record.at(edge).size()>1) sort(edge_to_vcf_record.at(edge).begin(),edge_to_vcf_record.at(edge).end());
        }
        new_edges.clear(); edge_instructions.clear();
        for (i=0; i<n_vcf_records; i++) {
            if (!vcf_record_to_edge_next.at(i).empty()) vcf_record_to_edge.at(i).insert(vcf_record_to_edge.at(i).end(),make_move_iterator(vcf_record_to_edge_next.at(i).begin()),make_move_iterator(vcf_record_to_edge_next.at(i).end()));
        }
        for (i=0; i<n_vcf_records; i++) vcf_record_to_edge_next.at(i).clear();
    }
    build_graph_closure_close_vcf_record_to_edge();
}


void VariantGraph::build_graph_closure_impl(const edge_t& e1, bool orientation, bool acyclic, vector<tuple<handle_t,handle_t,edge_t,edge_t>>& new_edges, vector<int32_t>& tmp_pos) {
    const handle_t& from = orientation?e1.first:graph.flip(e1.second);
    const handle_t& to = orientation?e1.second:graph.flip(e1.first);

    if (!graph.get_is_reverse(to)) {
        const handle_t& source = get_previous_reference_node(to);
        if (source!=to) {
            graph.follow_edges(source,false,[&](handle_t new_neighbor) {
                const edge_t e2 = graph.edge_handle(source,new_neighbor);
                if (build_graph_closure_should_create_edge(e1,e2,acyclic,tmp_pos) && !graph.has_edge(from,new_neighbor)) {
                    new_edges.emplace_back(from,new_neighbor,e1,e2);
                    return;
                }
            });
        }
    }
    else {
        const handle_t& source = get_next_reference_node(to);
        if (source!=to) {
            graph.follow_edges(source,true,[&](handle_t new_neighbor) {
                const edge_t e2 = graph.edge_handle(graph.flip(source),new_neighbor);
                if (build_graph_closure_should_create_edge(e1,e2,acyclic,tmp_pos) && !graph.has_edge(from,new_neighbor)) {
                    new_edges.emplace_back(from,new_neighbor,e1,e2);
                    return;
                }
            });
        }
    }
}


bool VariantGraph::build_graph_closure_should_create_edge(const edge_t& e1, const edge_t& e2, bool acyclic, vector<int32_t>& tmp_pos) {
    bool i1, d1;
    bool is_insertion_1, is_insertion_2, is_duplication_1, is_duplication_2;
    size_t p, q;

    // e2 must be non-ref
    if (!edge_to_vcf_record.contains(e2)) return false;

    // e1 and e2 cannot share a VCF record
    vector<size_t>& v1 = edge_to_vcf_record.at(e1);
    vector<size_t>& v2 = edge_to_vcf_record.at(e2);
    p=0; q=0;
    while (p<v1.size() && q<v2.size()) {
        if (v1.at(p)<v2.at(q)) p++;
        else if (v1.at(p)>v2.at(q)) q++;
        else return false;
    }

    // Collecting info on e1
    is_insertion_1=false; is_duplication_1=false; tmp_pos.clear();
    for (auto& record_id: edge_to_vcf_record.at(e1)) {
        VcfRecord& record = vcf_records.at(record_id);
        d1 = record.sv_type==VcfReader::TYPE_DUPLICATION || record.sv_type==VcfReader::TYPE_CNV;
        i1 = record.sv_type==VcfReader::TYPE_INSERTION || (acyclic && d1);
        if (i1) { is_insertion_1=true; tmp_pos.emplace_back(record.pos); }
        if (d1) is_duplication_1=true;
    }
    sort_and_compact_positions(tmp_pos);

    // INS/DUP logic
    if (!acyclic) {
        if (is_duplication_1) return false;
        for (auto& record_id: edge_to_vcf_record.at(e2)) {
            VcfRecord& record = vcf_records.at(record_id);
            if (record.sv_type==VcfReader::TYPE_DUPLICATION || record.sv_type==VcfReader::TYPE_CNV) return false;
        }
    }
    if (is_insertion_1) {
        for (auto& record_id: edge_to_vcf_record.at(e2)) {
            VcfRecord& record = vcf_records.at(record_id);
            is_duplication_2 = record.sv_type==VcfReader::TYPE_DUPLICATION || record.sv_type==VcfReader::TYPE_CNV;
            is_insertion_2 = record.sv_type==VcfReader::TYPE_INSERTION || (acyclic && is_duplication_2);
            if (is_insertion_2 && binary_search(tmp_pos.begin(),tmp_pos.end(),record.pos)) return false;
        }
    }
    return true;
}


void VariantGraph::build_graph_closure_update_vcf_record_to_edge(int32_t vcf_record, const edge_t& old_edge, const edge_t& new_edge, vector<vector<edge_t>>& vcf_record_to_edge_next) {
    bool found;
    size_t i, j;
    size_t first, n_edges;

    const vector<edge_t>& old_edges = vcf_record_to_edge.at(vcf_record);
    vector<edge_t>& new_edges = vcf_record_to_edge_next.at(vcf_record);
    n_edges=old_edges.size();
    first=0;
    for (i=0; i<n_edges; i++) {
        if (old_edges.at(i)!=null_edge) continue;
        found=false;
        for (j=first; j<i; j++) {
            if (old_edges.at(j)==old_edge) { found=true; break; }
        }
        if (found) {
            for (j=first; j<=i; j++) {
                if (old_edges.at(j)==old_edge) new_edges.emplace_back(new_edge);
                else new_edges.emplace_back(old_edges.at(j));
            }
        }
        first=i+1;
    }
}


void VariantGraph::build_graph_closure_close_vcf_record_to_edge() {
    bool found;
    int32_t i, j, k, h;
    size_t size, node_id, first, second, first_prime, second_prime;
    vector<edge_t> in_edges, out_edges, new_pairs;

    for (i=0; i<n_vcf_records; i++) {
        size=vcf_record_to_edge.at(i).size();
        if (size<3 || vcf_record_to_edge.at(i).at(2)!=null_edge) continue;
        // Determining the node that is traversed by the record
        first=graph.get_id(vcf_record_to_edge.at(i).at(0).first);
        second=graph.get_id(vcf_record_to_edge.at(i).at(0).second);
        first_prime=graph.get_id(vcf_record_to_edge.at(i).at(1).first);
        second_prime=graph.get_id(vcf_record_to_edge.at(i).at(1).second);
        if (first_prime==first || first_prime==second) node_id=first_prime;
        else if (second_prime==first || second_prime==second) node_id=second_prime;
        else throw runtime_error("ERROR: The following VCF record does not traverse a node: "+to_string(i));
        // Partitioning edges (in canonical form) into in- and out-edges WRT the traversed node.
        in_edges.clear(); out_edges.clear(); first=0;
        for (j=0; j<size; j++) {
            if (vcf_record_to_edge.at(i).at(j)==null_edge) {
                if (j!=first+2) throw runtime_error("ERROR: The following VCF record has a sequence of "+to_string(j-first)+" edges: "+to_string(i));
                first=j+1;
                continue;
            }
            handle_t& first_handle = vcf_record_to_edge.at(i).at(j).first;
            handle_t& second_handle = vcf_record_to_edge.at(i).at(j).second;
            if (graph.get_id(first_handle)==node_id) {
                if (!graph.get_is_reverse(first_handle)) out_edges.emplace_back(vcf_record_to_edge.at(i).at(j));
                else in_edges.emplace_back(vcf_record_to_edge.at(i).at(j));
            }
            else if (graph.get_id(second_handle)==node_id) {
                if (!graph.get_is_reverse(first_handle)) in_edges.emplace_back(vcf_record_to_edge.at(i).at(j));
                else out_edges.emplace_back(vcf_record_to_edge.at(i).at(j));
            }
            else throw runtime_error("ERROR: The following VCF record has a sequence of edges that does not traverse its designated node "+to_string(node_id)+": "+to_string(i));
        }
        // Adding every new pair of in-edge and out-edge
        new_pairs.clear();
        for (j=0; j<in_edges.size(); j++) {
            edge_t& edge1 = in_edges.at(j);
            for (k=0; k<out_edges.size(); k++) {
                edge_t& edge2 = out_edges.at(k);
                found=false;
                for (h=0; h<size; h+=3) {
                    edge_t& edge3 = vcf_record_to_edge.at(i).at(h);
                    edge_t& edge4 = vcf_record_to_edge.at(i).at(h+1);
                    if ((edge3==edge1 && edge4==edge2) || (edge3==edge2 && edge4==edge1)) { found=true; break; }
                }
                if (!found) { new_pairs.emplace_back(edge1); new_pairs.emplace_back(edge2); }
            }
        }
        size=new_pairs.size();
        for (j=0; j<size; j+=2) {
            vcf_record_to_edge.at(i).emplace_back(new_pairs.at(j));
            vcf_record_to_edge.at(i).emplace_back(new_pairs.at(j+1));
            vcf_record_to_edge.at(i).emplace_back(null_edge);
        }
    }
}


handle_t VariantGraph::get_previous_reference_node(const handle_t& node_handle) const {
    const nid_t node_id = graph.get_id(node_handle);
    const vector<handle_t>& handles_of_chromosome = node_handles.at(node_to_chromosome.at(node_id).first);
    const size_t n_handles = handles_of_chromosome.size();
    size_t i;

    for (i=1; i<n_handles; i++) {
        if (graph.get_id(handles_of_chromosome.at(i))==node_id && graph.has_edge(handles_of_chromosome.at(i-1),handles_of_chromosome.at(i))) return handles_of_chromosome.at(i-1);
    }
    return node_handle;
}


handle_t VariantGraph::get_next_reference_node(const handle_t& node_handle) const {
    const nid_t node_id = graph.get_id(node_handle);
    const vector<handle_t>& handles_of_chromosome = node_handles.at(node_to_chromosome.at(node_id).first);
    const auto n_handles = (int32_t)handles_of_chromosome.size();
    int32_t i;

    for (i=n_handles-2; i>=0; i--) {
        if (graph.get_id(handles_of_chromosome.at(i))==node_id && graph.has_edge(handles_of_chromosome.at(i),handles_of_chromosome.at(i+1))) return handles_of_chromosome.at(i+1);
    }
    return node_handle;
}


void VariantGraph::build(const string& chromosome, int32_t p, int32_t q, int32_t flank_length) {
    if (!chromosomes.contains(chromosome)) throw runtime_error("Invalid chromosome "+chromosome);
    main_chromosome=chromosome;
    const string& main_chromosome_sequence = chromosomes.at(main_chromosome);
    main_chromosome_length=(int32_t)main_chromosome_sequence.length();
    if (p<0 || p>=main_chromosome_length) throw runtime_error("Invalid p="+to_string(p)+" in "+chromosome+" of length "+to_string(main_chromosome_length));
    if (q<0 || q>main_chromosome_length) throw runtime_error("Invalid q="+to_string(q)+" in "+chromosome+" of length "+to_string(main_chromosome_length));
    if (p>q) throw runtime_error("Invalid p="+to_string(p)+" and q="+to_string(q)+" in "+chromosome+" of length "+to_string(main_chromosome_length));

    graph.clear();
    n_vcf_records=0; vcf_records.clear();
    bnd_ids.clear(); node_handles.clear();
    for (const auto& chromosome: chromosomes) node_handles.emplace(chromosome.first,vector<handle_t>());
    node_to_chromosome.clear(); insertion_handles.clear(); insertion_handles_set.clear();
    edge_to_vcf_record.clear(); vcf_record_to_edge.clear();

    bool previous_handle_exists;
    int32_t first_pos, last_pos;
    string tmp_buffer;
    handle_t previous_handle;
    vector<int32_t> first_positions;
    handle_t reference_handle;

    vector<handle_t>& handles = node_handles.at(chromosome);

    // Left flank
    first_pos=get_flank_boundary_left(chromosome,p,flank_length);
    if (first_pos==INT32_MAX) first_pos=0;
    if (p!=first_pos) {
        reference_handle=sequence2handle(main_chromosome_sequence.substr(first_pos,p-first_pos),tmp_buffer);
        handles.emplace_back(reference_handle);
        first_positions.emplace_back(first_pos);
        node_to_chromosome[graph.get_id(reference_handle)]=pair<string,int32_t>(chromosome,first_pos);
        previous_handle_exists=true; previous_handle=reference_handle;
    }
    else previous_handle_exists=false;

    // [p..q), if any.
    if (p<q) {
        reference_handle=sequence2handle(main_chromosome_sequence.substr(p,q-p),tmp_buffer);
        handles.emplace_back(reference_handle);
        first_positions.emplace_back(p);
        node_to_chromosome[graph.get_id(reference_handle)]=pair<string,int32_t>(chromosome,p);
        if (previous_handle_exists) graph.create_edge(previous_handle,reference_handle);
        previous_handle=reference_handle;
    }

    // Right flank
    if (q<main_chromosome_length) {
        last_pos=get_flank_boundary_right(chromosome,q,flank_length);
        if (last_pos==INT32_MAX) last_pos=(int32_t)(main_chromosome_length-1);
        reference_handle=sequence2handle(main_chromosome_sequence.substr(q,last_pos+1-q),tmp_buffer);
        handles.emplace_back(reference_handle);
        first_positions.emplace_back(q);
        node_to_chromosome[graph.get_id(reference_handle)]=pair<string,int32_t>(chromosome,q);
        if (previous_handle_exists) graph.create_edge(previous_handle,reference_handle);
    }

    chunk_first.clear(); chunk_first.emplace(chromosome,first_positions);
}


/**
 * Implemented with naive quadratic scans. Should be made faster.
 */
void VariantGraph::mark_redundant_records() {
    size_t i, j, k;
    size_t n_edges, n_sequences_i, n_sequences_j;

    for (auto& record: vcf_records) record.is_redundant=false;
    for (i=0; i<n_vcf_records; i++) {
        if (vcf_records.at(i).is_redundant) continue;
        n_edges=vcf_record_to_edge.at(i).size();
        n_sequences_i=0;
        for (k=0; k<n_edges; k++) {
            if (vcf_record_to_edge.at(i).at(k)==null_edge) n_sequences_i++;
        }
        for (j=i+1; j<n_vcf_records; j++) {
            if (vcf_records.at(j).is_redundant || n_edges!=vcf_record_to_edge.at(j).size()) continue;
            n_sequences_j=0;
            for (k=0; k<n_edges; k++) {
                if (vcf_record_to_edge.at(j).at(k)==null_edge) n_sequences_j++;
            }
            if (n_sequences_j!=n_sequences_i) continue;
            if (mark_redundant_records_impl(i,j,n_edges,n_sequences_i) && mark_redundant_records_impl(j,i,n_edges,n_sequences_i)) vcf_records.at(j).is_redundant=true;
        }
    }
}


bool VariantGraph::mark_redundant_records_impl(int32_t i, int32_t j, size_t n_edges, size_t n_sequences) {
    bool match;
    size_t k, h, m;
    size_t first_i, first_j;
    size_t n_sequences_matched;

    first_i=0; first_j=0; n_sequences_matched=0;
    for (k=1; k<n_edges; k++) {
        if (vcf_record_to_edge.at(i).at(k)!=null_edge) continue;
        for (h=1; h<n_edges; h++) {
            if (vcf_record_to_edge.at(j).at(h)!=null_edge) continue;
            if (h-first_j!=k-first_i) continue;
            match=true;
            for (m=0; m<k-first_i; m++) {
                if (vcf_record_to_edge.at(i).at(first_i+m)!=vcf_record_to_edge.at(j).at(first_j+m)) { match=false; break; }
            }
            if (!match) {
                match=true;
                for (m=0; m<k-first_i; m++) {
                    if (vcf_record_to_edge.at(i).at(first_i+m)!=vcf_record_to_edge.at(j).at(h-1-m)) { match=false; break; }
                }
            }
            if (match) n_sequences_matched++;
            first_j=h+1;
        }
        first_i=k+1;
    }
    return n_sequences_matched==n_sequences;
}


bool VariantGraph::would_graph_be_nontrivial(vector<VcfRecord>& records) {
    if (records.empty()) return false;
    pair<int32_t,int32_t> tmp_pair;
    for (auto& record: records) {
        record.get_reference_coordinates(false,tmp_pair);
        if (tmp_pair.first==INT32_MAX || tmp_pair.second==INT32_MAX) continue;
        if ( (record.sv_type==VcfReader::TYPE_INSERTION && !record.is_symbolic) ||
             record.sv_type==VcfReader::TYPE_DELETION ||
             record.sv_type==VcfReader::TYPE_INVERSION ||
             record.sv_type==VcfReader::TYPE_DUPLICATION ||
             record.sv_type==VcfReader::TYPE_REPLACEMENT ||
             record.sv_type==VcfReader::TYPE_CNV ||
             (record.sv_type==VcfReader::TYPE_BREAKEND && record.is_breakend_single()>=1 && !record.is_breakend_virtual(chromosomes)) ||
             record.sv_type==VcfReader::TYPE_SNP
             ) return true;
    }
    return false;
}


void VariantGraph::to_gfa(const path& gfa_path) const {
    ofstream outstream(gfa_path.string());
    if (!outstream.good() || !outstream.is_open()) throw runtime_error("ERROR: file could not be written: " + gfa_path.string());
    graph.for_each_handle([&](handle_t handle) { outstream << GFA_NODE_CHAR << GFA_SEPARATOR << graph.get_id(handle) << GFA_SEPARATOR << graph.get_sequence(handle) << GFA_LINE_END; });
    graph.for_each_edge([&](edge_t edge) { outstream << GFA_LINK_CHAR << GFA_SEPARATOR << graph.get_id(edge.first) << GFA_SEPARATOR << (graph.get_is_reverse(edge.first)?GFA_REV_CHAR:GFA_FWD_CHAR) << GFA_SEPARATOR << graph.get_id(edge.second) << GFA_SEPARATOR << (graph.get_is_reverse(edge.second)?GFA_REV_CHAR:GFA_FWD_CHAR) << GFA_SEPARATOR << GFA_OVERLAP_FIELD << GFA_LINE_END; });
    graph.for_each_path_handle([&](path_handle_t path) {
        outstream << GFA_PATH_CHAR << GFA_SEPARATOR << graph.get_path_name(path) << GFA_SEPARATOR;
        const step_handle_t source = graph.path_begin(path);
        const step_handle_t destination = graph.path_back(path);
        outstream << graph.get_id(graph.get_handle_of_step(source)) << (graph.get_is_reverse(graph.get_handle_of_step(source))?GFA_REV_CHAR:GFA_FWD_CHAR);
        step_handle_t from = source;
        while (graph.has_next_step(from)) {
            const step_handle_t to = graph.get_next_step(from);
            if (to==destination) break;
            outstream << GFA_PATH_SEPARATOR << graph.get_id(graph.get_handle_of_step(from)) << (graph.get_is_reverse(graph.get_handle_of_step(from))?GFA_REV_CHAR:GFA_FWD_CHAR);
            from=to;
        }
        outstream << GFA_SEPARATOR << GFA_OVERLAP_FIELD << GFA_LINE_END;
    });
    outstream.close();
}


vector<string> VariantGraph::load_gfa(const path& gfa_path) {
    bool is_node, is_link, is_path, orientation_from, orientation_to;
    char c;
    uint8_t n_fields;
    nid_t id_from, id_to;
    string buffer, node_sequence, path_name, path_encoding;
    vector<string> node_ids;
    ifstream instream;

    // Loading node IDs
    instream.open(gfa_path.string());
    if (!instream.good() || !instream.is_open()) throw runtime_error("ERROR: could not read file: "+gfa_path.string());
    n_fields=0;
    while (instream.get(c)) {
        if (c!=GFA_SEPARATOR && c!=GFA_LINE_END) { buffer.push_back(c); continue; }
        n_fields++;
        if (n_fields==1) {
            if (buffer.at(0)!=GFA_NODE_CHAR) {
                instream.ignore(STREAMSIZE_MAX,GFA_LINE_END);
                n_fields=0;
            }
        }
        else if (n_fields==2) node_ids.emplace_back(buffer);
        else {  // No other field is used
            if (c!=GFA_LINE_END) instream.ignore(STREAMSIZE_MAX,GFA_LINE_END);
            n_fields=0;
        }
        buffer.clear();
    }
    instream.close();
    if (node_ids.size()>1) sort(node_ids.begin(),node_ids.end());

    // Loading graph
    destroy_paths(); graph.clear();
    instream.open(gfa_path.string());
    if (!instream.good() || !instream.is_open()) throw runtime_error("ERROR: could not read file: "+gfa_path.string());
    is_node=false; is_link=false; is_path=false; orientation_from=false; orientation_to=false; id_from=INT32_MAX; id_to=INT32_MAX;
    n_fields=0;
    while (instream.get(c)) {
        if (c!=GFA_SEPARATOR && c!=GFA_LINE_END) { buffer.push_back(c); continue; }
        n_fields++;
        if (n_fields==1) {
            is_node=buffer.at(0)==GFA_NODE_CHAR;
            is_link=buffer.at(0)==GFA_LINK_CHAR;
            is_path=buffer.at(0)==GFA_PATH_CHAR;
            if (!is_node && !is_link && !is_path) {
                instream.ignore(STREAMSIZE_MAX,GFA_LINE_END);
                buffer.clear(); n_fields=0;
                continue;
            }
        }
        else if (n_fields==2) {
            if (is_path) path_name=buffer;
            else {
                const auto iterator = lower_bound(node_ids.begin(),node_ids.end(),buffer);
                id_from=distance(node_ids.begin(),iterator)+1;  // Node IDs cannot be zero
            }
        }
        else if (n_fields==3) {
            if (is_node) node_sequence=buffer;
            else if (is_link) orientation_from=buffer.at(0)==GFA_FWD_CHAR;
            else if (is_path) path_encoding=buffer;
        }
        else if (n_fields==4) {
            if (is_link) {
                const auto iterator = lower_bound(node_ids.begin(),node_ids.end(),buffer);
                id_to=distance(node_ids.begin(),iterator)+1;  // Node IDs cannot be zero
            }
        }
        else if (n_fields==5) {
            if (is_link) orientation_to=buffer.at(0)==GFA_FWD_CHAR;
        }
        else {  /* No other field is used */ }
        if (c==GFA_LINE_END) {
            if (is_node) graph.create_handle(node_sequence,id_from);
            else if (is_link) graph.create_edge(graph.get_handle(id_from,!orientation_from),graph.get_handle(id_to,!orientation_to));
            else if (is_path) load_gfa_path(path_encoding,node_ids,path_name,buffer);
            n_fields=0;
        }
        buffer.clear();
    }
    instream.close();
    if (!silent) cerr << "Loaded " << to_string(graph.get_node_count()) << " nodes, " << to_string(graph.get_edge_count()) << " edges, " << to_string(graph.get_path_count()) << " paths from the GFA file.\n";
    return node_ids;
}


path_handle_t VariantGraph::load_gfa_path(const string& path_encoding, const vector<string>& node_ids, const string& path_name, string& buffer) {
    const size_t LENGTH = path_encoding.length();
    char c;
    nid_t node_id;

    path_handle_t path = graph.create_path_handle(path_name);
    buffer.clear();
    for (size_t i=0; i<LENGTH; i++) {
        c=path_encoding.at(i);
        if (c==GFA_PATH_SEPARATOR) continue;
        else if (c!=GFA_FWD_CHAR && c!=GFA_REV_CHAR) { buffer.push_back(c); continue; }
        if (node_ids.empty()) node_id=stoi(buffer);
        else {
            const auto iterator = lower_bound(node_ids.begin(),node_ids.end(),buffer);
            node_id=distance(node_ids.begin(),iterator)+1;  // Node IDs cannot be zero
        }
        graph.append_step(path,graph.get_handle(node_id,c==GFA_REV_CHAR));
        buffer.clear();
    }
    return path;
}


void VariantGraph::destroy_paths() {
    vector<path_handle_t> path_handles;
    graph.for_each_path_handle([&](const path_handle_t& path){ path_handles.push_back(path); });
    for (auto& path_handle: path_handles) graph.destroy_path(path_handle);
}


void VariantGraph::mark_edge(const edge_t& query, size_t rank, const vector<edge_t>& edges, vector<size_t>& flags) const {
    const edge_t canonized_query = graph.edge_handle(query.first,query.second);
    const size_t N_EDGES = edges.size();

    for (size_t i=0; i<N_EDGES; i++) {
        if (edges.at(i)==canonized_query) flags.at(i)=rank;
    }
}


void VariantGraph::for_each_vcf_record_with_supporting_paths(const function<void(size_t id, const VcfRecord& record, const vector<string>& supporting_paths)>& callback) {
    bool all_present, is_increasing, is_decreasing;
    size_t i, j;
    size_t rank, n_edges, first;

    for (i=0; i<n_vcf_records; i++) printed.at(i)=false;
    for (i=0; i<n_vcf_records; i++) path_names.at(i).clear();
    for (i=0; i<n_vcf_records; i++) initialized.at(i)=false;
    graph.for_each_path_handle([&](path_handle_t path) {
        step_handle_t from = graph.path_begin(path);  // `graph` breaks circular paths
        const step_handle_t last = graph.path_back(path);  // `graph` breaks circular paths

        // Initializing `flags`.
        while (from!=last) {
            const step_handle_t to = graph.get_next_step(from);
            const edge_t canonized_edge = graph.edge_handle(graph.get_handle_of_step(from),graph.get_handle_of_step(to));
            if (!edge_to_vcf_record.contains(canonized_edge)) { from=to; continue; }
            for (const size_t& r: edge_to_vcf_record.at(canonized_edge)) initialized.at(r)=false;
            from=to;
        }
        from=graph.path_begin(path);
        while (from!=last) {
            const step_handle_t to = graph.get_next_step(from);
            const edge_t canonized_edge = graph.edge_handle(graph.get_handle_of_step(from),graph.get_handle_of_step(to));
            if (!edge_to_vcf_record.contains(canonized_edge)) { from=to; continue; }
            for (const size_t& r: edge_to_vcf_record.at(canonized_edge)) {
                if (printed.at(r) || initialized.at(r)) continue;
                n_edges=vcf_record_to_edge.at(r).size();
                flags.at(r).clear();
                for (i=0; i<n_edges; i++) flags.at(r).emplace_back(0);
                initialized.at(r)=true;
            }
            from=to;
        }

        // Marking `flags`.
        from=graph.path_begin(path); rank=1;
        while (from!=last) {
            const step_handle_t to = graph.get_next_step(from);
            const edge_t canonized_edge = graph.edge_handle(graph.get_handle_of_step(from),graph.get_handle_of_step(to));
            if (!edge_to_vcf_record.contains(canonized_edge)) { from=to; rank++; continue; }
            for (const size_t& r: edge_to_vcf_record.at(canonized_edge)) {
                if (printed.at(r)) continue;
                mark_edge(canonized_edge,rank,vcf_record_to_edge.at(r),flags.at(r));
            }
            from=to; rank++;
        }

        // Finding supported VCF records
        from=graph.path_begin(path);
        while (from!=last) {
            const step_handle_t to = graph.get_next_step(from);
            const edge_t canonized_edge = graph.edge_handle(graph.get_handle_of_step(from),graph.get_handle_of_step(to));
            if (!edge_to_vcf_record.contains(canonized_edge)) { from=to; continue; }
            for (const size_t& r: edge_to_vcf_record.at(canonized_edge)) {
                if (printed.at(r)) continue;
                n_edges=vcf_record_to_edge.at(r).size();
                first=0;
                for (i=0; i<n_edges; i++) {
                    if (vcf_record_to_edge.at(r).at(i)!=null_edge) continue;
                    all_present=true;
                    for (j=first; j<i; j++) {
                        if (flags.at(r).at(j)==0) { all_present=false; break; }
                    }
                    if (all_present) {
                        is_increasing=true; is_decreasing=true;
                        for (j=first+1; j<i; j++) {
                            if (flags.at(r).at(j)>flags.at(r).at(j-1)) is_decreasing=false;
                            if (flags.at(r).at(j)<flags.at(r).at(j-1)) is_increasing=false;
                        }
                        if (is_increasing || is_decreasing) {
                            printed.at(r)=true;
                            path_names.at(r).emplace_back(graph.get_path_name(path));
                            break;
                        }
                    }
                    first=i+1;
                }
            }
            from=to;
        }
    });

    for (i=0; i<n_vcf_records; i++) callback(i,vcf_records.at(i),path_names.at(i));
}


void VariantGraph::add_gaf_path_to_graph(const string& alignment_name, const vector <pair<string,bool> >& path){
    // Create a new path in the variant graph
    auto p = graph.create_path_handle(alignment_name);

    // Iterate the path steps and append each step to the prev step
    for (const auto& [step_name, is_reverse]: path){
        // Convert GAF name string back into nid, and construct handle of correct orientation
        nid_t id = stoll(step_name);
        auto h = graph.get_handle(id,is_reverse);
        graph.append_step(p,h);
    }
}


void VariantGraph::print_supported_vcf_records(ofstream& supported, ofstream& unsupported, bool print_all_records, const vector<string>& callers) {
    const size_t N_CALLERS = callers.size();
    size_t i, j;
    vector<vector<size_t>> caller_count;

    for (i=0; i<N_CALLERS; i++) caller_count.emplace_back(vector<size_t> {0,0,0,0,0});
    if (print_all_records) {
        for (i=0; i<n_vcf_records; i++) {
            VcfRecord& record = vcf_records.at(i);
            record.print(supported); supported << '\n';
            if (!callers.empty()) {
                if (record.sv_type==VcfReader::TYPE_INSERTION) increment_caller_count(record,callers,caller_count,0);
                else if (record.sv_type==VcfReader::TYPE_DUPLICATION) increment_caller_count(record,callers,caller_count,1);
                else if (record.sv_type==VcfReader::TYPE_DELETION) increment_caller_count(record,callers,caller_count,2);
                else if (record.sv_type==VcfReader::TYPE_INVERSION) increment_caller_count(record,callers,caller_count,3);
                else if (record.sv_type==VcfReader::TYPE_BREAKEND) increment_caller_count(record,callers,caller_count,4);
            }
        }
    }
    else {
        for_each_vcf_record_with_supporting_paths([&](size_t id, const VcfRecord& record, const vector<string>& supporting_paths) {
            if (supporting_paths.empty()) { record.print(unsupported); unsupported << '\n'; }
            else {
                record.print(supported); supported << '\n';
                if (!callers.empty()) {
                    if (record.sv_type==VcfReader::TYPE_INSERTION) increment_caller_count(record,callers,caller_count,0);
                    else if (record.sv_type==VcfReader::TYPE_DUPLICATION) increment_caller_count(record,callers,caller_count,1);
                    else if (record.sv_type==VcfReader::TYPE_DELETION) increment_caller_count(record,callers,caller_count,2);
                    else if (record.sv_type==VcfReader::TYPE_INVERSION) increment_caller_count(record,callers,caller_count,3);
                    else if (record.sv_type==VcfReader::TYPE_BREAKEND) increment_caller_count(record,callers,caller_count,4);
                }
            }
        });
    }
    if (!silent && !callers.empty()) {
        cerr << "Histogram: Caller | #INS printed | #DUP printed | #DEL printed | #INV printed | #BND printed\n";
        for (i=0; i<caller_count.size(); i++) {
            cerr << callers.at(i);
            for (j=0; j<caller_count.at(i).size(); j++) cerr << ',' << to_string(caller_count.at(i).at(j));
            cerr << '\n';
        }
    }
}


void VariantGraph::get_main_chromosome(string& result) {
    result = main_chromosome;
}


void VariantGraph::get_region_of_node(const handle_t& h, Region& result){
    const nid_t n = graph.get_id(h);
    const bool is_reverse = graph.get_is_reverse(h);
    const auto iter = node_to_chromosome.find(n);

    result.name.clear();
    result.start = -1;
    result.stop = -1;
    if (iter != node_to_chromosome.end()){
        const auto l = graph.get_length(h);

        result.name = iter->second.first;
        result.start = iter->second.second;
        result.stop = result.start + int32_t(l);

        if (is_reverse){
            std::swap(result.start, result.stop);
        }
    }
}


void VariantGraph::paths_to_vcf_records(ofstream& outstream) {
    const string DEFAULT_QUAL = "60";  // Arbitrary
    const string DEFAULT_FILTER = VcfReader::PASS_STR;
    const string DEFAULT_FORMAT = "GT";
    const string DEFAULT_SAMPLE = "0/1";
    const string SUFFIX_REPLACEMENT = VcfReader::VCF_SEPARATOR+DEFAULT_QUAL+VcfReader::VCF_SEPARATOR+DEFAULT_FILTER+VcfReader::VCF_SEPARATOR+VcfReader::VCF_MISSING_CHAR_1+VcfReader::VCF_SEPARATOR+DEFAULT_FORMAT+VcfReader::VCF_SEPARATOR+DEFAULT_SAMPLE+'\n';
    const string SUFFIX_BND = VcfReader::VCF_SEPARATOR+DEFAULT_QUAL+VcfReader::VCF_SEPARATOR+DEFAULT_FILTER+VcfReader::VCF_SEPARATOR+"SVTYPE=BND"+VcfReader::VCF_SEPARATOR+DEFAULT_FORMAT+VcfReader::VCF_SEPARATOR+DEFAULT_SAMPLE+'\n';
    char pos_base;
    int32_t pos;
    string old_sequence, new_sequence, alt_sequence;

    graph.for_each_path_handle([&](path_handle_t path) {
        const step_handle_t source = graph.path_begin(path);  // `graph` breaks circular paths
        const step_handle_t destination = graph.path_back(path);  // `graph` breaks circular paths
        const handle_t source_handle = graph.get_handle_of_step(source);
        const handle_t destination_handle = graph.get_handle_of_step(destination);
        const nid_t id1 = graph.get_id(source_handle);
        const nid_t id2 = graph.get_id(destination_handle);
        if (node_to_chromosome.contains(id1) && node_to_chromosome.contains(id2)) {
            const pair<string,int32_t> source_coordinate = node_to_chromosome.at(id1);
            const pair<string,int32_t> destination_coordinate = node_to_chromosome.at(id2);
            const bool is_source_reversed = graph.get_is_reverse(source_handle);
            const bool is_destination_reversed = graph.get_is_reverse(destination_handle);
            new_sequence.clear();
            if ( source_coordinate.first==destination_coordinate.first &&
                 ( (!is_source_reversed && !is_destination_reversed && source_coordinate.second<destination_coordinate.second) ||
                   (is_source_reversed && is_destination_reversed && source_coordinate.second>destination_coordinate.second)
                   )
                 ) {
                // Source and destination are on the same chromosome and in the same orientation. We represent this as a
                // replacement record.
                if (!is_source_reversed) {
                    pos=source_coordinate.second+(int32_t)graph.get_length(source_handle);  // One-based
                    pos_base=chromosomes.at(source_coordinate.first).at(pos-1);
                    old_sequence=chromosomes.at(source_coordinate.first).substr(pos,destination_coordinate.second-pos);
                    step_handle_t from = source;
                    while (graph.has_next_step(from)) {
                        const step_handle_t to = graph.get_next_step(from);
                        if (to==destination) break;
                        new_sequence.append(graph.get_sequence(graph.get_handle_of_step(to)));
                        from=to;
                    }
                }
                else {
                    pos=destination_coordinate.second+(int32_t)graph.get_length(destination_handle);  // One-based
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
                outstream << source_coordinate.first << VcfReader::VCF_SEPARATOR << pos << VcfReader::VCF_SEPARATOR << graph.get_path_name(path) << VcfReader::VCF_SEPARATOR << pos_base << old_sequence << VcfReader::VCF_SEPARATOR << pos_base << new_sequence << SUFFIX_REPLACEMENT;
            }
            else {
                // Source and destination are not on the same chromosome, or not in the same orientation. We represent
                // this as a BND record with inserted sequence.
                if (!is_source_reversed) {
                    pos=source_coordinate.second+(int32_t)graph.get_length(source_handle);  // One-based
                    pos_base=chromosomes.at(source_coordinate.first).at(pos-1);
                    step_handle_t from = source;
                    while (graph.has_next_step(from)) {
                        const step_handle_t to = graph.get_next_step(from);
                        if (to==destination) break;
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
                    outstream << source_coordinate.first << VcfReader::VCF_SEPARATOR << pos << VcfReader::VCF_SEPARATOR << graph.get_path_name(path) << VcfReader::VCF_SEPARATOR << pos_base << VcfReader::VCF_SEPARATOR << alt_sequence << SUFFIX_BND;
                }
                else {
                    pos=source_coordinate.second+1;  // One-based
                    pos_base=chromosomes.at(source_coordinate.first).at(pos-1);
                    step_handle_t from = source;
                    while (graph.has_next_step(from)) {
                        const step_handle_t to = graph.get_next_step(from);
                        if (to==destination) break;
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
                    outstream << source_coordinate.first << VcfReader::VCF_SEPARATOR << pos << VcfReader::VCF_SEPARATOR << graph.get_path_name(path) << VcfReader::VCF_SEPARATOR << pos_base << VcfReader::VCF_SEPARATOR << alt_sequence << SUFFIX_BND;
                }
            }
        }
        else if (node_to_chromosome.contains(id1)) {
            // We represent this as a "single BND" record with inserted sequence.
            const pair<string,int32_t> source_coordinate = node_to_chromosome.at(id1);
            const bool is_source_reversed = graph.get_is_reverse(source_handle);
            new_sequence.clear();
            if (!is_source_reversed) {
                pos=source_coordinate.second+(int32_t)graph.get_length(source_handle);  // One-based
                pos_base=chromosomes.at(source_coordinate.first).at(pos-1);
                step_handle_t from = source;
                while (graph.has_next_step(from)) {
                    const step_handle_t to = graph.get_next_step(from);
                    new_sequence.append(graph.get_sequence(graph.get_handle_of_step(to)));
                    if (to==destination) break;
                    from=to;
                }
                outstream << source_coordinate.first << VcfReader::VCF_SEPARATOR << pos << VcfReader::VCF_SEPARATOR << graph.get_path_name(path) << VcfReader::VCF_SEPARATOR << pos_base << VcfReader::VCF_SEPARATOR << pos_base << new_sequence << VcfReader::VCF_MISSING_CHAR_1 << SUFFIX_BND;
            }
            else {
                pos=source_coordinate.second+1;  // One-based
                pos_base=chromosomes.at(source_coordinate.first).at(pos-1);
                step_handle_t from = source;
                while (graph.has_next_step(from)) {
                    const step_handle_t to = graph.get_next_step(from);
                    new_sequence.append(graph.get_sequence(graph.get_handle_of_step(to)));
                    if (to==destination) break;
                    from=to;
                }
                reverse_complement(new_sequence);
                outstream << source_coordinate.first << VcfReader::VCF_SEPARATOR << pos << VcfReader::VCF_SEPARATOR << graph.get_path_name(path) << VcfReader::VCF_SEPARATOR << pos_base << VcfReader::VCF_SEPARATOR << VcfReader::VCF_MISSING_CHAR_1 << new_sequence << pos_base << SUFFIX_BND;
            }
        }
        else if (node_to_chromosome.contains(id2)) {
            // We represent this as a "single BND" record with inserted sequence.
            const pair<string,int32_t> destination_coordinate = node_to_chromosome.at(id2);
            const bool is_destination_reversed = graph.get_is_reverse(destination_handle);
            new_sequence.clear();
            if (!is_destination_reversed) {
                pos=destination_coordinate.second+1;  // One-based
                pos_base=chromosomes.at(destination_coordinate.first).at(pos-1);
                step_handle_t from = destination;
                while (graph.has_previous_step(from)) {
                    const step_handle_t to = graph.get_previous_step(from);
                    new_sequence.append(graph.get_sequence(graph.flip(graph.get_handle_of_step(to))));
                    if (to==source) break;
                    from=to;
                }
                reverse_complement(new_sequence);
                outstream << destination_coordinate.first << VcfReader::VCF_SEPARATOR << pos << VcfReader::VCF_SEPARATOR << graph.get_path_name(path) << VcfReader::VCF_SEPARATOR << pos_base << VcfReader::VCF_SEPARATOR << VcfReader::VCF_MISSING_CHAR_1 << new_sequence << pos_base << SUFFIX_BND;
            }
            else {
                pos=destination_coordinate.second+(int32_t)graph.get_length(destination_handle);  // One-based
                pos_base=chromosomes.at(destination_coordinate.first).at(pos-1);
                step_handle_t from = destination;
                while (graph.has_previous_step(from)) {
                    const step_handle_t to = graph.get_previous_step(from);
                    new_sequence.append(graph.get_sequence(graph.flip(graph.get_handle_of_step(to))));
                    if (to==source) break;
                    from=to;
                }
                outstream << destination_coordinate.first << VcfReader::VCF_SEPARATOR << pos << VcfReader::VCF_SEPARATOR << graph.get_path_name(path) << VcfReader::VCF_SEPARATOR << pos_base << VcfReader::VCF_SEPARATOR << pos_base << new_sequence << VcfReader::VCF_MISSING_CHAR_1 << SUFFIX_BND;
            }
        }
        else { /* NOP: impossible to decide chr and pos. */ }
    });
}


void VariantGraph::get_dangling_nodes(bool mode, unordered_set<nid_t>& out) const {
    if (!mode) out.clear();
    graph.for_each_handle([&](handle_t node) {
        const size_t left_degree = graph.get_degree(node,true);
        const size_t right_degree = graph.get_degree(node,false);
        if (left_degree==0 || right_degree==0) {
            if (mode) out.emplace(graph.get_id(node));
            else cerr << "Dangling node: " << to_string(graph.get_id(node)) << ":" << to_string(graph.get_is_reverse(node)) << "=" << graph.get_sequence(node) << " left_degree: " << to_string(left_degree) << " right_degree: " << to_string(right_degree) << '\n';
        }
    });
}


size_t VariantGraph::get_n_dangling_nodes() const {
    size_t out = 0;
    graph.for_each_handle([&](handle_t node) {
        if (graph.get_degree(node,true)==0 || graph.get_degree(node,false)==0) out++;
    });
    return out;
}


bool VariantGraph::is_dangling_node(const handle_t& node_handle) const {
    return graph.get_degree(node_handle,true)==0 || graph.get_degree(node_handle,false)==0;
}


bool VariantGraph::is_dangling_node(const nid_t& node_id) const {
    if (!graph.has_node(node_id)) return false;
    const handle_t node_handle = graph.get_handle(node_id);
    return graph.get_degree(node_handle,true)==0 || graph.get_degree(node_handle,false)==0;
}


bool VariantGraph::is_flanking_node(const handle_t& node_handle) const {
    if (!is_reference_node(node_handle) || !is_dangling_node(node_handle)) return false;
    return node_to_chromosome.at(graph.get_id(node_handle)).first==main_chromosome;
}


bool VariantGraph::is_flanking_node(const nid_t& node_id) const {
    if (!is_reference_node(node_id) || !is_dangling_node(node_id)) return false;
    return node_to_chromosome.at(node_id).first==main_chromosome;
}


size_t VariantGraph::get_n_interchromosomal_edges() const {
    size_t out = 0;
    graph.for_each_edge([&](edge_t edge) {
        const nid_t node1 = graph.get_id(edge.first);
        const nid_t node2 = graph.get_id(edge.second);
        if (node_to_chromosome.contains(node1) && node_to_chromosome.contains(node2) && node_to_chromosome.at(node1).first!=node_to_chromosome.at(node2).first) out++;
    });
    return out;
}


size_t VariantGraph::get_n_ins_edges() const {
    size_t out = 0;
    graph.for_each_edge([&](edge_t edge) {
        const nid_t node1 = graph.get_id(edge.first);
        const nid_t node2 = graph.get_id(edge.second);
        if (node_to_chromosome.contains(node1)!=node_to_chromosome.contains(node2)) out++;
    });
    return out;
}


bool VariantGraph::is_reference_node(const nid_t& node_id) const {
    return graph.has_node(node_id) && !insertion_handles_set.contains(graph.get_handle(node_id));
}


bool VariantGraph::is_reference_node(const handle_t& node_handle) const {
    return !insertion_handles_set.contains(node_handle) && !insertion_handles_set.contains(graph.flip(node_handle));
}


bool VariantGraph::is_reference_edge(const edge_t& edge) {
    if (!graph.has_edge(edge)) return false;
    const edge_t canonized_edge = graph.edge_handle(edge.first,edge.second);
    return !edge_to_vcf_record.contains(canonized_edge);
}


/**
 * Remark: the implementation is similar to `print_supported_vcf_records()`.
 */
void VariantGraph::get_vcf_records_with_edges_impl(const vector<edge_t>& edges, bool identical, vector<size_t>& out) {
    const size_t N_QUERY_EDGES = edges.size();
    bool all_present, is_increasing, is_decreasing;
    size_t i, j;
    size_t n_edges, rank, first;

    out.clear();

    // Canonizing edges
    tmp_edges.clear();
    for (const auto& edge: edges) tmp_edges.emplace_back(graph.edge_handle(edge.first,edge.second));

    // Initializing `printed` and `flags`.
    for (const auto& edge: tmp_edges) {
        if (!edge_to_vcf_record.contains(edge)) {
            if (identical) return;
            else continue;
        }
        for (const size_t& r: edge_to_vcf_record.at(edge)) initialized.at(r)=false;
    }
    for (const auto& edge: tmp_edges) {
        if (!edge_to_vcf_record.contains(edge)) continue;
        for (const size_t& r: edge_to_vcf_record.at(edge)) {
            if (initialized.at(r)) continue;
            printed.at(r)=false;
            n_edges=vcf_record_to_edge.at(r).size();
            flags.at(r).clear();
            for (i=0; i<n_edges; i++) flags.at(r).emplace_back(0);
            initialized.at(r)=true;
        }
    }

    // Marking `flags`.
    rank=1;
    for (const auto& edge: tmp_edges) {
        if (!edge_to_vcf_record.contains(edge)) { rank++; continue; }
        for (const size_t& r: edge_to_vcf_record.at(edge)) mark_edge(edge,rank,vcf_record_to_edge.at(r),flags.at(r));
        rank++;
    }

    // Finding supported VCF records
    for (const auto& edge: tmp_edges) {
        if (!edge_to_vcf_record.contains(edge)) continue;
        for (const size_t& r: edge_to_vcf_record.at(edge)) {
            if (printed.at(r)) continue;
            n_edges=vcf_record_to_edge.at(r).size();
            first=0;
            for (i=0; i<n_edges; i++) {
                if (vcf_record_to_edge.at(r).at(i)!=null_edge) continue;
                if (!identical || i-first==N_QUERY_EDGES) {
                    all_present=true;
                    for (j=first; j<i; j++) {
                        if (flags.at(r).at(j)==0) { all_present=false; break; }
                    }
                    if (all_present) {
                        is_increasing=true; is_decreasing=true;
                        for (j=first+1; j<i; j++) {
                            if (flags.at(r).at(j)>flags.at(r).at(j-1)) is_decreasing=false;
                            if (flags.at(r).at(j)<flags.at(r).at(j-1)) is_increasing=false;
                        }
                        if (is_increasing || is_decreasing) {
                            out.emplace_back(r);
                            printed.at(r)=true;
                            break;
                        }
                    }
                }
                first=i+1;
            }
        }
    }
    if (out.size()>1) sort(out.begin(),out.end());
    tmp_edges.clear();
}


void VariantGraph::get_vcf_records_with_edges(const vector<edge_t>& edges, vector<VcfRecord>& out) {
    out.clear();
    vector<size_t> records;
    get_vcf_records_with_edges_impl(edges,true,records);
    const size_t length = records.size();
    for (size_t i=0; i<length; i++) out.emplace_back(vcf_records.at(records.at(i)));
}


void VariantGraph::for_each_vcf_record(const function<void(size_t id, const vector<edge_t>& edges_of_the_record, const VcfRecord& record)>& callback) {
    for (size_t i=0; i<n_vcf_records; i++) {
        if (!vcf_record_to_edge.at(i).empty()) callback(i,vcf_record_to_edge.at(i),vcf_records.at(i));
    }
}


void VariantGraph::for_each_vcf_record(const vector<pair<string,bool>>& path, const function<void(size_t id, const vector<edge_t>& edges_of_the_record, const VcfRecord& record)>& callback) {
    const size_t LENGTH = path.size();
    bool is_reverse_from, is_reverse_to;
    size_t i, r, length;
    nid_t node_id_from, node_id_to;
    vector<size_t> ids;
    vector<edge_t> edges;

    for (i=1; i<LENGTH; i++) {
        node_id_from=stoi(path.at(i-1).first); is_reverse_from=path.at(i-1).second;
        node_id_to=stoi(path.at(i).first); is_reverse_to=path.at(i).second;
        edges.emplace_back(graph.edge_handle(graph.get_handle(node_id_from,is_reverse_from),graph.get_handle(node_id_to,is_reverse_to)));
    }
    get_vcf_records_with_edges_impl(edges,false,ids);
    edges.clear();
    length=ids.size();
    for (i=0; i<length; i++) {
        r=ids.at(i);
        callback(r,vcf_record_to_edge.at(r),vcf_records.at(r));
    }
}


path_handle_t VariantGraph::load_gaf_path(const string& path_encoding, const string& path_name, string& buffer) {
    const size_t LENGTH = path_encoding.length();
    bool current_orientation;
    char c;
    nid_t node_id;

    path_handle_t path = graph.create_path_handle(path_name);
    buffer.clear(); current_orientation=true;
    for (size_t i=0; i<LENGTH; i++) {
        c=path_encoding.at(i);
        if (c==GAF_FWD_CHAR || c==GAF_REV_CHAR) {
            node_id=stoi(buffer);
            graph.append_step(path,graph.get_handle(node_id,!current_orientation));
            current_orientation=c==GAF_FWD_CHAR;
            buffer.clear();
        }
        else buffer.push_back(c);
    }
    node_id=stoi(buffer);
    graph.append_step(path,graph.get_handle(node_id,!current_orientation));
    return path;
}


path_handle_t VariantGraph::load_gaf_path(vector<pair<string,bool>>& path, const string& path_name) {
    const size_t LENGTH = path.size();
    nid_t node_id;

    path_handle_t path_handle = graph.create_path_handle(path_name);
    for (size_t i=0; i<LENGTH; i++) {
        node_id=stoi(path.at(i).first);
        graph.append_step(path_handle,graph.get_handle(node_id,path.at(i).second));
    }
    return path_handle;
}


void VariantGraph::canonize_gaf_path(string& input_path, string& output_path, string& buffer) {
    const auto LENGTH = (int32_t)input_path.length();
    char c;
    int32_t i;

    output_path.clear(); buffer.clear();
    for (i=LENGTH-1; i>=0; i--) {
        c=input_path.at(i);
        if (c==GAF_FWD_CHAR || c==GAF_REV_CHAR) {
            output_path.push_back(c==GAF_FWD_CHAR?GAF_REV_CHAR:GAF_FWD_CHAR);
            std::reverse(buffer.begin(),buffer.end());
            output_path.append(buffer);
            buffer.clear();
        }
        else buffer.push_back(c);
    }
    if (output_path>input_path) output_path=input_path;
}


void VariantGraph::print_edge_histograms() const {
    const char SEPARATOR = ',';
    size_t n_edges = edge_to_vcf_record.size();
    if (n_edges==0) return;

    size_t i, last;
    vector<size_t> histogram;

    for (i=0; i<=n_vcf_records; i++) histogram.emplace_back(0);
    for (const auto& [key,value]: edge_to_vcf_record) histogram.at(value.size())++;
    cerr << "Histogram: X = n. VCF records | Y = n. non-reference edges supported by X VCF records\n";
    cerr << "Remark: Only VCF records in input to graph construction are counted.\n";
    for (last=n_vcf_records; last>=1; last--) {
        if (histogram.at(last)>0) break;
    }
    for (i=0; i<=last; i++) cerr << to_string(i) << SEPARATOR << to_string(histogram.at(i)) << '\n';

    histogram.reserve(n_edges+1);
    for (i=0; i<histogram.size(); i++) histogram.at(i)=0;
    for (i=histogram.size(); i<=n_edges; i++) histogram.emplace_back(0);
    for (i=0; i<n_vcf_records; i++) histogram.at(vcf_record_to_edge.at(i).size())++;
    cerr << "Histogram: X = n. non-reference edges | Y = n. VCF records creating X non-reference edges\n";
    cerr << "Remark: Only VCF records in input to graph construction are counted.\n";
    for (last=n_edges-1; last>=1; last--) { if (histogram.at(last)>0) break; }
    for (i=0; i<=last; i++) cerr << to_string(i) << SEPARATOR << to_string(histogram.at(i)) << '\n';

    n_edges=graph.get_edge_count();
    histogram.reserve(n_edges+1);
    for (i=0; i<histogram.size(); i++) histogram.at(i)=0;
    for (i=histogram.size(); i<=n_edges; i++) histogram.emplace_back(0);
    graph.for_each_handle([&](handle_t node) { histogram.at(max(graph.get_degree(node,false),graph.get_degree(node,true)))++; });
    cerr << "Histogram: X = max degree of a node | Y = n. nodes with max degree X\n";
    for (last=n_edges-1; last>=1; last--) { if (histogram.at(last)>0) break; }
    for (i=0; i<=last; i++) cerr << to_string(i) << SEPARATOR << to_string(histogram.at(i)) << '\n';
}


void VariantGraph::print_vcf_records_stats(uint8_t bin_size, const vector<string>& callers) const {
    const char SEPARATOR = ',';
    const size_t N_CALLERS = callers.size();
    size_t i, j;
    vector<size_t> lengths_ins, lengths_dup, lengths_del, lengths_inv;
    vector<vector<size_t>> caller_count;
    vector<pair<string,string>> bnds;
    string buffer;

    for (i=0; i<N_CALLERS; i++) caller_count.emplace_back(vector<size_t> {0,0,0,0,0,});
    for (const auto& record: vcf_records) {
        if (record.sv_type==VcfReader::TYPE_INSERTION) {
            lengths_ins.emplace_back(record.sv_length/bin_size);
            if (!callers.empty()) increment_caller_count(record,callers,caller_count,0);
        }
        else if (record.sv_type==VcfReader::TYPE_DUPLICATION) {
            lengths_dup.emplace_back(record.sv_length/bin_size);
            if (!callers.empty()) increment_caller_count(record,callers,caller_count,1);
        }
        else if (record.sv_type==VcfReader::TYPE_DELETION) {
            lengths_del.emplace_back(record.sv_length/bin_size);
            if (!callers.empty()) increment_caller_count(record,callers,caller_count,2);
        }
        else if (record.sv_type==VcfReader::TYPE_INVERSION) {
            lengths_inv.emplace_back(record.sv_length/bin_size);
            if (!callers.empty()) increment_caller_count(record,callers,caller_count,3);
        }
        else if (record.sv_type==VcfReader::TYPE_BREAKEND) {
            record.get_breakend_chromosome(buffer);
            if (record.chrom<buffer) bnds.emplace_back(record.chrom,buffer);
            else bnds.emplace_back(buffer,record.chrom);
            if (!callers.empty()) increment_caller_count(record,callers,caller_count,4);
        }
    }
    cerr << "Histogram: SV type | SV length (binned by " << to_string(bin_size) << "bp) | Number of occurrences\n";
    print_sv_lengths(lengths_ins,bin_size,VcfReader::INS_STR);
    print_sv_lengths(lengths_dup,bin_size,VcfReader::DUP_STR);
    print_sv_lengths(lengths_del,bin_size,VcfReader::DEL_STR);
    print_sv_lengths(lengths_inv,bin_size,VcfReader::INV_STR);
    cerr << "Histogram: BND | chrA | chrB | Number of occurrences\n";
    print_bnds(bnds);
    if (!callers.empty()) {
        cerr << "Histogram: Caller | #INS | #DUP | #DEL | #INV | #BND\n";
        for (i=0; i<caller_count.size(); i++) {
            cerr << callers.at(i);
            for (j=0; j<caller_count.at(i).size(); j++) cerr << SEPARATOR << to_string(caller_count.at(i).at(j));
            cerr << '\n';
        }
    }
}


void VariantGraph::print_sv_lengths(vector<size_t>& lengths, uint8_t bin_size, const string& prefix) {
    const char SEPARATOR = ',';
    size_t previous, previous_count;
    const size_t SIZE = lengths.size();

    if (SIZE==0) return;
    if (SIZE>1) sort(lengths.begin(),lengths.end());
    previous=lengths.at(0); previous_count=1;
    for (size_t i=1; i<SIZE; i++) {
        if (lengths.at(i)!=previous) {
            cerr << prefix << SEPARATOR << to_string(previous*bin_size) << SEPARATOR << to_string(previous_count) << '\n';
            previous=lengths.at(i); previous_count=1;
        }
        else previous_count++;
    }
    cerr << prefix << SEPARATOR << to_string(previous*bin_size) << SEPARATOR << to_string(previous_count) << '\n';
}


void VariantGraph::print_bnds(vector<pair<string,string>>& bnds) {
    const char SEPARATOR = ',';
    const size_t SIZE = bnds.size();
    size_t previous_count;

    if (SIZE==0) return;
    if (SIZE>1) sort(bnds.begin(),bnds.end());
    pair<string,string>& previous = bnds.at(0);
    previous_count=1;
    for (size_t i=1; i<SIZE; i++) {
        if (bnds.at(i)!=previous) {
            cerr << "BND" << SEPARATOR << previous.first << SEPARATOR << previous.second << SEPARATOR << to_string(previous_count) << '\n';
            previous=bnds.at(i); previous_count=1;
        }
        else previous_count++;
    }
    cerr << "BND" << SEPARATOR << previous.first << SEPARATOR << previous.second << SEPARATOR << to_string(previous_count) << '\n';
}


void VariantGraph::increment_caller_count(const VcfRecord& record, const vector<string>& callers, vector<vector<size_t>>& caller_count, uint8_t column) {
    const size_t N_CALLERS = callers.size();
    string lower = record.id;
    lowercase_string(lower);
    for (size_t i=0; i<N_CALLERS; i++) {
        if (lower.find(callers.at(i))!=string::npos) { caller_count.at(i).at(column)++; return; }
    }
}


void VariantGraph::print_graph_signature(size_t max_steps, const path& output_path) const {
    const char SEPARATOR = ',';
    size_t i;
    string buffer;
    vector<string> signatures, tokens;

    graph.for_each_handle([&](handle_t handle) {
        tokens.clear();
        print_signature_impl(handle,graph.get_sequence(handle),0,max_steps,tokens);
        const handle_t rev_handle = graph.flip(handle);
        print_signature_impl(rev_handle,graph.get_sequence(rev_handle),0,max_steps,tokens);
        sort(tokens.begin(),tokens.end());
        buffer.clear();
        for (i=0; i<tokens.size(); i++) { buffer.append(tokens.at(i)); buffer.push_back(SEPARATOR); }
        signatures.emplace_back(buffer);
    });
    sort(signatures.begin(),signatures.end());
    ofstream output_file(output_path.string());
    if (!output_file.good() || !output_file.is_open()) throw runtime_error("ERROR: file could not be written: " + output_path.string());
    for (i=0; i<signatures.size(); i++) output_file << signatures.at(i) << '\n';
    output_file.close();
}


void VariantGraph::print_signature_impl(const handle_t& handle, const string& path, size_t steps_performed, size_t max_steps, vector<string>& strings) const {
    if (steps_performed==max_steps || graph.get_degree(handle,false)==0) {
        strings.emplace_back(path);
        return;
    }
    graph.follow_edges(handle,false,[&](handle_t neighbor_handle) {
        string new_path = path;
        new_path.append(graph.get_sequence(neighbor_handle));
        print_signature_impl(neighbor_handle,new_path,steps_performed+1,max_steps,strings);
    });
}


void VariantGraph::load_edge_record_map(const vector<pair<edge_t,size_t>>& map, size_t n_vcf_records) {
    size_t i;

    edge_to_vcf_record.clear(); vcf_record_to_edge.clear();
    vcf_record_to_edge.reserve(n_vcf_records);
    for (i=0; i<n_vcf_records; i++) vcf_record_to_edge.emplace_back();
    for (const auto& pair: map) {
        const edge_t canonized_edge = graph.edge_handle(pair.first.first,pair.first.second);
        if (edge_to_vcf_record.contains(canonized_edge)) edge_to_vcf_record.at(canonized_edge).emplace_back(pair.second);
        else edge_to_vcf_record[canonized_edge]={pair.second};
        vcf_record_to_edge.at(pair.second).emplace_back(canonized_edge);
    }
    for (i=0; i<n_vcf_records; i++) vcf_record_to_edge.at(i).emplace_back(null_edge);
}


unordered_map<edge_t,vector<size_t>> VariantGraph::get_edge_record_map() {
    unordered_map<edge_t,vector<size_t>> out;
    out=edge_to_vcf_record;
    return out;
}


void VariantGraph::node_to_chromosome_clear() { node_to_chromosome.clear(); }


void VariantGraph::node_to_chromosome_insert(const string& node_label, const vector<string>& node_labels, const string& chromosome, int32_t position) {
    const auto iterator = lower_bound(node_labels.begin(),node_labels.end(),node_label);
    const nid_t node_id = distance(node_labels.begin(),iterator)+1;  // Node IDs cannot be zero
    node_to_chromosome[node_id]=pair<string,int32_t>(chromosome,position);
}


void VariantGraph::insertion_handles_set_clear() { insertion_handles_set.clear(); }


void VariantGraph::insertion_handles_set_insert(const string& node_label, const vector<string>& node_labels) {
    const auto iterator = lower_bound(node_labels.begin(),node_labels.end(),node_label);
    const nid_t node_id = distance(node_labels.begin(),iterator)+1;  // Node IDs cannot be zero
    insertion_handles_set.insert(graph.get_handle(node_id));
}

}