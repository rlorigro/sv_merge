#include "path_optimizer_mathopt.hpp"
#include "SimpleAlignment.hpp"
#include "TransitiveMap.hpp"
#include "interval_tree.hpp"
#include "VariantGraph.hpp"
#include "debug_mode.hpp"
#include "VcfReader.hpp"
#include "fasta.hpp"
#include "Timer.hpp"
#include "CLI11.hpp"
#include "fetch.hpp"
#include "misc.hpp"
#include "gaf.hpp"
#include "bed.hpp"

using lib_interval_tree::interval_tree_t;

#include "bdsg/hash_graph.hpp"

using bdsg::HandleGraph;

#include <filesystem>

using std::filesystem::path;
using std::filesystem::create_directories;


#include <unordered_map>
#include <exception>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <thread>
#include <limits>
#include <cmath>
#include <tuple>

using std::this_thread::sleep_for;
using std::numeric_limits;
using std::unordered_map;
using std::runtime_error;
using std::exception;
using std::ifstream;
using std::ofstream;
using std::thread;
using std::atomic;
using std::mutex;
using std::log10;
using std::tuple;
using std::cerr;
using std::max;
using std::min;
using std::cref;
using std::ref;

#include "cpptrace/from_current.hpp"

using namespace sv_merge;

#include "bindings/cpp/WFAligner.hpp"
using namespace wfa;


void cross_align_sample_reads(TransMap& transmap, int64_t score_threshold, const string& sample_name){
    WFAlignerGapAffine aligner(4,6,2,WFAligner::Alignment,WFAligner::MemoryHigh);

    vector <tuple <int64_t,int64_t,float> > edges_to_add;
    string seq_a;
    string seq_b;

    auto sample_id = transmap.get_id(sample_name);
    transmap.for_each_read_of_sample(sample_id, [&](int64_t a){
        transmap.for_each_read_of_sample(sample_id, [&](int64_t b){
            if (transmap.has_edge(a,b)){
                return;
            }

            transmap.get_sequence(a,seq_a);
            transmap.get_sequence(b,seq_b);

            if (llabs(int64_t(seq_a.size()) - int64_t(seq_b.size())) > score_threshold){
//                cerr << "skipping edge: " << a << ',' << b << " l: " << int64_t(seq_a.size()) << ',' << int64_t(seq_b.size()) << '\n';
                return;
            }

            if (a != b) {
                aligner.alignEnd2End(seq_a, seq_b);
                auto score = aligner.getAlignmentScore();

                // WFA creates negative scores for "distance", we make it positive again and rescale it as a percent
                size_t scaled_score = 100 - (100*size_t(-score)) / min(seq_a.size(), seq_b.size());

                // Avoid adding while iterating
                edges_to_add.emplace_back(a,b,scaled_score);

//                cerr << transmap.get_node(a).name << ',' << transmap.get_node(b).name << ' ' << int64_t(seq_a.size()) << ',' << int64_t(seq_b.size()) << ' ' << score << ' ' << scaled_score << '\n';
            }
        });
    });

    for (const auto& [a,b,score]: edges_to_add){
        transmap.add_edge(a,b,score);
    }
}


void write_alignment_debug_info(
        const VariantGraph& variant_graph,
        const string& sample_name,
        const string& read_name,
        const string& path_name,
        const string& path_sequence,
        const string& read_sequence,
        const string& cigar,
        const interval_t& path_evaluation_window,
        float nonmatch_portion,
        path output_dir
        ){

    vector <pair <string,bool> > p;
    GafAlignment::parse_string_as_path(path_name, p);

    vector<interval_t> ref_intervals;
    vector<interval_t> empty_intervals;

    // Get the corresponding bdsg handles from the handlegraph for each step in the path and concatenate their seqs
    for (const auto& [node_name, reversal]: p){
        nid_t n = stoll(node_name);
        auto h = variant_graph.graph.get_handle(n, reversal);
        auto l = int32_t(variant_graph.graph.get_length(h));

        if (ref_intervals.empty()){
            ref_intervals.push_back({0,l});
        }
        else{
            ref_intervals.push_back({ref_intervals.back().second, ref_intervals.back().second + l});
        }
    }

    path out_path = output_dir / (sample_name + ".txt");

    if (not exists(output_dir)) {
        create_directories(output_dir);
    }

    ofstream file(out_path, ofstream::app);

    string r;
    string r_all;
    string x;
    string x_all;
    string q;
    string q_all;
    string bounds;

    SimpleAlignment alignment(path_sequence, read_sequence, cigar);

    auto [a,b] = path_evaluation_window;

    file << read_name << " " << path_name << ' ' << '\n';
    file << cigar << '\n';

    interval_t prev_interval = {-1,-1};
    for_cigar_interval_in_alignment(false, alignment, ref_intervals, empty_intervals,
        [&](const CigarInterval &intersection, const interval_t &interval) {
            // use the cigars to construct a formatted alignment
            get_formatted_sequence_of_cigar_interval(
                    intersection,
                    path_sequence,
                    read_sequence,
                    r,
                    q,
                    x
            );

            file << cigar_code_to_char[intersection.code] << ' ' << intersection.get_op_length() << '\n';

            r_all += r;
            x_all += x;
            q_all += q;

            if (prev_interval.first != -1){
                if (interval != prev_interval){
                    bounds.back() = '|';
                }
            }

            bounds += string(x.size(), ' ');

            prev_interval = interval;
        },
        [&](const CigarInterval &intersection, const interval_t &interval) {
            return;
        }
    );

    file << "path length: " << path_sequence.size() << '\n';
    file << "read length: " << read_sequence.size() << '\n';
    file << "nonmatch portion: " << nonmatch_portion << '\n';
    file << bounds << '\n';
    file << r_all << '\n';
    file << x_all << '\n';
    file << q_all << '\n';
    file << '\n';
}


void align_read_to_path(
        const string& read_sequence,
        const string& path_sequence,
        WFAlignerGapAffine& aligner,
        interval_t path_evaluation_window,
        string& cigar_result,
        float& non_match_portion_result
        ){

    // TODO: Add length cutoff before aligning
    aligner.alignEnd2End(path_sequence, read_sequence);

    // Extract indel edit distance from cigar
    cigar_result.clear();
    cigar_result = aligner.getCIGAR(true);
    string length_token;

    int64_t n_non_match = 0;
    int64_t n = 0;

    SimpleAlignment alignment(path_sequence, read_sequence, cigar_result);

    auto [a,b] = path_evaluation_window;
    vector<interval_t> ref_intervals = {{a,b}};
    vector<interval_t> empty_intervals;

    // Accumulate indels and total length for non-flanking path intervals.
    // We don't want to keep track of indels in the flanking regions. However, sometimes the aligner squishes
    // them into the flanking regions, so we have a small buffer inside the flanks to catch these.
    for_cigar_interval_in_alignment(false, alignment, ref_intervals, empty_intervals,
        [&](const CigarInterval& intersection, const interval_t& interval) {
            auto l = intersection.get_op_length();

            // Count non-matches
            if (intersection.code != cigar_char_to_code['=']){
                n_non_match += l;
            }

            n += l;

        },

        // Do nothing for query intervals
        [&](const CigarInterval& intersection, const interval_t& interval){return;}
    );

    non_match_portion_result = float(n_non_match) / float(n);
}


/**
 * Method to align ALL reads to ALL paths
 * @param transmap will be filled with new edges, reads with no edge are deleted !!!
 * @param variant_graph used to infer path sequences for all paths
 * @param min_read_hap_identity threshold to accept an edge
 * @param flank_length used to filter the region that cigars are computed
 * @param output_dir for debug/logging
 */
void align_reads_vs_paths(TransMap& transmap, const VariantGraph& variant_graph, float min_read_hap_identity, int32_t flank_length, path output_dir=""){
    // TODO: test these params??
    WFAlignerGapAffine aligner(4,6,2,WFAligner::Alignment,WFAligner::MemoryHigh);

    // First extract the sequences of the paths
    unordered_map<string,string> path_sequences;
    unordered_map<string,pair<int32_t,int32_t> > path_flank_coords;

    // Find the windows that will be evaluated (exclude most of the flanks)
    transmap.for_each_path([&](const string& name, int64_t id){
        vector <pair <string,bool> > path;
        GafAlignment::parse_string_as_path(name, path);

        string& sequence = path_sequences[name];

        int32_t total_length = 0;
        int32_t buffer_size = min(60, flank_length - 1);

        path_flank_coords[name] = {-1,-1};

        // Get the corresponding bdsg handles from the handlegraph for each step in the path and concatenate their seqs
        for (const auto& [node_name, reversal]: path){
            nid_t n = stoll(node_name);
            auto h = variant_graph.graph.get_handle(n, reversal);
            auto l = int32_t(variant_graph.graph.get_length(h));

            if (node_name == path.front().first and variant_graph.is_dangling_node(n)){
                path_flank_coords[name].first = min(l, flank_length - buffer_size);
            }

            if (node_name == path.back().first and variant_graph.is_dangling_node(n)){
                path_flank_coords[name].second = total_length + l - min(l, flank_length - buffer_size);
            }

            sequence += variant_graph.graph.get_sequence(h);

            total_length += l;
        }

        if (path_flank_coords[name].first == -1){
            path_flank_coords[name].first = 0;
        }
        if (path_flank_coords[name].second == -1){
            path_flank_coords[name].second = sequence.size();
        }
    });

    // Force all path sequences to uppercase
    // TODO: remove when VariantGraph supports uppercasing
    for (auto& [name, seq]: path_sequences){
        for (size_t i=0; i<seq.size(); i++){
            seq[i] = toupper(seq[i]);
        }
    }

    unordered_map<string, unordered_set<int64_t> > deduplicated_reads;

    vector <tuple <int64_t,int64_t,float> > edges_to_add;
    vector <int64_t> reads_to_be_removed;
    vector<interval_t> empty_intervals;
    float non_match_portion;
    string cigar;
    string read_sequence;

    // Deduplicate by hashing strings
    transmap.for_each_read([&](const string& read_name, int64_t id) {
        transmap.get_sequence(id, read_sequence);
        deduplicated_reads[read_sequence].insert(id);
    });

    for (const auto& [read_sequence, read_ids]: deduplicated_reads){
        bool has_alignment = false;

        for (const auto& [path_name, path_sequence]: path_sequences) {
            if (read_sequence.empty() or path_sequence.empty()){
                continue;
            }

            auto path_evaluation_window = path_flank_coords[path_name];

            align_read_to_path(read_sequence, path_sequence, aligner, path_evaluation_window, cigar, non_match_portion);

            if ((1.0f - non_match_portion) > min_read_hap_identity) {
                // Store the permil score as a rounded int, add 0.001 (1 permil) to avoid NaNs in the cost function
                auto int_score = int64_t(round(non_match_portion*1000)) + 1;

                // For all reads that had this identical sequence, update the transmap edges
                for (auto id: read_ids) {
                    // Avoid adding while iterating
                    edges_to_add.emplace_back(id, transmap.get_id(path_name), int_score);
                    has_alignment = true;
                }
            }

            // Write out the full alignments for debugging
            if (HAPESTRY_DEBUG) {
                for (auto id: read_ids) {
                    auto read_name = transmap.get_node(id).name;

                    string sample_name = transmap.get_sample_of_read(read_name);
                    write_alignment_debug_info(
                            variant_graph,
                            sample_name,
                            read_name,
                            path_name,
                            path_sequence,
                            read_sequence,
                            cigar,
                            path_evaluation_window,
                            non_match_portion,
                            output_dir
                        );
                }
            }
        }

        if (not has_alignment) {
            for (auto id: read_ids){
                reads_to_be_removed.emplace_back(id);;
            }
        }
    }

    for (auto [a,b,score]: edges_to_add){
        transmap.add_edge(a,b,score);
    }

    // Doing this because later on in hapestry it becomes annoying if some reads have no edges
    for (auto id: reads_to_be_removed) {
        transmap.remove_node(id);
    }
}


/**
 * Only used for debugging purposes
 * @param transmap
 * @param variant_graph
 * @param min_read_hap_identity
 * @param flank_length
 * @param output_dir
 */
void align_read_path_edges_of_transmap(TransMap& transmap, const VariantGraph& variant_graph, float min_read_hap_identity, int32_t flank_length, path output_dir=""){
    // TODO: test these params??
    WFAlignerGapAffine aligner(4,6,2,WFAligner::Alignment,WFAligner::MemoryHigh);

    // First extract the sequences of the paths
    unordered_map<string,string> path_sequences;
    unordered_map<string,pair<int32_t,int32_t> > path_flank_coords;

    transmap.for_each_path([&](const string& name, int64_t id){
        vector <pair <string,bool> > path;
        GafAlignment::parse_string_as_path(name, path);

        string& sequence = path_sequences[name];

        int32_t total_length = 0;
        int32_t buffer_size = 30;

        path_flank_coords[name] = {-1,-1};

        // Get the corresponding bdsg handles from the handlegraph for each step in the path and concatenate their seqs
        for (const auto& [node_name, reversal]: path){
            nid_t n = stoll(node_name);
            auto h = variant_graph.graph.get_handle(n, reversal);
            auto l = int32_t(variant_graph.graph.get_length(h));

            if (node_name == path.front().first and variant_graph.is_dangling_node(n)){
                path_flank_coords[name].first = min(l, flank_length - buffer_size);
            }

            if (node_name == path.back().first and variant_graph.is_dangling_node(n)){
                path_flank_coords[name].second = total_length + l - min(l, flank_length - buffer_size);
            }

            sequence += variant_graph.graph.get_sequence(h);

            total_length += l;
        }

        if (path_flank_coords[name].first == -1){
            path_flank_coords[name].first = 0;
        }
        if (path_flank_coords[name].second == -1){
            path_flank_coords[name].second = total_length;
        }
    });

    vector <tuple <int64_t,int64_t,float> > edges_to_add;
    vector<interval_t> empty_intervals;
    float nonmatch_portion;
    string read_sequence;
    string cigar;

    transmap.for_each_read([&](const string& read_name, int64_t id){
        transmap.get_sequence(id, read_sequence);

        transmap.for_each_path_of_read(id, [&](int64_t path_id){
            string path_name = transmap.get_node(path_id).name;
            auto& path_sequence = path_sequences[path_name];

            if (read_sequence.empty() or path_sequence.empty()){
                return;
            }

            auto path_evaluation_window = path_flank_coords[path_name];

            align_read_to_path(read_sequence, path_sequence, aligner, path_evaluation_window, cigar, nonmatch_portion);

            if ((1.0f - nonmatch_portion) > min_read_hap_identity) {
                // Store the permil score as a rounded int, add 0.001 (1 permil) to avoid NaNs in the cost function
                auto int_score = int64_t(round(nonmatch_portion*1000)) + 1;

                // Avoid adding while iterating
                edges_to_add.emplace_back(id, transmap.get_id(path_name), int_score);
            }

            // Write out the full alignments for debugging
            if (HAPESTRY_DEBUG) {
                string sample_name = transmap.get_sample_of_read(read_name);
                write_alignment_debug_info(
                        variant_graph,
                        sample_name,
                        read_name,
                        path_name,
                        path_sequence,
                        read_sequence,
                        cigar,
                        path_evaluation_window,
                        nonmatch_portion,
                        output_dir
                    );
            }
        });
    });

    for (const auto& [a,b,score]: edges_to_add){
        transmap.add_edge(a,b,score);
    }
}


void write_region_subsequences_to_file(const TransMap& t, const path& output_fasta) {
    string s;

    ofstream file(output_fasta);

    t.for_each_read([&](const string& name, int64_t id){
        t.get_sequence(id,s);

        if (s.empty()){
            return;
        }

        file << '>' << name << '\n';
        file << s << '\n';
    });
}


void write_deduplicated_region_subsequences_to_file(const TransMap& t, const path& output_fasta) {
    unordered_map<string,int64_t> counter;
    string s;

    // Decompress and count the reads
    t.for_each_read([&](const string& name, int64_t id){
        t.get_sequence(id,s);

        if (s.empty()){
            return;
        }

        counter[s]++;
    });

    ofstream file(output_fasta);

    // Write arbitrary read names to the FASTA, and keep track of their count (in case needed for coverage later)
    size_t i = 0;
    for (const auto& [seq, n]: counter) {
        file << '>' << i << '_' << n << '\n';
        file << seq << '\n';

        i++;
    }
}


void write_region_subsequences_to_file_thread_fn(
        const unordered_map<Region,TransMap>& region_transmaps,
        const vector<Region>& regions,
        const path& output_dir,
        const path& filename,
        const int32_t flank_length,
        bool deduplicate_reads,
        atomic<size_t>& job_index
){
    size_t i = job_index.fetch_add(1);

    while (i < regions.size()){
        const auto& region = regions.at(i);
        const auto& t = region_transmaps.at(region);

        path output_subdir = output_dir / region.to_unflanked_string('_', flank_length);

        create_directories(output_subdir);

        path output_fasta = output_subdir / filename;

        if (deduplicate_reads) {
            write_deduplicated_region_subsequences_to_file(t, output_fasta);
        }
        else {
            write_region_subsequences_to_file(t, output_fasta);
        }

        i = job_index.fetch_add(1);
    }
}


void get_path_coverages(path gaf_path, const VariantGraph& variant_graph, unordered_map<string,int64_t>& path_coverage){
    // Despite removing duplicate variants, it is still possible for a duplicate path sequence to exist with a different
    // combination of variants. In this case, we reassign any duplicates arbitrarily to the first path iterated

    // This object stores path_sequence -> path_name, and is used to reassign paths to the first path with the same
    // sequence
    unordered_map<string,string> path_reassignment;

    // Accumulate coverage of spanning paths in the GAF
    for_alignment_in_gaf(gaf_path, [&](GafAlignment& alignment){
        auto& path = alignment.get_path();
        if (path.size() < 2){
            return;
        }

        nid_t id_front = stoll(path.front().first);
        nid_t id_back = stoll(path.back().first);

        // Check if spanning
        auto l = variant_graph.is_dangling_node(id_front);
        auto r = variant_graph.is_dangling_node(id_back);
        bool is_spanning = l and r;

//        cerr << alignment.get_path_string() << ' ' << is_spanning << '\n';

        if (not is_spanning){
            return;
        }

        auto path_name = alignment.get_path_string();

        string sequence;

        // Get the corresponding bdsg handles from the handlegraph for each step in the path and concatenate their seqs
        for (const auto& [node_name, reversal]: path) {
            nid_t n = stoll(node_name);
            auto h = variant_graph.graph.get_handle(n, reversal);
            sequence += variant_graph.graph.get_sequence(h);
        }

        // Check if this path has been seen before
        auto result = path_reassignment.find(sequence);
        bool is_new = (result == path_reassignment.end());

        if (is_new){
            // Add the path to the reassignment map and increment its coverage
            path_reassignment.emplace(sequence, alignment.get_path_string());
            path_coverage[path_name]++;
        }
        else{
            // Arbitrarily reassign the path to the first path with the same sequence
            path_coverage[result->second]++;
        }
    });
}


void write_solution_as_hap_vcf(
        VariantGraph& variant_graph,
        const TransMap& transmap,
        const vector<string>& sample_names,
        const unordered_map<string,string>& ref_sequences,
        const path& output_path,
        const Region& region,
        int32_t flank_length
        ){

    // Open the VCF output file
    ofstream vcf_file(output_path);

    if (not (vcf_file.is_open() and vcf_file.good())){
        throw runtime_error("ERROR: could not write file: " + output_path.string());
    }

    // Write the VCF header
    vcf_file << "##fileformat=VCFv4.2" << '\n';
    vcf_file << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (const auto& sample_name: sample_names){
        vcf_file << '\t' << sample_name;
    }
    vcf_file << '\n';

    size_t unflanked_length = region.stop-region.start-(2*flank_length);
    size_t unflanked_start = region.start + flank_length;

    const string& ref_sequence = ref_sequences.at(region.name);

    // Since there may be NO ref sequence between the flank nodes, we must anchor on the last base of the left flank
    const char anchor_base = ref_sequence[unflanked_start - 1];

    const string ref_allele = anchor_base + ref_sequence.substr(unflanked_start, unflanked_length);

    // cerr << "ref_allele: " << ref_allele << '\n';

    unordered_map<string,string> path_sequences;

    // Some of the paths may include BND operations that do not terminate on the flanks of this region, so here we opt
    // to just skip them
    unordered_set<string> invalid_paths;

    // There should be only one path that is exactly the reference path
    string ref_path_name;

    unordered_set<string> samples_not_found;
    unordered_map<string, unordered_map<int64_t, array<bool,2> > > sample_gts;

    // Get an ordered representation of ALL the paths that each sample has. Use this to construct consistent GTs
    transmap.for_each_sample([&](const string& sample_name, int64_t sample_id) {
        set<int64_t> path_ids;

        // First find an ordered assignment of phases, which is sorted by path id
        transmap.for_each_path_of_sample(sample_name, [&](const string& p_name, int64_t p){
            path_ids.emplace(p);
        });

        if (path_ids.size() > 2){
            string s;
            for (auto p: path_ids){
                s += transmap.get_node(p).name + " ";
            }
            throw runtime_error("ERROR: more than two paths found for sample " + sample_name + ": " + s);
        }

        // If a path has
        if (path_ids.empty()) {
            return;
        }

        for (const auto& p: path_ids) {
            sample_gts[sample_name][p][0] = *path_ids.begin() == p;
            sample_gts[sample_name][p][1] = *path_ids.rbegin() == p;
        }
    });

    transmap.for_each_path([&](const string& path_name, int64_t path_id) {
        string var_id = region.name + "_" + to_string(unflanked_start) + "_" + to_string(path_id);

        string sequence;

        vector <pair <string,bool> > path;
        GafAlignment::parse_string_as_path(path_name, path);

        bool is_only_ref = true;

        // cerr << "path_name: " << path_name << '\n';

        bool has_left_flank = false;
        bool has_right_flank = false;

        // Get the corresponding bdsg handles from the handlegraph for each step in the path and concatenate their seqs
        for (size_t i=0; i<path.size(); i++){
            auto n = stoll(path[i].first);
            auto r = path[i].second;
            auto h = variant_graph.graph.get_handle(n, r);

            bool is_ref_node = variant_graph.is_reference_node(n);

            // cerr << n << ' ' << sequence << '\n';

            if (i > 0) {
                auto n_prev = stoll(path[i-1].first);
                auto r_prev = path[i-1].second;
                auto h_prev = variant_graph.graph.get_handle(n_prev, r_prev);

                auto e = variant_graph.graph.edge_handle(h_prev, h);

                if (not variant_graph.is_reference_edge(e)) {
                    // cerr << "edge is NOT ref" << '\n';
                    is_only_ref = false;
                }
            }

            // If a node is encountered that is a tip, and also is not reference, then we are likely in a BND
            if (variant_graph.is_dangling_node(n)){
                if (is_ref_node) {
                    if (i == 0) {
                        has_left_flank = true;
                    }
                    else if (i == path.size() - 1) {
                        has_right_flank = true;
                    }
                }
                else {
                    invalid_paths.emplace(path_name);
                    // Break callback function, do not write to hap VCF
                    return;
                }
            }

            sequence += variant_graph.graph.get_sequence(h);
        }

        if (is_only_ref) {
            if (ref_path_name.empty()) {
                ref_path_name = path_name;
            }
            else if (ref_path_name != path_name) {
                cerr << "WARNING: multiple paths are ref-only nodes: " + ref_path_name + " != " + path_name << '\n';
            }

            // Break callback function, do not write to hap VCF because you cannot substitute a REF allele with itself
            return;
        }

        // TODO: infer an actual quality from the read-hap edge weights...

        if (has_left_flank) {
            // trim left, but keep last base of the flank as anchor base
            sequence = sequence.substr(flank_length-1);
        }

        if (has_right_flank) {
            // trim right
            sequence = sequence.substr(0,sequence.size()-flank_length);
        }

        // #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	[SAMPLES]
        vcf_file << region.name << '\t';
        vcf_file << unflanked_start - 1 << '\t';
        vcf_file << var_id << '\t';
        vcf_file << ref_allele << '\t';
        vcf_file << sequence << '\t';
        vcf_file << '.' << '\t';
        vcf_file << '.' << '\t';
        vcf_file << '.' << '\t';
        vcf_file << "GT" << '\t';

        for (size_t i=0; i<sample_names.size(); i++) {
            const auto& sample_name = sample_names[i];
            auto result = sample_gts.find(sample_name);

            if (result == sample_gts.end()) {
                samples_not_found.emplace(sample_name);
                vcf_file << "0|0" << ((i < sample_names.size() - 1) ? "\t" : "");
                continue;
            }

            const auto& path_gts = result->second;
            auto result2 = path_gts.find(path_id);

            if (result2 == path_gts.end()) {
                vcf_file << "0|0" << ((i < sample_names.size() - 1) ? "\t" : "");
                continue;
            }

            const auto& gt = result2->second;
            vcf_file << gt[0] << '|' << gt[1] << ((i < sample_names.size() - 1) ? "\t" : "");
        }

        vcf_file << '\n';
    });

    if (not samples_not_found.empty()) {
        string s;

        for (const auto& x: samples_not_found) {
            s += x + " ";
        }
        cerr << "WARNING: for hap VCF region " << region.to_string() << " transmap is missing VCF samples " << s << '\n';
    }
}


void write_solution_to_vcf(
        VariantGraph& variant_graph,
        const TransMap& transmap,
        const vector<string>& sample_names,
        const path& output_path,
        const Region& region
        ){

    // TODO: consider just directly overwriting the vector<string> in the VcfRecords stored by VariantGraph
    // Generate a vector of genotypes for each sample, where the vector indexes correspond to the VariantGraph indexes
    unordered_map <string, vector <array<int8_t,2> > > sample_genotypes;

    cerr << output_path << '\n';

    string samples_not_found;

    // TODO: find a way to deal with regions for which the transmap was cleared by the optimizer and no calls were made
    // (for now they get 0|0 genotypes)
    for (const auto& sample_name: sample_names){
        // Initialize the vectors with arrays of {0,0}
        sample_genotypes[sample_name] = vector <array<int8_t,2> >(variant_graph.vcf_records.size(), {0,0});

        // Get the sample ID
        auto [success,sample_id] = transmap.try_get_id(sample_name);

        // Maybe some info was missing for a given sample, in which case give it 0/0 gt
        if (not success){
            samples_not_found += " " + sample_name;
            continue;
        }

        transmap.for_each_phased_variant_of_sample(sample_id, [&](const string& var_name, int64_t _, bool phase){
            // Reconstruct the variant ID from the name
            size_t v = stoull(var_name.substr(1));

            // If there is a transitive edge from sample->path->variant, set the genotype to 1 (with phase)
            sample_genotypes[sample_name][v][phase] = 1;
        });
    }

    if (not samples_not_found.empty()) {
        cerr << "WARNING: region " << region.to_string() << " transmap is missing VCF samples" << samples_not_found << '\n';
    }

    // Open the VCF output file
    ofstream vcf_file(output_path);

    if (not (vcf_file.is_open() and vcf_file.good())){
        throw runtime_error("ERROR: could not write file: " + output_path.string());
    }

    // Write the VCF header
    vcf_file << "##fileformat=VCFv4.2" << '\n';
    vcf_file << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (const auto& sample_name: sample_names){
        vcf_file << '\t' << sample_name;
    }
    vcf_file << '\n';

    // TODO: infer an actual quality from the read-hap edge weights...

    for (size_t v=0; v<variant_graph.vcf_records.size(); v++){
        auto& record = variant_graph.vcf_records.at(v);

        record.genotypes.clear();
        record.format = "GT";

        bool has_nonzero_gt = false;

        // Iterate all the samples and accumulate GTs for the variant object
        for (const auto& sample_name: sample_names){
            string variant_name = "v" + to_string(v);

            // Get the genotype for this sample
            auto& genotype = sample_genotypes.at(sample_name).at(v);

            // Update the record genotypes
            record.genotypes.emplace_back(to_string(int(genotype[0])) + "|" + to_string(int(genotype[1])));

            if (genotype[0] != 0 or genotype[1] != 0){
                has_nonzero_gt = true;
            }
        }

        // Write the record to the VCF
        if (has_nonzero_gt){
            record.print(vcf_file);
            vcf_file << '\n';
        }
    }
}


void write_sample_path_divergence(VariantGraph& variant_graph, const TransMap& transmap, path output_dir){
    unordered_map <int64_t, unordered_map <int64_t, vector<float> > > sample_path_divergence;

    transmap.for_each_sample([&](const string& sample_name, int64_t sample_id){
        transmap.for_each_read_of_sample(sample_id, [&](int64_t read_id){
            transmap.for_each_path_of_read(read_id, [&](int64_t path_id){
                auto [_, score] = transmap.try_get_edge_weight(path_id,read_id);
                sample_path_divergence[sample_id][path_id].push_back(score);
            });
        });
    });

    path divergence_path = output_dir / "sample_path_divergence.csv";
    ofstream divergence_file(divergence_path);

    divergence_file << "sample,path,avg_score,variant_ids\n";
    for (const auto& [sample_id, path_scores]: sample_path_divergence){
        for (const auto& [path_id, scores]: path_scores){
            float avg = 0;
            for (const auto& score: scores){
                avg += score;
            }
            avg /= scores.size();

            string var_ids;

            transmap.for_each_variant_of_path(path_id, [&](int64_t variant_id){
                string var_name = transmap.get_node(variant_id).name;
                size_t v = stoull(var_name.substr(1));
                var_ids += variant_graph.vcf_records.at(v).id + ' ';
            });

            divergence_file << transmap.get_node(sample_id).name << ',' << transmap.get_node(path_id).name << ',' << avg << ',' << var_ids << '\n';
        }
    }
}


void merge_thread_fn(
        unordered_map<Region,vector<VcfRecord> >& region_records,
        const unordered_map<string,vector<interval_t> >& contig_tandems,
        unordered_map<Region,TransMap>& region_transmaps,
        const unordered_map<string,string>& ref_sequences,
        const vector<Region>& regions,
        const VcfReader& vcf_reader,
        const path& output_dir,
        const HapestryConfig& hapestry_config,
        const OptimizerConfig& opt_config,
        atomic<size_t>& job_index
) {
    // For brevity unpack some of these vars
    auto flank_length = hapestry_config.flank_length;
    auto min_sv_length = hapestry_config.min_sv_length;
    auto graphaligner_timeout = hapestry_config.graphaligner_timeout;

    Timer t;
    Timer t_total;

    // TODO: make this a parameter
    int64_t min_path_coverage = 1;

    unordered_map<string, interval_tree_t<int32_t> > contig_interval_trees;

    for (const auto &[contig, intervals]: contig_tandems) {
        for (const auto &interval: intervals) {
            contig_interval_trees[contig].insert({interval.first, interval.second});
        }
    }

    size_t i = job_index.fetch_add(1);

    path vcf;
    vcf_reader.get_file_path(vcf);
    string vcf_name_prefix = get_vcf_name_prefix(vcf);

    // Thread jobs are regions
    while (i < regions.size()) {
        t.reset();
        t_total.reset();

        const auto &region = regions.at(i);
        auto& transmap = region_transmaps.at(region);

        path subdir = output_dir / region.to_unflanked_string('_', flank_length);

        auto records = region_records.at(region);

        create_directories(subdir);

        path gfa_path = subdir / "graph.gfa";
        path fasta_filename = subdir / "sequences.fasta";

        t.reset();
        VariantGraph variant_graph(ref_sequences, contig_tandems, min_sv_length);

        // Check if the region actually contains any usable variants, and use corresponding build() functions
        if (variant_graph.would_graph_be_nontrivial(records)) {
            variant_graph.build(records, int32_t(flank_length), numeric_limits<int32_t>::max(),
                                region.start + flank_length, region.stop - flank_length, false);
        } else {
            // VariantGraph assumes that the flank length needs to be added to the region
            variant_graph.build(region.name, region.start + flank_length, region.stop - flank_length, flank_length);
        }
        write_time_log(subdir, "variant_graph", t, true);

        cerr << "WRITING GFA to file: " << gfa_path << '\n';
        variant_graph.to_gfa(gfa_path);

        // Write a simple csv for viewing in Bandage
        if (not hapestry_config.skip_nonessential_logs){
            path nodes_csv = subdir / "nodes.csv";

            ofstream nodes_file(nodes_csv);

            nodes_file << "name,is_ref,color\n";
            variant_graph.graph.for_each_handle([&](const handle_t& h){
                nid_t id = variant_graph.graph.get_id(h);
                bool is_ref = variant_graph.is_reference_node(id);

                string color = is_ref ? "#6495ED" : "#BEBEBE";

                nodes_file << id << ',' << is_ref << ',' << color << '\n';
            });

            nodes_file.close();
        }

        path gaf_path = subdir / "alignments.gaf";

        auto name_prefix = get_vcf_name_prefix(vcf);

        path fasta_path = subdir / fasta_filename;

        string command = "GraphAligner"
                         " --seeds-mum-count " "-1"
                         " --seeds-mxm-windowsize " "0"
                         " -b " "10"
                         " --max-cluster-extend " "10"
                         " --multimap-score-fraction  " "1"
                         " -t " "1"
                         " -a " + gaf_path.string() +
                         " -g " + gfa_path.string() +
                         " -f " + fasta_path.string();

        t.reset();

        // Run GraphAligner and check how long it takes, if it times out
        bool success = run_command(command, false, float(graphaligner_timeout));

        write_time_log(subdir, "graphaligner", t, success);

        // Skip remaining steps for this region/tool if alignment failed and get the next job index for the thread
        if (not success) {
            cerr << "WARNING: Command timed out: " << command << '\n';
            write_time_log(subdir, "window_total", t_total, false);

            transmap = {};
            i = job_index.fetch_add(1);
            continue;
        }

        // Iterate the alignments and accumulate their coverages
        unordered_map<string,int64_t> path_coverage;
        get_path_coverages(gaf_path, variant_graph, path_coverage);

        // Add the paths to the TransMap
        for (const auto& [path_name, coverage]: path_coverage) {
            if (coverage >= min_path_coverage){
                transmap.add_path(path_name);
            }
            else{
                cerr << "WARNING: skipping path with low coverage: " << path_name << '\n';
            }
        }

        t.reset();

        // Align all reads to all candidate paths and update transmap
        align_reads_vs_paths(transmap, variant_graph, opt_config.min_read_hap_identity, flank_length, subdir / "pre_optimization_alignments");

        write_time_log(subdir, "align_reads_to_paths", t, true);

        if (not hapestry_config.skip_nonessential_logs){
            // Write the full transmap to CSV (in the form of annotated edges) BEFORE RESCALING WEIGHTS!
            path output_csv = subdir / "reads_to_paths.csv";
            transmap.write_edge_info_to_csv(output_csv, variant_graph, hapestry_config.obscure_sample_names);
        }

        if (opt_config.rescale_weights) {
            t.reset();
            // Rescale the edge weights for each read as a quadratic function of the difference from the best weight
            rescale_weights_as_quadratic_best_diff(transmap, 0, 2000);
            write_time_log(subdir, "rescale_weights", t, true);
        }

        if (not opt_config.skip_solve){
            CPPTRACE_TRY {
                TerminationReason termination_reason;

                if (opt_config.samplewise) {
                    // Split out all the paths per sample so that every sample is an independent connected component
                    unordered_map<string,string> hapmap;
                    transmap.detangle_sample_paths(hapmap);

                    termination_reason = optimize(transmap, opt_config, subdir, true);

                    // Reverse the detangling step
                    transmap.retangle_sample_paths(hapmap);
                }
                else {

                    if (opt_config.use_compression) {
                        termination_reason = optimize_compressed(transmap, opt_config, subdir, true);
                    }
                    else {
                        termination_reason = optimize(transmap, opt_config, subdir, true);
                    }
                }

                // Handle timeout case (do nothing)
                if (termination_reason == TerminationReason::kNoSolutionFound or termination_reason == TerminationReason::kFeasible) {
                    cerr << "WARNING: solver timed out: " << region.to_unflanked_string(':',flank_length) << '\n';
                }
                // Handle failure case (ERROR)
                else if (termination_reason != TerminationReason::kOptimal) {
                    throw runtime_error("ERROR: solver failed with reason " + termination_reason_to_string(termination_reason) + " at: " + region.to_unflanked_string(':',flank_length));
                }
                // Normal operation, solver succeeded
                else {
                    vector<int64_t> unused_paths;

                    transmap.for_each_path([&](const string& path_name, int64_t path_id) {
                        bool used = false;
                        transmap.for_each_sample_of_path(path_name, [&](const string& sample_name, int64_t sample_id) {
                            used = true;
                        });

                        if (not used) {
                            unused_paths.push_back(path_id);
                        }
                    });

                    for (const auto& p: unused_paths) {
                        transmap.remove_node(p);
                    }

                    // Add all the variant nodes to the transmap using a simple name based on the variantgraph ID which
                    // likely does not conflict with existing names
                    for (size_t v=0; v<variant_graph.vcf_records.size(); v++){
                        transmap.add_variant("v" + to_string(v));
                    }

                    // Construct mapping of paths to variants
                    transmap.for_each_path([&](const string& path_name, int64_t path_id){
                        // Convert the path to a vector
                        vector <pair <string,bool> > path;

                        GafAlignment::parse_string_as_path(path_name, path);

                        variant_graph.for_each_vcf_record(path, [&](size_t id, const vector<edge_t>& edges_of_the_record, const VcfRecord& record){
                            string var_name = "v" + to_string(id);
                            transmap.add_edge(var_name, path_name);
                        });
                    });

                    // Write the solution to a VCF
                    path output_path = subdir / "solution.vcf";

                    write_solution_to_vcf(variant_graph, transmap, vcf_reader.sample_ids, output_path, region);

                    output_path = subdir / "solution_haps.vcf";

                    if (hapestry_config.write_hap_vcf) {
                        write_solution_as_hap_vcf(
                            variant_graph,
                            transmap,
                            vcf_reader.sample_ids,
                            ref_sequences,
                            output_path,
                            region,
                            flank_length
                        );
                    }
                }
            } CPPTRACE_CATCH(const std::exception& e) {
                cerr << "ERROR caught at " << region.to_unflanked_string(':',flank_length) << '\n';
                std::cerr<<"Exception: "<<e.what()<<std::endl;
                cpptrace::from_current_exception().print_with_snippets();
            }
        }

        if (HAPESTRY_DEBUG){
            if (transmap.get_is_compressed()) {
                cerr << "WARNING: cannot iterate read->sample mappings for (de)compressed transmap, skipping debug output..." << '\n';
            }
            else {
                write_sample_path_divergence(variant_graph, transmap, subdir);
                align_read_path_edges_of_transmap(transmap, variant_graph, opt_config.min_read_hap_identity, flank_length, subdir / "post_optimization_alignments");
            }
        }

        // Attempt to manually deallocate transmap for this region
        transmap = {};

        write_time_log(subdir, "window_total", t_total, true);

        i = job_index.fetch_add(1);
    }
}


/**
 * @param min_sv_length only variants that affect at least this number of bps are merged; shorter variants are used to
 * build graphs and haplotypes, but they are not merged or printed in output.
 */
void merge_variants(
        const unordered_map <string, interval_tree_t<int32_t> >& contig_interval_trees,
        const unordered_map<string,vector<interval_t> >& contig_tandems,
        unordered_map<Region,TransMap>& region_transmaps,
        const unordered_map<string,string>& ref_sequences,
        const vector<Region>& regions,
        const path& vcf,
        size_t n_threads,
        const HapestryConfig& config,
        const OptimizerConfig& optimizer_config,
        const path& output_dir,
        Timer& t
){
    cerr << t << "Loading VCF records into windows" << '\n';

    unordered_map<Region,vector<VcfRecord> > region_records;
    region_records.reserve(regions.size());

    // Load records for this VCF
    VcfReader vcf_reader(vcf);
    vcf_reader.min_qual = numeric_limits<float>::min();
    vcf_reader.min_sv_length = 1;  // Every record is loaded and used to build the graph, including SNPs/small indels.
    vcf_reader.progress_n_lines = 100'000;
    coord_t record_coord;

    cerr << "Reading VCF... " << '\n';
    vcf_reader.for_record_in_vcf([&](VcfRecord& r){
        // TODO: allow breakends in evaluation
        if (r.sv_type == VcfReader::TYPE_BREAKEND){
            cerr << "WARNING: skipping breakend"  << '\n';
            return;
        }

        r.get_reference_coordinates(false, record_coord);

        auto result = contig_interval_trees.find(r.chrom);

        // First make sure there are actually some windows in this contig (might fail with small satellite contigs)
        if (result == contig_interval_trees.end()){
            return;
        }

        // For each overlapping region, put the VcfRecord in that region
        result->second.overlap_find_all({record_coord.first, record_coord.second}, [&](auto iter){
            coord_t unflanked_window = {iter->low() + config.flank_length, iter->high() - config.flank_length};

            // Skip large events in the population
            // TODO: address these as breakpoints in the VariantGraph and avoid constructing windows as intervals
            // for very large events
            if (record_coord.second - record_coord.first > config.interval_max_length){
                return true;
            }

            // Check if this record exceeds the region
            if (record_coord.first < unflanked_window.first or record_coord.second > unflanked_window.second){
                cerr << "WARNING: skipping record that exceeds the un-flanked window. Record: " << record_coord.first << ',' << record_coord.second << " window: " << unflanked_window.first << ',' << unflanked_window.second << '\n';
                return true;
            }

            Region region(r.chrom, iter->low(), iter->high());
            region_records[region].emplace_back(r);
            return true;
        });
    });

    // Before moving on, make sure every region has at least an empty vector
    for (const auto& r: regions){
        if (region_records.find(r) == region_records.end()){
            region_records[r] = {};
        }
    }

    cerr << t << "Peak memory usage: " << get_peak_memory_usage() << '\n';
    cerr << t << "Launching threads for graph alignment and optimization: " << vcf << '\n';

    // Convert VCFs to graphs and run graph aligner and build/solve the model for each region
    {
        // Thread-related variables
        atomic<size_t> job_index = 0;
        vector<thread> threads;

        threads.reserve(n_threads);

        // Launch threads
        for (size_t n=0; n<n_threads; n++) {
            try {
                cerr << "launching: " << n << '\n';
                threads.emplace_back(merge_thread_fn,
                     std::ref(region_records),
                     std::cref(contig_tandems),
                     std::ref(region_transmaps),
                     std::cref(ref_sequences),
                     std::cref(regions),
                     std::cref(vcf_reader),
                     std::cref(output_dir),
                     std::cref(config),
                     std::cref(optimizer_config),
                     std::ref(job_index)
                );
            } catch (const exception &e) {
                throw e;
            }
        }

        // Wait for threads to finish
        for (auto &n: threads) {
            n.join();
        }
    }
}


void concatenate_output_vcfs(
        const path& out_vcf,
        const path& vcf,
        const path& output_dir,
        vector<Region>& regions,
        int32_t flank_length,
        const string& fail_log_suffix = "failed",
        const string& vcf_prefix = "solution"
        ) {

    ofstream out_file(out_vcf);

    if (not (out_file.is_open() and out_file.good())){
        throw runtime_error("ERROR: could not write file: " + out_vcf.string());
    }

    ifstream input_vcf(vcf);

    if (not (input_vcf.is_open() and input_vcf.good())){
        throw runtime_error("ERROR: could not write file: " + vcf.string());
    }

    // Just copy over the header lines
    string line;
    while (getline(input_vcf, line)){
        if (line.starts_with("##")){
            out_file << line << '\n';
        }
        else{
            // The next non-## line MUST be the sample line
            if (line.starts_with("#")){
                // Write header lines for the INFO tags that we added here
                out_file << VcfReader::INFO_REDUNDANT_HEADER << '\n';

                // Copy the input sample line to the output VCF
                out_file << line << '\n';
            }
            else {
                throw runtime_error("ERROR: cannot find sample line in input VCF: " + vcf.string());
            }

            // Stop iteration
            break;
        }
    }

    input_vcf.close();

    path fail_regions_bed_path = output_dir / ("windows_" + fail_log_suffix + ".bed");
    ofstream fail_regions_file(fail_regions_bed_path);

    path fail_regions_flanked_bed_path = output_dir / ("windows_" + fail_log_suffix + "_flanked.bed");
    ofstream fail_regions_flanked_file(fail_regions_flanked_bed_path);

    if (not (fail_regions_file.is_open() and fail_regions_file.good())) {
        throw runtime_error("ERROR: could not write file: " + fail_regions_bed_path.string());
    }

    if (not (fail_regions_flanked_file.is_open() and fail_regions_flanked_file.good())) {
        throw runtime_error("ERROR: could not write file: " + fail_regions_flanked_bed_path.string());
    }

    // Copy over the mutable parts of the header and then the main contents of the filtered VCF
    // Will be empty if no regions were processed
    for (size_t i=0; i<regions.size(); i++){
        const auto& region = regions[i];
        path sub_vcf = output_dir / region.to_unflanked_string('_', flank_length) / (vcf_prefix + ".vcf");

        // Skip if the file does not exist
        if (not exists(sub_vcf)){
            // Log that this region did not contain any solution
            fail_regions_file << region.to_unflanked_bed(flank_length) << '\n';
            fail_regions_flanked_file << region.to_bed() << '\n';
            continue;
        }

        ifstream file(sub_vcf);
        while (getline(file, line)){
            if (not line.starts_with('#')){
                out_file << line << '\n';
            }
        }
    }
}


void hapestry(
        path vcf,
        path output_dir,
        path windows_bed,                // Override the interval graph if this is provided
        path tandem_bed,
        path bam_csv,
        path ref_fasta,
        int32_t n_threads,
        const HapestryConfig& hapestry_config,
        const OptimizerConfig& optimizer_config
){
    Timer t;

    if (std::filesystem::exists(output_dir)){
        throw runtime_error("ERROR: output dir exists already: " + output_dir.string());
    }
    else{
        std::filesystem::create_directories(output_dir);
    }

    output_dir = std::filesystem::weakly_canonical(output_dir);
    tandem_bed = std::filesystem::weakly_canonical(tandem_bed);
    bam_csv = std::filesystem::weakly_canonical(bam_csv);
    ref_fasta = std::filesystem::weakly_canonical(ref_fasta);
    vcf = std::filesystem::weakly_canonical(vcf);

    cerr << t << "Initializing" << '\n';

    vector <pair <string,path> > bam_paths;
    unordered_map<string,string> ref_sequences;
    vector<Region> regions;

    cerr << t << "Peak memory usage: " << get_peak_memory_usage() << '\n';
    cerr << t << "Loading reference sequences" << '\n';

    // Load all chromosome sequences (in case of BND)
    for_sequence_in_fasta_file(ref_fasta, [&](const Sequence& s){
        ref_sequences[s.name] = s.sequence;
    });

    cerr << t << "Peak memory usage: " << get_peak_memory_usage() << '\n';
    cerr << "Reading tandem BED" << '\n';

    unordered_map<string,vector<interval_t> > contig_tandems;
    interval_t interval;
    for_region_in_bed_file(tandem_bed, [&](const Region& r){
        interval.first = r.start;
        interval.second = r.stop;
        contig_tandems[r.name].emplace_back(interval);
    });

    if (windows_bed.empty()){
        cerr << t << "Constructing windows from VCFs and tandem BED" << '\n';
        path bed_log_path = output_dir / "windows_omitted.bed";
        construct_windows_from_vcf_and_bed(ref_sequences, contig_tandems, {vcf}, hapestry_config.flank_length, hapestry_config.interval_max_length, hapestry_config.min_sv_length, regions, bed_log_path, false);
    }
    else {
        cerr << t << "Reading BED file" << '\n';
        load_windows_from_bed(windows_bed, regions);
    }

    if (not regions.empty()){
        // This is only used while loading VCFs to find where each record belongs
        unordered_map <string, interval_tree_t<int32_t> > contig_interval_trees;

        // Log which windows were used
        path bed_output_path = output_dir / "windows.bed";
        path bed_flanked_output_path = output_dir / "windows_flanked.bed";
        ofstream output_bed_file(bed_output_path);
        ofstream output_bed_flanked_file(bed_flanked_output_path);

        cerr << t << "Peak memory usage: " << get_peak_memory_usage() << '\n';
        cerr << t << "Flanking windows and writing BED" << '\n';

        // Add flanks, place the regions in the interval tree, and log the windows
        for (auto& r: regions) {
            output_bed_file << r.to_bed() << '\n';

            r.start -= hapestry_config.flank_length;
            r.stop += hapestry_config.flank_length;

            contig_interval_trees[r.name].insert({r.start, r.stop});
            output_bed_flanked_file << r.to_bed() << '\n';
        }
        output_bed_file.close();
        output_bed_flanked_file.close();

        cerr << t << "Peak memory usage: " << get_peak_memory_usage() << '\n';
        cerr << t << "Forgetting unused reference contigs" << '\n';
        vector<string> unused_contigs;

        for (const auto& [name, seq]: ref_sequences) {
            if (not contig_interval_trees.contains(name)) {
                unused_contigs.emplace_back(name);
            }
        }

        for (const auto& name: unused_contigs) {
            ref_sequences.erase(name);
        }

        cerr << t << "Peak memory usage: " << get_peak_memory_usage() << '\n';
        cerr << t << "Fetching reads for all windows" << '\n';

        // The container to store all fetched reads and their relationships to samples/paths
        unordered_map<Region,TransMap> region_transmaps;

        auto max_length = size_t(float(hapestry_config.interval_max_length) * 2.5);

        if (hapestry_config.bam_not_hardclipped){
            cerr << "Fetching from NON-hardclipped BAMs" << '\n';
            fetch_reads(
                    t,
                    regions,
                    bam_csv,
                    n_threads,
                    region_transmaps,
                    true,
                    hapestry_config.force_unique_reads,
                    true,
                    false,
                    hapestry_config.flank_length
            );
        }
        else{
            cerr << "Fetching from HARDCLIPPED BAMs" << '\n';
            fetch_reads_from_clipped_bam(
                    t,
                    regions,
                    bam_csv,
                    n_threads,
                    max_length,
                    hapestry_config.flank_length,
                    region_transmaps,
                    true,
                    false,
                    false,
                    hapestry_config.force_unique_reads,
                    true,
                    hapestry_config.flank_length
            );
        }

        cerr << t << "Peak memory usage: " << get_peak_memory_usage() << '\n';
        cerr << t << "Writing sequences to disk" << '\n';

        path fasta_filename = "sequences.fasta";
        path staging_dir = output_dir;

        // Dump sequences into each region directory
        {
            // Thread-related variables
            atomic<size_t> job_index = 0;
            vector<thread> threads;

            threads.reserve(n_threads);

            // Launch threads
            for (size_t n=0; n<n_threads; n++) {
                try {
                    cerr << "launching: " << n << '\n';
                    threads.emplace_back(write_region_subsequences_to_file_thread_fn,
                                         std::cref(region_transmaps),
                                         std::cref(regions),
                                         std::cref(staging_dir),
                                         std::cref(fasta_filename),
                                         hapestry_config.flank_length,
                                         true,
                                         std::ref(job_index)
                    );
                } catch (const exception &e) {
                    throw e;
                }
            }

            // Wait for threads to finish
            for (auto &n: threads) {
                n.join();
            }
        }


        // Generate GFAs/GAFs/CSVs and folder structure for every VCF * every region
        // TODO: create option to use /dev/shm/ as staging dir. Absolutely must delete the /dev/shm/ copy or warn the user at termination

        merge_variants(
                contig_interval_trees,
                contig_tandems,
                region_transmaps,
                ref_sequences,
                regions,
                vcf,
                n_threads,
                hapestry_config,
                optimizer_config,
                staging_dir,
                t
        );

        auto vcf_prefix = get_vcf_name_prefix(vcf);
        path out_vcf = output_dir / ("merged.vcf");

        cerr << t << "Writing output VCF: " << out_vcf << '\n';
        concatenate_output_vcfs(out_vcf, vcf, output_dir, regions, hapestry_config.flank_length);

        // Additionally write a window-length haplotype substitution VCF if the user requested one
        if (hapestry_config.write_hap_vcf) {
            out_vcf = output_dir / ("merged_hap.vcf");
            cerr << t << "Writing output VCF: " << out_vcf << '\n';
            concatenate_output_vcfs(out_vcf, vcf, output_dir, regions, hapestry_config.flank_length, "failed_haps", "solution_haps");
        }

        if (hapestry_config.skip_nonessential_logs) {
            // Delete all the intermediate VCFs (this should really be done with a proper DB but for now we just remove them after concatenating them)
            for (size_t i = 0; i < regions.size(); i++) {
                const auto &region = regions[i];
                path sub_vcf = output_dir / region.to_unflanked_string('_', hapestry_config.flank_length) / (vcf_prefix + ".vcf");

                // Remove if the VCF was generated
                if (exists(sub_vcf)) {
                    std::filesystem::remove(sub_vcf);
                    continue;
                }
            }
        }

    }
    else{
        cerr << "WARNING: No regions to process" << '\n';
    }

    cerr << t << "Peak memory usage: " << get_peak_memory_usage() << '\n';
    cerr << t << "Done" << '\n';
}


int main (int argc, char* argv[]){
    path output_dir;
    path windows_bed;
    path tandem_bed;
    string bam_csv;
    path ref_fasta;
    path vcf;
    int32_t n_threads = 1;
    OptimizerConfig optimizer_config;
    HapestryConfig hapestry_config;
    bool use_gurobi = false;

    CLI::App app{"App description"};

    app.add_option(
            "--n_threads",
            n_threads,
            "Maximum number of threads to use");

    app.add_option(
            "--output_dir",
            output_dir,
            "Path to output directory which must not exist yet")
            ->required();

    app.add_option(
            "--windows",
            windows_bed,
            "Path to BED file containing windows to merge (inferred automatically if not provided)");

    app.add_option(
            "--tandems",
            tandem_bed,
            "Path to BED file containing tandem track which will inform how to aggregate variants in windows");

    app.add_option(
            "--vcf",
            vcf,
            "Path to VCF file containing variants to annotate (must be biallelic/normed)");

    app.add_option(
            "--bam_csv",
            bam_csv,
            "Simple headerless CSV file with the format [sample_name],[bam_path]")
            ->required();

    app.add_option(
            "--ref",
            ref_fasta,
            "Path to reference sequence FASTA file that corresponds to VCF")
            ->required();

    app.add_option(
            "--flank_length",
            hapestry_config.flank_length,
            "How much flanking sequence to use when fetching and aligning reads")
            ->required();

    app.add_option(
            "--interval_max_length",
            hapestry_config.interval_max_length,
            "How long a window can be in bp before it is skipped")
            ->required();

    app.add_option(
            "--min_sv_length",
            hapestry_config.min_sv_length,
            "Only variants that affect at least this number of bps are merged. Shorter variants are used to build graphs and haplotypes, but they are not merged or printed in output.")
            ->required();

    app.add_option(
            "--graphaligner_timeout",
            hapestry_config.graphaligner_timeout,
            "Abort graphaligner after this many seconds, and do not compute the remaining steps for that window");

    app.add_option(
            "--solver_timeout",
            optimizer_config.timeout_sec,
            "Abort the optimizer after this many seconds, use 0 for no limit");

    app.add_option(
            "--min_read_hap_identity",
            optimizer_config.min_read_hap_identity,
            "Minimum alignment identity to consider a read-hap assignment in the optimization step")
            ->required();

    app.add_option(
            "--d_weight",
            optimizer_config.d_weight,
            "Scalar coefficient to apply to d_norm^2 in the optimization step, used to priotize or deprioritize the distance term");

    app.add_flag("--skip_solve", optimizer_config.skip_solve, "Invoke this to skip the optimization step. CSVs for each optimization input will still be written.");

    app.add_flag("--rescale_weights", optimizer_config.rescale_weights, "Invoke this to use quadratic difference-from-best match rescaling for read-to-path edges.");

    app.add_flag("--debug", HAPESTRY_DEBUG, "Invoke this to add more logging and output");

    app.add_flag("--use_gurobi", use_gurobi, "Invoke this to use Gurobi instead of SCIP. License must be in conventional location.");

    app.add_flag("!--no_sum_constraints", optimizer_config.use_sum_constraints, "Use individual constraints instead of a sum for 'any of' implications (strongly recommended NOT to use sums for SCIP, RAM+CPU explodes with sums)");

    app.add_flag("--force_unique_reads", hapestry_config.force_unique_reads, "Invoke this to add append each read name with the sample name so that inter-sample read collisions cannot occur");

    app.add_flag("--skip_nonessential_logs", hapestry_config.skip_nonessential_logs, "Invoke this to skip logs: reads_to_paths.csv, solution.csv, nodes.csv");

    app.add_flag("--bam_not_hardclipped", hapestry_config.bam_not_hardclipped, "Invoke this if you expect your BAMs NOT to contain ANY hardclips. Saves time on iterating.");

    app.add_flag("--write_hap_vcf", hapestry_config.write_hap_vcf, "Invoke this to write an additional VCF which only writes solutions in the form of full window length haplotypes (replacement/substitution operations).");

    app.add_flag("--samplewise", optimizer_config.samplewise, "Use samplewise solver instead of global solver");

    app.add_flag("--compress", optimizer_config.use_compression, "Use reversible compression methods to simplify the TransMap (problem input) before solving with LP solver");

    app.add_option("--compress_quantum", optimizer_config.compress_quantum, "Constant by which to quantize the weight comparisons. By default (quantum=0) edges may have a maximum of up to 1000 different values. E.g. quantum=2, possible values = 500");

    app.add_flag("--prune_with_d_min", optimizer_config.prune_with_d_min, "Use d_min solution to remove all edges not used");

    app.add_flag("--obscure_sample_names_from_csv", hapestry_config.obscure_sample_names, "Don't write sample names to reads_to_paths.csv. Instead, write an arbitrarily determined integer ID.");

    try{
        app.parse(argc, argv);
    }
    catch (const CLI::ParseError &e) {
        return app.exit(e);
    }

    if (HAPESTRY_DEBUG){
        cerr << DEBUG_BANNER;
    }

    // We don't want to make the user try to guess the names of these variables, so this is a simple boolean flag
    if (use_gurobi) {
        optimizer_config.solver_type = SolverType::kGurobi;
    }

    hapestry(
            vcf,
            output_dir,
            windows_bed,
            tandem_bed,
            bam_csv,
            ref_fasta,
            n_threads,
            hapestry_config,
            optimizer_config
    );

    return 0;
}
