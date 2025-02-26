#include "Region.hpp"
#include "VcfReader.hpp"
#include "bed.hpp"
#include "CLI11.hpp"
#include "gaf.hpp"

using sv_merge::VcfRecord;
using sv_merge::VcfReader;
using sv_merge::Region;
using sv_merge::GafAlignment;

#include <stdexcept>
#include <filesystem>

using std::filesystem::path;
using std::filesystem::directory_iterator;
using std::filesystem::exists;
using std::filesystem::create_directory;
using std::numeric_limits;
using std::streamsize;
using std::string;
using std::vector;
using std::unordered_map;
using std::set;
using std::ifstream;
using std::runtime_error;
using std::sort;
using std::reverse;
using std::to_string;
using std::max;
using std::min;
using std::ofstream;
using std::pair;


using namespace sv_merge;


/**
 * Reads window directories and prints the values of all windows
 */
class Counts {
public:
    static const char SUBDIR_SEPARATOR_1 = '_';
    static const char SUBDIR_SEPARATOR_2 = '-';

    /**
     * Allocates output arrays
     *
     * @param coverage_threshold fraction above which a node or haplotype is considered fully covered;
     * @param n_windows an estimate on the number of windows that will be processed; used just to allocate space;
     * @param log_* keeps track of every window that satisfies at least one of these low-support thresholds.
     */
    explicit Counts(const vector<string>& tools, const string& truth_tool, double coverage_threshold, size_t max_hap_length, double log_nodes_fully_covered_leq, double log_edges_covered_leq, double log_vcf_records_supported_leq, double log_alignment_identity_leq, double log_hap_coverage_leq, size_t min_sv_length, size_t n_haplotypes_expected, size_t n_windows = 1e6):
        TOOLS(tools),
        TRUTH_TOOL(truth_tool),
        N_TOOLS(tools.size()),
        coverage_threshold(coverage_threshold),
        max_hap_length(max_hap_length),
        log_nodes_fully_covered_leq(log_nodes_fully_covered_leq),
        log_edges_covered_leq(log_edges_covered_leq),
        log_vcf_records_supported_leq(log_vcf_records_supported_leq),
        log_alignment_identity_leq(log_alignment_identity_leq),
        log_hap_coverage_leq(log_hap_coverage_leq),
        min_sv_length(min_sv_length),
        N_HAPLOTYPES_EXPECTED(n_haplotypes_expected)
    {
        size_t i;

        n_alignments.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { n_alignments.emplace_back(); n_alignments.at(i).reserve(n_windows); }

        n_haplotypes.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { n_haplotypes.emplace_back(); n_haplotypes.at(i).reserve(n_windows); }
        n_nonref_haplotypes.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { n_nonref_haplotypes.emplace_back(); n_nonref_haplotypes.at(i).reserve(n_windows); }
        n_haplotype_clusters.reserve(n_windows);

        haplotype_coverage_avg.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { haplotype_coverage_avg.emplace_back(); haplotype_coverage_avg.at(i).reserve(n_windows); }
        alignment_identity_avg.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { alignment_identity_avg.emplace_back(); alignment_identity_avg.at(i).reserve(n_windows); }

        nonref_haplotype_coverage_avg.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { nonref_haplotype_coverage_avg.emplace_back(); nonref_haplotype_coverage_avg.at(i).reserve(n_windows); }
        nonref_alignment_identity_avg.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { nonref_alignment_identity_avg.emplace_back(); nonref_alignment_identity_avg.at(i).reserve(n_windows); }

        cluster_coverage_avg.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { cluster_coverage_avg.emplace_back(); cluster_coverage_avg.at(i).reserve(n_windows); }
        cluster_alignment_identity_avg.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { cluster_alignment_identity_avg.emplace_back(); cluster_alignment_identity_avg.at(i).reserve(n_windows); }

        nonref_nodes.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { nonref_nodes.emplace_back(); nonref_nodes.at(i).reserve(n_windows); }
        nonref_nodes_small.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { nonref_nodes_small.emplace_back(); nonref_nodes_small.at(i).reserve(n_windows); }
        nonref_nodes_large.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { nonref_nodes_large.emplace_back(); nonref_nodes_large.at(i).reserve(n_windows); }
        nonref_nodes_fully_covered.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { nonref_nodes_fully_covered.emplace_back(); nonref_nodes_fully_covered.at(i).reserve(n_windows); }
        nonref_nodes_fully_covered_small.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { nonref_nodes_fully_covered_small.emplace_back(); nonref_nodes_fully_covered_small.at(i).reserve(n_windows); }
        nonref_nodes_fully_covered_large.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { nonref_nodes_fully_covered_large.emplace_back(); nonref_nodes_fully_covered_large.at(i).reserve(n_windows); }
        nonref_nodes_partially_covered.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { nonref_nodes_partially_covered.emplace_back(); nonref_nodes_partially_covered.at(i).reserve(n_windows); }
        nonref_nodes_bps.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { nonref_nodes_bps.emplace_back(); nonref_nodes_bps.at(i).reserve(n_windows); }
        nonref_nodes_bps_covered.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { nonref_nodes_bps_covered.emplace_back(); nonref_nodes_bps_covered.at(i).reserve(n_windows); }
        nonref_edges.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { nonref_edges.emplace_back(); nonref_edges.at(i).reserve(n_windows); }
        nonref_edges_covered.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { nonref_edges_covered.emplace_back(); nonref_edges_covered.at(i).reserve(n_windows); }
        nonref_node_length_vs_coverage.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { nonref_node_length_vs_coverage.emplace_back(); nonref_node_length_vs_coverage.at(i).reserve(n_windows); }
        nonref_node_length_vs_identity.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { nonref_node_length_vs_identity.emplace_back(); nonref_node_length_vs_identity.at(i).reserve(n_windows); }
        supported_vcf_records.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { supported_vcf_records.emplace_back(); supported_vcf_records.at(i).reserve(n_windows); }
        unsupported_vcf_records.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { unsupported_vcf_records.emplace_back(); unsupported_vcf_records.at(i).reserve(n_windows); }

        supported_del.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { supported_del.emplace_back(); supported_del.at(i).reserve(n_windows); }
        unsupported_del.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { unsupported_del.emplace_back(); unsupported_del.at(i).reserve(n_windows); }
        supported_ins.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { supported_ins.emplace_back(); supported_ins.at(i).reserve(n_windows); }
        unsupported_ins.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { unsupported_ins.emplace_back(); unsupported_ins.at(i).reserve(n_windows); }
        supported_inv.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { supported_inv.emplace_back(); supported_inv.at(i).reserve(n_windows); }
        unsupported_inv.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { unsupported_inv.emplace_back(); unsupported_inv.at(i).reserve(n_windows); }
        supported_dup.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { supported_dup.emplace_back(); supported_dup.at(i).reserve(n_windows); }
        unsupported_dup.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { unsupported_dup.emplace_back(); unsupported_dup.at(i).reserve(n_windows); }
        sv_length_supported.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { sv_length_supported.emplace_back(); sv_length_supported.at(i).reserve(n_windows); }
        sv_length_unsupported.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { sv_length_unsupported.emplace_back(); sv_length_unsupported.at(i).reserve(n_windows); }

        recall_numerator.reserve(N_TOOLS); recall_denominator.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { recall_numerator.emplace_back(); recall_denominator.emplace_back(); }
    };


    /**
     * Adds to the lists of measures all the values in the directory of a window.
     *
     * @param directory assumed structure (every one of these files must be present, unless otherwise noted):
     *
     * DIRECTORY
     * ├── `CLUSTERS_FILE`: for every cluster (row): `cluster_id,haplotype_ìds`, where `haplotype_ìds` is a space-
     * │   separated list; this file might be absent, and some haplotypes might be missing; `cluster_id` might be a GAF
     * │   path in some graph (not necessarily in the local graph), but the procedure does not try to parse it.
     * ├── TOOL_ID_1
     * │   ├── `NODES_FILE`: for every node (rows): `id,length,is_reference (0/1),coverage,identity`.
     * │   ├── `HAPLOTYPES_FILE`: for every haplotype (row): `id,length,is_reference (0/1),coverage,identity`.
     * │   ├── `OTHER_COUNTS_FILE`: `total_n_alignments,n_nonreference_edges,n_nonreference_edges_covered`.
     * │   ├── `SUPPORTED_VCF_FILE`: contains only calls that are fully supported by some alignment; can be empty;
     * │   └── `UNSUPPORTED_VCF_FILE`: contains only calls that are not fully supported by any alignment; can be empty.
     * ├── TOOL_ID_2
     * │   ├── ...
     * ... ...
     *
     * @return FALSE if the optional clusters file is missing from `directory`.
     */
    bool load_window(const path& directory) {
        bool out;
        char c;
        size_t i;
        size_t field;
        path input_file;
        ifstream file;

        n_windows++;

        // Haplotype clusters (optional)
        input_file=directory/CLUSTERS_FILE;
        file.clear(); file.open(input_file);
        on_init_clusters();
        if (!file.good() || !file.is_open()) out=false;
        else {
            out=true;
            file.ignore(STREAMSIZE_MAX,LINE_DELIMITER);  // Skipping CSV header
            field=0;
            while (file.get(c)) {
                if (c==LINE_DELIMITER) { on_field_end_clusters(++field,tmp_buffer_1); on_line_end_clusters(tmp_buffer_2); tmp_buffer_1.clear(); field=0; }
                else if (c==CSV_DELIMITER) { on_field_end_clusters(++field,tmp_buffer_1); tmp_buffer_1.clear(); }
                else tmp_buffer_1.push_back(c);
            }
            if (!tmp_buffer_1.empty()) { on_field_end_clusters(++field,tmp_buffer_1); on_line_end_clusters(tmp_buffer_2); }
            file.close();
        }
        on_window_end_clusters();

        for (i=0; i<N_TOOLS; i++) {
            // Node counts
            input_file=directory/TOOLS.at(i)/NODES_FILE;
            file.clear(); file.open(input_file);
            if (!file.good() || !file.is_open()) throw runtime_error("ERROR: could not read file: "+input_file.string());
            on_init_nodes();
            file.ignore(STREAMSIZE_MAX,LINE_DELIMITER);  // Skipping CSV header
            field=0;
            while (file.get(c)) {
                if (c==LINE_DELIMITER) { on_field_end_nodes(++field,tmp_buffer_1); on_line_end_nodes(i); tmp_buffer_1.clear(); field=0; }
                else if (c==CSV_DELIMITER) { on_field_end_nodes(++field,tmp_buffer_1); tmp_buffer_1.clear(); }
                else tmp_buffer_1.push_back(c);
            }
            if (!tmp_buffer_1.empty()) { on_field_end_nodes(++field,tmp_buffer_1); on_line_end_nodes(i); }
            file.close();
            on_window_end_nodes(i);

            // Haplotype counts
            input_file=directory/TOOLS.at(i)/HAPLOTYPES_FILE;
            file.clear(); file.open(input_file);
            if (!file.good() || !file.is_open()) throw runtime_error("ERROR: could not read file: "+input_file.string());
            on_init_haplotypes();
            file.ignore(STREAMSIZE_MAX,LINE_DELIMITER);  // Skipping CSV header
            field=0;
            while (file.get(c)) {
                if (c==LINE_DELIMITER) { on_field_end_haplotypes(++field,tmp_buffer_1); on_line_end_haplotypes(input_file); tmp_buffer_1.clear(); field=0; }
                else if (c==CSV_DELIMITER) { on_field_end_haplotypes(++field,tmp_buffer_1); tmp_buffer_1.clear(); }
                else tmp_buffer_1.push_back(c);
            }
            if (!tmp_buffer_1.empty()) { on_field_end_haplotypes(++field,tmp_buffer_1); on_line_end_haplotypes(input_file); }
            file.close();
            on_window_end_haplotypes(i,input_file);

            // Edge counts
            input_file=directory/TOOLS.at(i)/EDGE_COUNTS_FILE;
            file.clear(); file.open(input_file);
            if (!file.good() || !file.is_open()) throw runtime_error("ERROR: could not read file: "+input_file.string());
            file.ignore(STREAMSIZE_MAX,LINE_DELIMITER);  // Skipping CSV header
            field=0;
            while (file.get(c)) {
                if (c==LINE_DELIMITER) { on_field_end_edge_counts(++field,tmp_buffer_1); on_window_end_edge_counts(i); tmp_buffer_1.clear(); field=0; }
                else if (c==CSV_DELIMITER) { on_field_end_edge_counts(++field,tmp_buffer_1); tmp_buffer_1.clear(); }
                else tmp_buffer_1.push_back(c);
            }
            if (!tmp_buffer_1.empty()) { on_field_end_edge_counts(++field,tmp_buffer_1); on_window_end_edge_counts(i); }
            file.close();

            // VCF counts
            input_file=directory/TOOLS.at(i)/SUPPORTED_VCF_FILE;
            on_init_vcf_counts();
            VcfReader reader1(input_file);
            reader1.progress_n_lines=0;
            reader1.for_record_in_vcf([&](VcfRecord& record) {
                vcf_counts_supported.at(0)++;
                if (record.sv_length>0 && record.sv_length<INT32_MAX) sv_length_supported.at(i).emplace_back(record.sv_length);
                if (record.sv_type==VcfReader::TYPE_DELETION) vcf_counts_supported.at(1)++;
                else if (record.sv_type==VcfReader::TYPE_INSERTION) vcf_counts_supported.at(2)++;
                else if (record.sv_type==VcfReader::TYPE_INVERSION) vcf_counts_supported.at(3)++;
                else if (record.sv_type==VcfReader::TYPE_DUPLICATION) vcf_counts_supported.at(4)++;
            });
            input_file=directory/TOOLS.at(i)/UNSUPPORTED_VCF_FILE;
            VcfReader reader2(input_file);
            reader2.progress_n_lines=0;
            reader2.for_record_in_vcf([&](VcfRecord& record) {
                vcf_counts_unsupported.at(0)++;
                if (record.sv_length>0 && record.sv_length<INT32_MAX) sv_length_unsupported.at(i).emplace_back(record.sv_length);
                if (record.sv_type==VcfReader::TYPE_DELETION) vcf_counts_unsupported.at(1)++;
                else if (record.sv_type==VcfReader::TYPE_INSERTION) vcf_counts_unsupported.at(2)++;
                else if (record.sv_type==VcfReader::TYPE_INVERSION) vcf_counts_unsupported.at(3)++;
                else if (record.sv_type==VcfReader::TYPE_DUPLICATION) vcf_counts_unsupported.at(4)++;
            });
            on_window_end_vcf_counts(i);
        }

        return out;
    }


    /**
     * Computes a version of recall based on non-ref haplotypes rather than on VCF records. Assume that one of the VCFs
     * represents the true SVs in a window, and call "true graph" the corresponding graph T. A haplotype is considered
     * non-reference iff it does not map to the ref walk in T. Every distinct walk in T taken by some non-ref haplotype
     * contributes one to the denominator of recall (non-ref haplotypes with distinct IDs might map to the same walk).
     * Given another VCF graph G, if a non-ref haplotype maps well to a walk in G, its walk in T is considered to be
     * found in G. The numerator of recall for that VCF is incremented by the number of distinct walks in T that are
     * found in G.
     *
     * @param directory assumed to have the same structure as in `load_window()`.
     */
    void haplotype_recall(const path& directory, const string& truth_tool) {
        const float MIN_COVERAGE = 0.9;  // Arbitrary
        const float MIN_IDENTITY = 0.9;  // Arbitrary
        bool is_ref;
        char c;
        size_t i;
        size_t denominator, field, length, reference_path_length;
        int node_id, previous_id, direction;
        string hap_id, reference_path, canonized_path;
        path input_file;
        ifstream file;

        // Finding ref nodes in the true graph
        ref_node_ids.clear();
        input_file=directory/truth_tool/NODES_FILE;
        file.clear(); file.open(input_file);
        if (!file.good() || !file.is_open()) throw runtime_error("ERROR: could not read file: "+input_file.string());
        file.ignore(STREAMSIZE_MAX,LINE_DELIMITER);  // Skipping CSV header
        field=0; tmp_buffer_1.clear();
        while (file.get(c)) {
            if (c==LINE_DELIMITER) { tmp_buffer_1.clear(); field=0; }
            else if (c==CSV_DELIMITER) {
                if (field==0) node_id=stoi(tmp_buffer_1);
                else if (field==2 && stoi(tmp_buffer_1)==1) ref_node_ids.emplace(node_id);
                field++; tmp_buffer_1.clear();
            }
            else tmp_buffer_1.push_back(c);
        }
        file.close();

        // Finding the reference path, which is assumed to be the longest walk that uses only reference nodes with
        // increasing IDs (this might not exist, since no haplotype might be reference in the window).
        reference_path=""; reference_path_length=0;
        input_file=directory/truth_tool/ALIGNMENTS_FILE;
        file.clear(); file.open(input_file);
        if (!file.good() || !file.is_open()) throw runtime_error("ERROR: could not read file: "+input_file.string());
        field=0; tmp_buffer_1.clear(); tmp_buffer_2.clear();
        while (file.get(c)) {
            if (c==LINE_DELIMITER) { tmp_buffer_1.clear(); field=0; }
            else if (c==TSV_DELIMITER) {
                if (field==0) hap_id=tmp_buffer_1;
                else if (field==5) {
                    GafAlignment::parse_string_as_path(tmp_buffer_1,tmp_path);
                    is_ref=true; previous_id=-1; direction=-1;
                    for (auto& [str, is_reverse]: tmp_path) {
                        node_id=stoi(str);
                        if (!ref_node_ids.contains(node_id)) { is_ref=false; break; }
                        if (previous_id!=-1 && direction==-1) {
                            if (node_id==previous_id+1) direction=1;
                            else if (node_id==previous_id-1) direction=0;
                            else { is_ref=false; break; }
                        }
                        else if ((direction==1 && node_id!=previous_id+1) || (direction==0 && node_id!=previous_id-1)) { is_ref=false; break; }
                        previous_id=node_id;
                    }
                    length=tmp_path.size();
                    if (is_ref && length>reference_path_length) {
                        reference_path_length=length;
                        VariantGraph::canonize_gaf_path(tmp_buffer_1,reference_path,tmp_buffer_2);
                    }
                    break;
                }
                field++; tmp_buffer_1.clear();
            }
            else tmp_buffer_1.push_back(c);
        }
        file.close();

        // Finding non-ref haplotypes using the true graph
        nonref_haps.clear(); canonized_paths.clear(); hap2path.clear();
        input_file=directory/truth_tool/ALIGNMENTS_FILE;
        file.clear(); file.open(input_file);
        if (!file.good() || !file.is_open()) throw runtime_error("ERROR: could not read file: "+input_file.string());
        field=0; tmp_buffer_1.clear(); tmp_buffer_2.clear();
        while (file.get(c)) {
            if (c==LINE_DELIMITER) { tmp_buffer_1.clear(); field=0; }
            else if (c==TSV_DELIMITER) {
                if (field==0) hap_id=tmp_buffer_1;
                else if (field==5) {
                    VariantGraph::canonize_gaf_path(tmp_buffer_1,canonized_path,tmp_buffer_2);
                    if (canonized_path!=reference_path) {
                        canonized_paths.emplace(canonized_path);
                        nonref_haps.emplace(hap_id);
                        hap2path.emplace(hap_id,canonized_path);
                    }
                }
                field++; tmp_buffer_1.clear();
            }
            else tmp_buffer_1.push_back(c);
        }
        file.close();
        denominator=canonized_paths.size();
        if (denominator==0) {
            for (i=0; i<N_TOOLS; i++) {
                recall_numerator.at(i).emplace_back(0);
                recall_denominator.at(i).emplace_back(0);
            }
            return;
        }

        // Incrementing recall numerator and denominator
        tmp_buffer_1.clear();
        for (i=0; i<N_TOOLS; i++) {
            if (TOOLS.at(i)==truth_tool) continue;
            covered_paths.clear();
            input_file=directory/TOOLS.at(i)/HAPLOTYPES_FILE;
            file.clear(); file.open(input_file);
            if (!file.good() || !file.is_open()) throw runtime_error("ERROR: could not read file: "+input_file.string());
            file.ignore(STREAMSIZE_MAX,LINE_DELIMITER);  // Skipping CSV header
            field=0;
            while (file.get(c)) {
                if (c==LINE_DELIMITER) {
                    identity=stof(tmp_buffer_1);
                    if (nonref_haps.contains(hap_id) && coverage>=MIN_COVERAGE && identity>=MIN_IDENTITY) covered_paths.emplace(hap2path.at(hap_id));
                    tmp_buffer_1.clear(); field=0;
                }
                else if (c==CSV_DELIMITER) {
                    if (field==0) hap_id=tmp_buffer_1;
                    else if (field==4) coverage=stof(tmp_buffer_1);
                    tmp_buffer_1.clear(); field++;
                }
                else tmp_buffer_1.push_back(c);
            }
            if (!tmp_buffer_1.empty()) {
                identity=stof(tmp_buffer_1);
                if (nonref_haps.contains(hap_id) && coverage>=MIN_COVERAGE && identity>=MIN_IDENTITY) covered_paths.emplace(hap2path.at(hap_id));
            }
            file.close();
            recall_numerator.at(i).emplace_back(covered_paths.size());
            recall_denominator.at(i).emplace_back(denominator);
        }
    }


    /**
     * Prints every measure as a separate file in `output_dir`. Every file contains a list numbers (one number per
     * window in `windows_to_print`) sorted in non-decreasing order. Only fractional values whose denominator is nonzero
     * are printed.
     *
     * Directory structure:
     *
     * OUTPUT_DIR
     * ├── TOOL_ID_1
     * │   ├── measure_1.txt
     * │   ├── measure_2.txt
     * │   ├── ...
     * ├── TOOL_ID_2
     * │   ├── ...
     * ... ...
     */
    void print_measures(const path& output_dir, const vector<size_t>& windows_to_print) {
        const size_t N_WINDOWS = windows_to_print.size();
        const string SUFFIX = ".txt";
        size_t i, j;
        int32_t truth_tool;
        vector<size_t> tmp_vector_1;
        vector<double> tmp_vector_2;
        vector<vector<ofstream>> out;
        const vector<string> FILE_NAMES = {
                "n_alignments",
                "n_haplotypes",
                "n_nonref_haplotypes",  // Fraction
                "n_haplotype_clusters",
                "haplotype_coverage_avg","nonref_haplotype_coverage_avg","cluster_coverage_avg",
                "alignment_identity_avg","nonref_alignment_identity_avg","cluster_alignment_identity_avg",
                "nonref_nodes",
                "nonref_nodes_fully_covered","nonref_nodes_partially_covered",  // Fractions
                "nonref_nodes_bps",
                "nonref_nodes_bps_covered",  // Fraction
                "nonref_edges",
                "nonref_edges_covered",  // Fraction
                "supported_vcf_records","unsupported_vcf_records","supported_del","unsupported_del","supported_ins","unsupported_ins","supported_inv","unsupported_inv","supported_dup","unsupported_dup",  // Fractions
                "nonref_node_length_vs_coverage","nonref_node_length_vs_identity",  // Pairs
                "sv_length_supported","sv_length_unsupported",
                "nonref_nodes_fully_covered_small","nonref_nodes_fully_covered_large",  // Fractions
                "recall"
        };
        tmp_vector_1.reserve(N_WINDOWS); tmp_vector_2.reserve(N_WINDOWS);

        create_directory(output_dir);
        truth_tool=-1;
        for (i=0; i<N_TOOLS; i++) {
            const string tool_name = TOOLS.at(i);
            if (tool_name==TRUTH_TOOL) truth_tool=i;
            create_directory(output_dir/tool_name);
            out.emplace_back();
            for (j=0; j<FILE_NAMES.size(); j++) {
                out.at(i).emplace_back(output_dir/tool_name/(FILE_NAMES.at(j)+SUFFIX));
                if (!out.at(i).at(j).good() || !out.at(i).at(j).is_open()) throw runtime_error("ERROR: cannot create files in directory "+(output_dir/tool_name).string());
            }
        }
        print_measures_impl(n_alignments,windows_to_print,0,out,false,tmp_vector_1);
        print_measures_impl(n_haplotypes,windows_to_print,1,out,false,tmp_vector_1);
        print_measures_impl_normalized(n_nonref_haplotypes,n_haplotypes,windows_to_print,2,out,true,tmp_vector_2);
        print_measures_impl(n_haplotype_clusters,windows_to_print,3,out,false,tmp_vector_1);
        print_measures_impl(haplotype_coverage_avg,windows_to_print,4,out,false,tmp_vector_2);
        print_measures_impl(nonref_haplotype_coverage_avg,windows_to_print,5,out,false,tmp_vector_2);
        print_measures_impl(cluster_coverage_avg,windows_to_print,6,out,false,tmp_vector_2);
        print_measures_impl(alignment_identity_avg,windows_to_print,7,out,false,tmp_vector_2);
        print_measures_impl(nonref_alignment_identity_avg,windows_to_print,8,out,false,tmp_vector_2);
        print_measures_impl(cluster_alignment_identity_avg,windows_to_print,9,out,false,tmp_vector_2);
        print_measures_impl(nonref_nodes,windows_to_print,10,out,false,tmp_vector_1);
        print_measures_impl_normalized(nonref_nodes_fully_covered,nonref_nodes,windows_to_print,11,out,true,tmp_vector_2);
        print_measures_impl_normalized(nonref_nodes_partially_covered,nonref_nodes,windows_to_print,12,out,true,tmp_vector_2);
        print_measures_impl(nonref_nodes_bps,windows_to_print,13,out,false,tmp_vector_1);
        print_measures_impl_normalized(nonref_nodes_bps_covered,nonref_nodes_bps,windows_to_print,14,out,true,tmp_vector_2);
        print_measures_impl(nonref_edges,windows_to_print,15,out,false,tmp_vector_1);
        print_measures_impl_normalized(nonref_edges_covered,nonref_edges,windows_to_print,16,out,true,tmp_vector_2);
        print_measures_impl_binormalized(supported_vcf_records,supported_vcf_records,unsupported_vcf_records,windows_to_print,17,out,true,tmp_vector_2);
        print_measures_impl_binormalized(unsupported_vcf_records,supported_vcf_records,unsupported_vcf_records,windows_to_print,18,out,true,tmp_vector_2);
        print_measures_impl_binormalized(supported_del,supported_del,unsupported_del,windows_to_print,19,out,true,tmp_vector_2);
        print_measures_impl_binormalized(unsupported_del,supported_del,unsupported_del,windows_to_print,20,out,true,tmp_vector_2);
        print_measures_impl_binormalized(supported_ins,supported_ins,unsupported_ins,windows_to_print,21,out,true,tmp_vector_2);
        print_measures_impl_binormalized(unsupported_ins,supported_ins,unsupported_ins,windows_to_print,22,out,true,tmp_vector_2);
        print_measures_impl_binormalized(supported_inv,supported_inv,unsupported_inv,windows_to_print,23,out,true,tmp_vector_2);
        print_measures_impl_binormalized(unsupported_inv,supported_inv,unsupported_inv,windows_to_print,24,out,true,tmp_vector_2);
        print_measures_impl_binormalized(supported_dup,supported_dup,unsupported_dup,windows_to_print,25,out,true,tmp_vector_2);
        print_measures_impl_binormalized(unsupported_dup,supported_dup,unsupported_dup,windows_to_print,26,out,true,tmp_vector_2);
        print_measures_impl_pair(nonref_node_length_vs_coverage,27,out);
        print_measures_impl_pair(nonref_node_length_vs_coverage,28,out);
        print_measures_impl_basic(sv_length_supported,29,out);
        print_measures_impl_basic(sv_length_unsupported,30,out);
        print_measures_impl_normalized(nonref_nodes_fully_covered_small,nonref_nodes_small,windows_to_print,31,out,true,tmp_vector_2);
        print_measures_impl_normalized(nonref_nodes_fully_covered_large,nonref_nodes_large,windows_to_print,32,out,true,tmp_vector_2);
        if (truth_tool!=-1) {
            print_measures_impl_normalized(truth_tool, recall_numerator, recall_denominator, windows_to_print, 33, out);
        }
    }


    /**
     * Stores a list of anomalous windows for each tool.
     *
     * @param directories all input directories examined.
     */
    void log_anomalous_windows(const path& output_dir, const vector<Region>& directories) {
        const char TSV_SEPARATOR = '\t';
        for (auto& [tool_id, windows]: logged_windows) {
            sort(windows.begin(),windows.end());
            const auto iterator = unique(windows.begin(),windows.end());
            windows.resize(distance(windows.begin(),iterator));
            ofstream out(output_dir/TOOLS.at(tool_id)/LOGGED_WINDOWS_FILE);
            for (auto w: windows) out << directories.at(w).name+TSV_SEPARATOR+to_string(directories.at(w).start)+TSV_SEPARATOR+to_string(directories.at(w).stop) << '\n';
            out.close();
        }
    }


    /**
     * Stores a list of windows with too few haplotypes for each tool.
     *
     * @param directories all input directories examined.
     */
    void log_few_haps_windows(const path& output_dir, const vector<Region>& directories) {
        const char TSV_SEPARATOR = '\t';
        int64_t i;
        int64_t window_id;

        for (auto& [tool_id, windows]: few_haps_windows) {
            ofstream out(output_dir/TOOLS.at(tool_id)/FEW_HAPS_WINDOWS_FILE);
            for (i=0; i<windows.size(); i+=2) {
                window_id=windows.at(i);
                out << directories.at(window_id).name+TSV_SEPARATOR+to_string(directories.at(window_id).start)+TSV_SEPARATOR+to_string(directories.at(window_id).stop) << TSV_SEPARATOR << windows.at(i+1) << '\n';
            }
            out.close();
        }
    }


    /**
     * @return TRUE iff every tool completed evaluation in `directory`, according to the execution log file.
     */
    static bool was_execution_successful(const path& directory) {
        char c;
        size_t field;
        string buffer;
        ifstream file;

        bool success;
        string name;
        size_t h, m, s, ms;

        auto on_field_end = [&]() {
            field++;
            switch (field) {
                case 1: name=buffer; break;
                case 2: h=stoul(buffer); break;
                case 3: m=stoul(buffer); break;
                case 4: s=stoul(buffer); break;
                case 5: ms=stoul(buffer); break;
                case 6: success=buffer=="1"; break;
            }
            buffer.clear();
        };
        file.open(directory/EXECUTION_FILE);
        file.ignore(STREAMSIZE_MAX,LINE_DELIMITER);  // Skipping CSV header
        field=0;
        while (file.get(c)) {
            if (c==LINE_DELIMITER) {
                on_field_end();
                if (!success) return false;
                field=0;
            }
            else if (c==CSV_DELIMITER) on_field_end();
            else buffer.push_back(c);
        }
        if (!buffer.empty()) {
            on_field_end();
            if (!success) return false;
        }
        file.close();
        return true;
    }


private:
    inline static const string NODES_FILE = "nodes.csv";
    inline static const string HAPLOTYPES_FILE = "haps.csv";
    inline static const string CLUSTERS_FILE = "clusters.csv";
    inline static const string EDGE_COUNTS_FILE = "edges.csv";
    inline static const string SUPPORTED_VCF_FILE = "supported.vcf";
    inline static const string UNSUPPORTED_VCF_FILE = "unsupported.vcf";
    inline static const string LOGGED_WINDOWS_FILE = "anomalous_windows.bed";
    inline static const string FEW_HAPS_WINDOWS_FILE = "few_haps_windows.bed";
    inline static const string EXECUTION_FILE = "log.csv";
    inline static const string ALIGNMENTS_FILE = "alignments.gaf";
    static const char GAF_FWD_CHAR = '>';
    static const char GAF_REV_CHAR = '<';
    static const char CSV_DELIMITER = ',';
    static const char TSV_DELIMITER = '\t';
    static const char LINE_DELIMITER = '\n';
    static const char CLUSTER_DELIMITER = ' ';
    static const size_t STREAMSIZE_MAX = numeric_limits<streamsize>::max();

    const vector<string>& TOOLS;
    const string TRUTH_TOOL;
    const size_t N_TOOLS;
    const double coverage_threshold;
    const size_t max_hap_length;
    const double log_nodes_fully_covered_leq, log_edges_covered_leq, log_vcf_records_supported_leq, log_alignment_identity_leq, log_hap_coverage_leq;
    const size_t min_sv_length;
    const size_t N_HAPLOTYPES_EXPECTED;

    size_t n_windows;

    /**
     * Output measures: alignments.
     */
    vector<vector<size_t>> n_alignments;

    /**
     * Output measures: haplotypes.
     */
    vector<vector<size_t>> n_haplotypes, n_nonref_haplotypes;
    vector<size_t> n_haplotype_clusters;

    /**
     * Output measures: haplotype averages. Every haplotype contributes equally.
     */
    vector<vector<double>> haplotype_coverage_avg, alignment_identity_avg;

    /**
     * Output measures: haplotype averages. Every haplotype contributes equally. Haplotypes marked as reference do not
     * contribute.
     */
    vector<vector<double>> nonref_haplotype_coverage_avg, nonref_alignment_identity_avg;

    /**
     * Output measures: haplotype clusters. Assume that every haplotype is assigned to a cluster. The following measures
     * are computed by averaging over all the haplotypes in the same cluster, and then by computing the average over all
     * clusters. Every cluster contributes equally to the measure, regardless of the number of haplotypes it contains.
     */
    vector<vector<double>> cluster_coverage_avg, cluster_alignment_identity_avg;

    /**
     * Output measures: graph.
     */
    vector<vector<size_t>> nonref_nodes, nonref_nodes_small, nonref_nodes_large, nonref_nodes_fully_covered, nonref_nodes_fully_covered_small, nonref_nodes_fully_covered_large, nonref_nodes_partially_covered;
    vector<vector<size_t>> nonref_nodes_bps, nonref_nodes_bps_covered;
    vector<vector<size_t>> nonref_edges, nonref_edges_covered;
    vector<vector<pair<size_t,double>>> nonref_node_length_vs_coverage;
    vector<vector<pair<size_t,double>>> nonref_node_length_vs_identity;

    /**
     * Output measures: VCF.
     */
    vector<vector<size_t>> supported_vcf_records, unsupported_vcf_records;
    vector<vector<size_t>> supported_del, unsupported_del, supported_ins, unsupported_ins, supported_inv, unsupported_inv, supported_dup, unsupported_dup;
    vector<vector<size_t>> sv_length_supported, sv_length_unsupported;

    /**
     * Output measures: recall.
     */
     vector<vector<size_t>> recall_numerator, recall_denominator;

    /**
     * Output log: row=toolID, column=window.
     */
    unordered_map<size_t,vector<size_t>> logged_windows, few_haps_windows;

    /**
     * Reused temporary space: stores the current line of a CSV file.
     */
    bool is_ref, is_flank;
    size_t n_alignments_in_window, n_edges, n_nonref_edges, n_edges_covered, n_nonref_edges_covered;
    int64_t length;
    double coverage, identity;
    string name, names, cluster_id, color;
    string tmp_buffer_1, tmp_buffer_2;

    /**
     * Reused temporary space: cumulates counts over the current file.
     */
    vector<size_t> node_counts;

    size_t n_clusters;
    unordered_map<string,string> hap2cluster;
    unordered_map<string,size_t> cluster_size;
    unordered_map<string,vector<double>> cluster_counts;

    size_t n_haps, n_nonref_haps;
    vector<double> hap_counts, nonref_hap_counts;

    vector<size_t> vcf_counts_supported, vcf_counts_unsupported;

    set<size_t> ref_node_ids;
    set<string> nonref_haps;
    vector<pair<string,bool>> tmp_path;
    unordered_map<string,string> hap2path;
    set<string> covered_paths;
    set<string> canonized_paths;


    /**
     * `node_counts` has the following content:
     * nonref_nodes
     * nonref_nodes_fully_covered
     * nonref_nodes_partially_covered
     * nonref_nodes_bps
     * nonref_nodes_bps_covered
     * nonref_nodes_small
     * nonref_nodes_large
     * nonref_nodes_fully_covered_small
     * nonref_nodes_fully_covered_large
     */
    void on_init_nodes() {
        node_counts={0,0,0,0,0,0,0,0,0};
    }

    void on_field_end_nodes(size_t field, const string& buffer) {
        switch (field) {
            case 1: name=buffer; break;
            case 2: length=stol(buffer); break;
            case 3: is_ref=stoi(buffer)==1; break;
            case 4: is_flank=stoi(buffer)==1; break;
            case 5: coverage=stod(buffer); break;
            case 6: identity=stod(buffer); break;
            case 7: color=buffer; break;
            default: throw runtime_error("ERROR: node CSV field not recognized: "+to_string(field));
        }
    }

    void on_line_end_nodes(size_t tool_id) {
        if (is_ref) return;
        node_counts.at(0)++;
        if (coverage>=coverage_threshold) {
            node_counts.at(1)++;
            if (length<min_sv_length) node_counts.at(7)++;
            else node_counts.at(8)++;
        }
        else node_counts.at(2)++;
        node_counts.at(3)+=length;
        node_counts.at(4)+=(size_t)(coverage*length);
        if (length<min_sv_length) node_counts.at(5)++;
        else node_counts.at(6)++;
        nonref_node_length_vs_coverage.at(tool_id).emplace_back(length,coverage);
        nonref_node_length_vs_identity.at(tool_id).emplace_back(length,identity);
    }

    void on_window_end_nodes(size_t tool_id) {
        nonref_nodes.at(tool_id).emplace_back(node_counts.at(0));
        nonref_nodes_fully_covered.at(tool_id).emplace_back(node_counts.at(1));
        nonref_nodes_partially_covered.at(tool_id).emplace_back(node_counts.at(2));
        nonref_nodes_bps.at(tool_id).emplace_back(node_counts.at(3));
        nonref_nodes_bps_covered.at(tool_id).emplace_back(node_counts.at(4));
        nonref_nodes_small.at(tool_id).emplace_back(node_counts.at(5));
        nonref_nodes_large.at(tool_id).emplace_back(node_counts.at(6));
        nonref_nodes_fully_covered_small.at(tool_id).emplace_back(node_counts.at(7));
        nonref_nodes_fully_covered_large.at(tool_id).emplace_back(node_counts.at(8));

        if (node_counts.at(0)>0 && node_counts.at(1)<=log_nodes_fully_covered_leq*node_counts.at(0)) {
            if (logged_windows.contains(tool_id)) logged_windows.at(tool_id).push_back(nonref_nodes.at(tool_id).size()-1);
            else logged_windows[tool_id]={nonref_nodes.at(tool_id).size()-1};
        }
    }


    void on_init_clusters() {
        hap2cluster.clear(); cluster_size.clear();
    }

    void on_field_end_clusters(size_t field, const string& buffer) {
        switch (field) {
            case 1: cluster_id=buffer; break;
            case 2: names=buffer; break;
            default: throw runtime_error("ERROR: clusters CSV field not recognized: "+to_string(field));
        }
    }

    void on_line_end_clusters(string& buffer) {
        char c;
        size_t i;

        buffer.clear();
        for (i=0; i<names.size(); i++) {
            c=names.at(i);
            if (c==CLUSTER_DELIMITER) { hap2cluster.emplace(buffer,cluster_id); buffer.clear(); }
            else buffer.push_back(c);
        }
        hap2cluster.emplace(buffer,cluster_id);
    }

    void on_window_end_clusters() {
        n_clusters=0;
        for (auto& [haplotype_id,cluster_id]: hap2cluster) {
            if (!cluster_size.contains(cluster_id)) {
                n_clusters++;
                cluster_size[cluster_id]=1;
            }
            else cluster_size.at(cluster_id)++;
        }
        n_haplotype_clusters.emplace_back(n_clusters);
    }


    /**
     * `hap_counts` and `nonref_hap_counts` have the following content:
     * coverage_sum
     * alignment_identity_sum
     */
    void on_init_haplotypes() {
        n_haps=0; n_nonref_haps=0;
        hap_counts={0,0}; nonref_hap_counts={0,0};
        cluster_counts.clear();
    }

    void on_field_end_haplotypes(size_t field, const string& buffer) {
        switch (field) {
            case 1: name=buffer; break;
            case 2: length=stol(buffer); break;
            case 3: is_ref=stoi(buffer)==1; break;
            case 4: is_flank=stoi(buffer)==1; break;
            case 5: coverage=buffer.starts_with("nan")||buffer.ends_with("nan")?0.0:stod(buffer); break;
            case 6: identity=buffer.starts_with("nan")||buffer.ends_with("nan")?0.0:stod(buffer); break;
            default: throw runtime_error("ERROR: haplotypes CSV field not recognized: "+to_string(field));
        }
    }

    void on_line_end_haplotypes(const path& input_file) {
        if (length>max_hap_length) throw runtime_error("ERROR: haplotype anomalously long:"+to_string(length)+"bps, file "+input_file.string());
        if (length<0) throw runtime_error("ERROR: haplotype with negative length:"+to_string(length)+"bps, file "+input_file.string());
        n_haps++;
        hap_counts.at(0)+=coverage;
        hap_counts.at(1)+=identity;
        if (!is_ref) {
            n_nonref_haps++;
            nonref_hap_counts.at(0)+=coverage;
            nonref_hap_counts.at(1)+=identity;
        }
        if (!hap2cluster.contains(name)) return;  // Haplotypes without a cluster do not contribute to clustering stats
        const string id = hap2cluster.at(name);
        if (cluster_counts.contains(id)) {
            cluster_counts.at(id).at(0)+=coverage;
            cluster_counts.at(id).at(1)+=identity;
        }
        else cluster_counts[id]={coverage,identity};
    }

    void on_window_end_haplotypes(size_t tool_id, const path& input_file) {
        n_haplotypes.at(tool_id).emplace_back(n_haps);
        n_nonref_haplotypes.at(tool_id).emplace_back(n_nonref_haps);
        haplotype_coverage_avg.at(tool_id).emplace_back(n_haps>0?hap_counts.at(0)/((double)n_haps):-1);
        alignment_identity_avg.at(tool_id).emplace_back(n_haps>0?hap_counts.at(1)/((double)n_haps):-1);
        nonref_haplotype_coverage_avg.at(tool_id).emplace_back(n_nonref_haps>0?nonref_hap_counts.at(0)/((double)n_nonref_haps):-1);
        nonref_alignment_identity_avg.at(tool_id).emplace_back(n_nonref_haps>0?nonref_hap_counts.at(1)/((double)n_nonref_haps):-1);
        double sum1 = 0; double sum2 = 0;
        for (const auto& [id, size]: cluster_size) {
            sum1+=cluster_counts.at(id).at(0)/((double)size);
            sum2+=cluster_counts.at(id).at(1)/((double)size);
        }
        cluster_coverage_avg.at(tool_id).emplace_back(n_clusters>0?sum1/n_clusters:-1);
        cluster_alignment_identity_avg.at(tool_id).emplace_back(n_clusters>0?sum2/n_clusters:-1);
        if (n_haps>0 && (hap_counts.at(0)<=log_hap_coverage_leq*n_haps || hap_counts.at(1)<=log_alignment_identity_leq*n_haps)) {
            if (logged_windows.contains(tool_id)) logged_windows.at(tool_id).push_back(n_haplotypes.at(tool_id).size()-1);
            else logged_windows[tool_id]={n_haplotypes.at(tool_id).size()-1};
        }
        if (n_haps<N_HAPLOTYPES_EXPECTED) {
            if (few_haps_windows.contains(tool_id)) {
                few_haps_windows.at(tool_id).push_back(n_haplotypes.at(tool_id).size()-1);
                few_haps_windows.at(tool_id).push_back(n_haps);
            }
            else { few_haps_windows[tool_id]={n_haplotypes.at(tool_id).size()-1,n_haps}; }
        }
    }


    void on_field_end_edge_counts(size_t field, const string& buffer) {
        switch (field) {
            case 1: n_alignments_in_window=stoi(buffer); break;
            case 2: n_edges=stoi(buffer); break;
            case 3: n_edges_covered=stoi(buffer); break;
            case 4: n_nonref_edges=stoi(buffer); break;
            case 5: n_nonref_edges_covered=stoi(buffer); break;
            default: throw runtime_error("ERROR: edges CSV field not recognized: "+to_string(field));
        }
    }

    void on_window_end_edge_counts(size_t tool_id) {
        n_alignments.at(tool_id).emplace_back(n_alignments_in_window);
        nonref_edges.at(tool_id).emplace_back(n_nonref_edges);
        nonref_edges_covered.at(tool_id).emplace_back(n_nonref_edges_covered);

        if (n_nonref_edges>0 && n_nonref_edges_covered<=log_edges_covered_leq*n_nonref_edges) {
            if (logged_windows.contains(tool_id)) logged_windows.at(tool_id).push_back(n_alignments.at(tool_id).size()-1);
            else logged_windows[tool_id]={n_alignments.at(tool_id).size()-1};
        }
    }


    /**
     * `vcf_counts_supported` and `vcf_counts_unsupported` have the following content:
     * total_records
     * DEL
     * INS
     * INV
     * DUP
     */
    void on_init_vcf_counts() {
        vcf_counts_supported={0,0,0,0,0};
        vcf_counts_unsupported={0,0,0,0,0};
    }

    void on_window_end_vcf_counts(size_t tool_id) {
        supported_vcf_records.at(tool_id).emplace_back(vcf_counts_supported.at(0));
        supported_del.at(tool_id).emplace_back(vcf_counts_supported.at(1));
        supported_ins.at(tool_id).emplace_back(vcf_counts_supported.at(2));
        supported_inv.at(tool_id).emplace_back(vcf_counts_supported.at(3));
        supported_dup.at(tool_id).emplace_back(vcf_counts_supported.at(4));

        unsupported_vcf_records.at(tool_id).emplace_back(vcf_counts_unsupported.at(0));
        unsupported_del.at(tool_id).emplace_back(vcf_counts_unsupported.at(1));
        unsupported_ins.at(tool_id).emplace_back(vcf_counts_unsupported.at(2));
        unsupported_inv.at(tool_id).emplace_back(vcf_counts_unsupported.at(3));
        unsupported_dup.at(tool_id).emplace_back(vcf_counts_unsupported.at(4));

        const size_t total = vcf_counts_supported.at(0)+vcf_counts_unsupported.at(0);
        if (total>0 && vcf_counts_supported.at(0)<=log_vcf_records_supported_leq*total) {
            if (logged_windows.contains(tool_id)) logged_windows.at(tool_id).push_back(supported_vcf_records.at(tool_id).size()-1);
            else logged_windows[tool_id]={supported_vcf_records.at(tool_id).size()-1};
        }
    }


    /**
     * Prints the same array of values for every tool
     */
    template<class T> void print_measures_impl(const vector<T>& source, const vector<size_t>& windows_to_print, size_t column, vector<vector<ofstream>>& out, bool nonzeros_only, vector<T> tmp_vector) const {
        tmp_vector.clear();
        for (auto& value: windows_to_print) tmp_vector.emplace_back(source.at(value));
        sort(tmp_vector.begin(),tmp_vector.end());
        const size_t n_elements = tmp_vector.size();
        for (size_t i=0; i<N_TOOLS; i++) {
            for (size_t j=0; j<n_elements; j++) {
                if (!nonzeros_only || tmp_vector.at(j)!=0) out.at(i).at(column) << to_string(tmp_vector.at(j)) << '\n';
            }
            out.at(i).at(column).close();
        }
    }


    template<class T> void print_measures_impl(const vector<vector<T>>& source, const vector<size_t>& windows_to_print, size_t column, vector<vector<ofstream>>& out, bool nonzeros_only, vector<T> tmp_vector) const {
        size_t n_elements;
        for (size_t i=0; i<N_TOOLS; i++) {
            tmp_vector.clear();
            for (auto& value: windows_to_print) tmp_vector.emplace_back(source.at(i).at(value));
            sort(tmp_vector.begin(),tmp_vector.end());
            n_elements=tmp_vector.size();
            for (size_t j=0; j<n_elements; j++) {
                if (!nonzeros_only || tmp_vector.at(j)!=0) out.at(i).at(column) << to_string(tmp_vector.at(j)) << '\n';
            }
            out.at(i).at(column).close();
        }
    }


    /**
     * @param nonzero_denom_only FALSE: prints a zero for every fraction with zero denominator.
     */
    template<class T> void print_measures_impl_normalized(const vector<vector<T>>& numerator, const vector<vector<size_t>>& denominator, const vector<size_t>& windows_to_print, size_t column, vector<vector<ofstream>>& out, bool nonzero_denom_only, vector<double> tmp_vector) const {
        size_t n_elements;
        for (size_t i=0; i<N_TOOLS; i++) {
            tmp_vector.clear();
            for (auto& value: windows_to_print) {
                const T denom = denominator.at(i).at(value);
                if (denom!=0 || !nonzero_denom_only) tmp_vector.emplace_back(denom!=0?((double)(numerator.at(i).at(value)))/denom:0);
            }
            sort(tmp_vector.begin(),tmp_vector.end());
            n_elements=tmp_vector.size();
            for (size_t j=0; j<n_elements; j++) out.at(i).at(column) << to_string(tmp_vector.at(j)) << '\n';
            out.at(i).at(column).close();
        }
    }


    void print_measures_impl_normalized(size_t exclude_tool, const vector<vector<size_t>>& numerator, const vector<vector<size_t>>& denominator, const vector<size_t>& windows_to_print, size_t column, vector<vector<ofstream>>& out) const {
        for (size_t i=0; i<N_TOOLS; i++) {
            if (i==exclude_tool) continue;
            size_t num = 0; size_t denom = 0;
            for (auto& value: windows_to_print) {
                num+=numerator.at(i).at(value);
                denom+=denominator.at(i).at(value);
            }
            out.at(i).at(column) << to_string(num) << "," << to_string(denom) << "," << (denom==0?"0":to_string(((double)num)/denom)) << '\n';
            out.at(i).at(column).close();
        }
    }


    /**
     * @param nonzero_denom_only FALSE: prints a zero for every fraction with zero denominator.
     */
    template<class T> void print_measures_impl_binormalized(const vector<vector<T>>& numerator, const vector<vector<size_t>>& denominator1, const vector<vector<size_t>>& denominator2, const vector<size_t>& windows_to_print, size_t column, vector<vector<ofstream>>& out, bool nonzero_denom_only, vector<double> tmp_vector) const {
        size_t n_elements;
        for (size_t i=0; i<N_TOOLS; i++) {
            tmp_vector.clear();
            for (auto& value: windows_to_print) {
                const T denom = denominator1.at(i).at(value)+denominator2.at(i).at(value);
                if (denom!=0 || !nonzero_denom_only) tmp_vector.emplace_back(denom!=0?((double)(numerator.at(i).at(value)))/denom:0);
            }
            sort(tmp_vector.begin(),tmp_vector.end());
            n_elements=tmp_vector.size();
            for (size_t j=0; j<n_elements; j++) out.at(i).at(column) << to_string(tmp_vector.at(j)) << '\n';
            out.at(i).at(column).close();
        }
    }


    void print_measures_impl_pair(const vector<vector<pair<size_t,double>>>& source, size_t column, vector<vector<ofstream>>& out) const {
        for (size_t i=0; i<N_TOOLS; i++) {
            for (auto& measure: source.at(i)) out.at(i).at(column) << to_string(measure.first) << ',' << to_string(measure.second) << '\n';
            out.at(i).at(column).close();
        }
    }


    void print_measures_impl_basic(vector<vector<size_t>>& source, size_t column, vector<vector<ofstream>>& out) const {
        for (size_t i=0; i<N_TOOLS; i++) {
            if (source.at(i).size()>1) sort(source.at(i).begin(),source.at(i).end());
            for (auto& measure: source.at(i)) out.at(i).at(column) << to_string(measure) << '\n';
            out.at(i).at(column).close();
        }
    }

};


/**
 * @param directories window coordinates associated with each call to `load_directory()`, in the order in which the
 * calls were issued; format: [start..stop);
 * @param intervals in the same order along the reference as the calls to `load_directory()`; format: [start..stop);
 * @param min_covered_fraction the function returns only windows that intersect with some element of `intervals`; if
 * this value is nonzero, the procedure returns only windows with at least this fraction of bases covered by
 * elements of `intervals`.
 */
void get_windows_to_print(const vector<Region>& directories, const vector<Region>& intervals, double min_covered_fraction, vector<size_t>& out) {
    const size_t N_DIRECTORIES = directories.size();
    const size_t N_INTERVALS = intervals.size();
    size_t i, j, first_j_for_next_i;
    size_t surface, current_intervals_size;
    int32_t first, last, current_directory_start, current_directory_stop;
    string current_directory_name;
    Region current_directory, current_interval;
    vector<Region> current_intervals;

    out.clear();
    j=0; first_j_for_next_i=UINT64_MAX;
    for (i=0; i<N_DIRECTORIES; i++) {
        current_directory=directories.at(i);
        current_directory_name=current_directory.name;
        current_directory_start=current_directory.start;
        current_directory_stop=current_directory.stop;
        current_intervals.clear();
        while (j<N_INTERVALS) {
            current_interval=intervals.at(j);
            if (first_j_for_next_i==UINT64_MAX && i<N_DIRECTORIES-1 && current_interval.name==directories.at(i+1).name && current_interval.stop>directories.at(i+1).start && current_interval.start<directories.at(i+1).stop) first_j_for_next_i=j;
            if (current_interval.name>current_directory_name || (current_interval.name==current_directory_name && current_interval.start>=current_directory_stop)) break;
            if (current_interval.name<current_directory_name || (current_interval.name==current_directory_name && current_interval.stop<=current_directory_start)) { j++; continue; }
            current_intervals.emplace_back(current_interval);
            j++;
        }
        if (!current_intervals.empty()) {
            surface=0;
            first=current_intervals.at(0).start; last=current_intervals.at(0).stop-1;
            current_intervals_size=current_intervals.size();
            for (j=1; j<current_intervals_size; j++) {
                if (current_intervals.at(j).start>last) {
                    surface+=min(last,current_directory_stop-1)-max(first,current_directory_start)+1;
                    first=current_intervals.at(j).start; last=current_intervals.at(j).stop-1;
                }
                else last=max(last,current_intervals.at(j).stop-1);
            }
            surface+=min(last,current_directory_stop-1)-max(first,current_directory_start)+1;
            if (surface>=min_covered_fraction*(current_directory_stop-current_directory_start)) out.emplace_back(i);
        }
        if (first_j_for_next_i!=UINT64_MAX) j=first_j_for_next_i;
        first_j_for_next_i=UINT64_MAX;
    }
}




int main (int argc, char* argv[]) {
    const string SUBDIR_PREFIX = "chr";
    const string SUBDIR_ALL_WINDOWS = "all_windows";
    const size_t PROGRESS_N_DIRS = 100;  // Arbitrary

    bool all_cluster_files_present;
    char c;
    size_t i;
    size_t length, n_directories;
    int32_t first, last;
    string chromosome, current_dir, buffer;
    vector<size_t> windows_to_print;
    vector<string> unsuccessful_directories;
    vector<Region> directories, intervals;

    // Parsing the input
    CLI::App app{"Evaluates VCFs using assemblies aligned to a variation graph"};
    path INPUT_DIR, OUTPUT_DIR;
    vector<string> TOOLS;
    vector<path> BED_FILES;
    double ALIGNMENT_COVERAGE_THRESHOLD = 0.95;
    double BED_COVERAGE_THRESHOLD = 0.1;
    size_t MAX_HAP_LENGTH = 100000;
    double LOG_NODES_FULLY_COVERED_LEQ = 0.8;
    double LOG_EDGES_COVERED_LEQ = 0.8;
    double LOG_VCF_RECORDS_SUPPORTED_LEQ = 0.8;
    double LOG_ALIGNMENT_IDENTITY_LEQ = 0.8;
    double LOG_HAP_COVERAGE_LEQ = 0.97;
    size_t MIN_SV_LENGTH = 50;
    size_t N_HAPLOTYPES_EXPECTED = 47*2;  // From HPRC Y1
    string TRUTH_TOOL = "";
    app.add_option("--input_dir",INPUT_DIR,"Input directory, with one subdirectory per window.")->required();
    app.add_option("--output_dir",OUTPUT_DIR,"Output directory. Must not already exist.")->required();
    app.add_option("--tools",TOOLS,"List of tools to be evaluated. These names are matched to subdirectories of the input directory.")->expected(1,-1)->required();
    app.add_option("--beds",BED_FILES,"List of BED files to select windows for evaluation. No file = Run the evaluation over all windows. BED files can contain overlapping intervals and might not be sorted.")->expected(1,-1);
    app.add_option("--min_alignment_coverage",ALIGNMENT_COVERAGE_THRESHOLD,"Count a node or haplotype as fully covered iff at least this fraction of it is covered by alignments (default: 0.95).")->capture_default_str();
    app.add_option("--min_bed_coverage",BED_COVERAGE_THRESHOLD,"Use a window for evaluation iff at least this fraction of it is covered by BED intervals. 0=Iff even a single basepair of the window is covered by BED intervals.")->capture_default_str();
    app.add_option("--max_hap_length",MAX_HAP_LENGTH,"Haplotypes longer than this number of basepairs are considered errors and they make the program halt immediately.")->capture_default_str();
    app.add_option("--log_nodes_fully_covered",LOG_NODES_FULLY_COVERED_LEQ,"Stores in a file the coordinates of every input directory whose fraction of nodes fully covered is at most this.")->capture_default_str();
    app.add_option("--log_edges_covered",LOG_EDGES_COVERED_LEQ,"Stores in a file the coordinates of every input directory whose fraction of covered edges is at most this.")->capture_default_str();
    app.add_option("--log_vcf_supported",LOG_VCF_RECORDS_SUPPORTED_LEQ,"Stores in a file the coordinates of every input directory whose fraction of supported VCF records is at most this.")->capture_default_str();
    app.add_option("--log_identity",LOG_ALIGNMENT_IDENTITY_LEQ,"Stores in a file the coordinates of every input directory whose avg. haplotype identity is at most this.")->capture_default_str();
    app.add_option("--log_hap_coverage",LOG_HAP_COVERAGE_LEQ,"Stores in a file the coordinates of every input directory whose avg. haplotype coverage is at most this.")->capture_default_str();
    app.add_option("--truth_id",TRUTH_TOOL,"Assumes that this directory contains the true calls. Used only for computing a haplotype-based recall-like measure.")->capture_default_str();
    app.add_option("--min_sv_length",MIN_SV_LENGTH,"Smallest length for a variant to be considered an SV (default="+to_string(MIN_SV_LENGTH)+").")->capture_default_str();
    app.add_option("--n_haplotypes_expected",N_HAPLOTYPES_EXPECTED,"Expected number of haplotypes (default="+to_string(N_HAPLOTYPES_EXPECTED)+").")->capture_default_str();

    app.parse(argc,argv);
    if (exists(OUTPUT_DIR)) throw runtime_error("ERROR: the output directory already exists: "+OUTPUT_DIR.string());

    INPUT_DIR = std::filesystem::weakly_canonical(INPUT_DIR);
    OUTPUT_DIR = std::filesystem::weakly_canonical(OUTPUT_DIR);

    for (auto& b: BED_FILES){
        b = std::filesystem::weakly_canonical(b);
    }

    Counts counts(TOOLS,TRUTH_TOOL,ALIGNMENT_COVERAGE_THRESHOLD,MAX_HAP_LENGTH,LOG_NODES_FULLY_COVERED_LEQ,LOG_EDGES_COVERED_LEQ,LOG_VCF_RECORDS_SUPPORTED_LEQ,LOG_ALIGNMENT_IDENTITY_LEQ,LOG_HAP_COVERAGE_LEQ,MIN_SV_LENGTH,N_HAPLOTYPES_EXPECTED);

    // Sorting all directories by coordinate
    auto region_comparator = [](const Region& a, const Region& b) {
        if (a.name<b.name) return true;
        else if (a.name>b.name) return false;
        if (a.start<b.start) return true;
        else if (a.start>b.start) return false;
        return a.stop<b.stop;
    };
    for (const auto& entry: directory_iterator(INPUT_DIR)) {
        if (!entry.is_directory()) continue;
        current_dir.clear(); current_dir.append(entry.path().stem().string());
        if (!current_dir.starts_with(SUBDIR_PREFIX)) continue;
        if (!Counts::was_execution_successful(entry.path())) {
            unsuccessful_directories.push_back(current_dir);
            continue;
        }
        buffer=""; length=current_dir.length(); first=-1;
        for (i=0; i<length; i++) {
            c=current_dir.at(i);
            if (c==Counts::SUBDIR_SEPARATOR_1) { chromosome.clear(); chromosome.append(buffer); buffer.clear(); }
            else if (c==Counts::SUBDIR_SEPARATOR_2) { first=stoi(buffer); buffer.clear(); }
            else buffer.push_back(c);
        }
        last=stoi(buffer);
        directories.emplace_back(chromosome,first,last);
    }
    n_directories=directories.size();
    if (n_directories>1) sort(directories.begin(),directories.end(),region_comparator);

    // Collecting counts in canonical order
    all_cluster_files_present=true;
    for (i=0; i<n_directories; i++) {
        current_dir.clear();
        current_dir.append(directories.at(i).name+Counts::SUBDIR_SEPARATOR_1+to_string(directories.at(i).start)+Counts::SUBDIR_SEPARATOR_2+to_string(directories.at(i).stop));
        all_cluster_files_present&=counts.load_window(INPUT_DIR/current_dir);
        if (!TRUTH_TOOL.empty()) counts.haplotype_recall(INPUT_DIR/current_dir,TRUTH_TOOL);
        if ((i+1)%PROGRESS_N_DIRS==0) cerr << "Loaded " << to_string(i+1) << " directories\n";
    }
    cerr << "Loaded " << to_string(n_directories) << " directories\n";
    if (!all_cluster_files_present) cerr << "WARNING: some cluster files are missing.\n";

    // Processing each BED file in canonical order
    create_directory(OUTPUT_DIR);
    for (i=0; i<n_directories; i++) windows_to_print.emplace_back(i);
    counts.print_measures(OUTPUT_DIR/SUBDIR_ALL_WINDOWS,windows_to_print);
    counts.log_anomalous_windows(OUTPUT_DIR/SUBDIR_ALL_WINDOWS,directories);
    counts.log_few_haps_windows(OUTPUT_DIR/SUBDIR_ALL_WINDOWS,directories);
    for (auto& bed_file: BED_FILES) {
        intervals.clear();
        for_region_in_bed_file(bed_file,[&](const Region &r) { intervals.emplace_back(r); });
        if (intervals.size()>1) sort(intervals.begin(),intervals.end(),region_comparator);
        get_windows_to_print(directories,intervals,BED_COVERAGE_THRESHOLD,windows_to_print);
        counts.print_measures(OUTPUT_DIR/bed_file.stem(),windows_to_print);
    }

    // Reporting
    if (!unsuccessful_directories.empty()) {
        cerr << "Skipped " << to_string(unsuccessful_directories.size()) << " directories (out of " << to_string(n_directories) << " total directories) where the evaluation did not complete:\n";
        for (const auto& directory: unsuccessful_directories) cerr << directory << '\n';
    }
}
