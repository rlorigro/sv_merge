#include "Filesystem.hpp"
#include "Region.hpp"
#include "VcfReader.hpp"
#include "CLI11.hpp"
#include "misc.hpp"

using ghc::filesystem::path;
using ghc::filesystem::directory_iterator;
using sv_merge::VcfRecord;
using sv_merge::VcfReader;
using sv_merge::Region;

#include <stdexcept>

using std::numeric_limits;
using std::streamsize;
using std::string;
using std::vector;
using std::unordered_map;
using std::ifstream;
using std::runtime_error;
using std::sort;


using namespace sv_merge;


class Counts {
public:
    /**
     * @param coverage_threshold fraction above which a node or haplotype is considered fully covered;
     * @param n_windows an estimate on the number of windows that will be processed; used just to allocate space.
     */
    explicit Counts(vector<string>& tools, double coverage_threshold, size_t n_haplotypes, size_t n_nonref_haplotypes, size_t n_haplotype_clusters, const vector<size_t>& haplotype_cluster_size, const unordered_map<string,size_t>& haplotype2cluster, size_t n_windows = 1e6):
        tools(tools),
        N_TOOLS(tools.size()),
        coverage_threshold(coverage_threshold),
        n_haplotypes(n_haplotypes),
        n_nonref_haplotypes(n_nonref_haplotypes),
        n_haplotype_clusters(n_haplotype_clusters),
        haplotype_cluster_size(haplotype_cluster_size),
        haplotype2cluster(haplotype2cluster)
    {
        size_t i;

        n_alignments.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { n_alignments.emplace_back(); n_alignments.at(i).reserve(n_windows); }

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
        nonref_nodes_fully_covered.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) { nonref_nodes_fully_covered.emplace_back(); nonref_nodes_fully_covered.at(i).reserve(n_windows); }
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
    }


    /**
     * Adds to the lists of measures all the values in the directory of a window.
     *
     * @param directory assumed structure:
     * ├── TOOL_ID_1
     * │   ├── `NODES_FILE`: for every node (rows): `id,length,is_reference (0/1),coverage,identity`.
     * │   ├── `HAPLOTYPES_FILE`: for every haplotype(row): `id,length,is_reference (0/1),coverage,identity`.
     * │   ├── `OTHER_COUNTS_FILE`: `total_n_alignments,n_nonreference_edges,n_nonreference_edges_covered`.
     * │   ├── `SUPPORTED_VCF_FILE`: contains only calls that are fully supported by some alignment.
     * │   └── `UNSUPPORTED_VCF_FILE`: contains only calls that are not fully supported by any alignment.
     * ├── TOOL_ID_2
     * │   ├── ...
     * ... ...
     */
    void load_window(const path& directory) {
        char c;
        size_t i, j;
        size_t field;
        string buffer;
        path input_file;
        ifstream file;

        n_windows++;
        for (i=0; i<N_TOOLS; i++) {
            // Node counts
            input_file=directory/tools.at(i)/NODES_FILE;
            file.clear(); file.open(input_file);
            if (!file.good() || !file.is_open()) throw runtime_error("ERROR: could not read file: "+input_file.string());
            on_init_nodes();
            file.ignore(STREAMSIZE_MAX,LINE_DELIMITER);  // Skipping CSV header
            field=0;
            while (file.get(c)) {
                if (c==LINE_DELIMITER) { on_line_end_nodes(); buffer.clear(); }
                else if (c==CSV_DELIMITER) { on_field_end_nodes(++field,buffer); buffer.clear(); }
                else buffer.push_back(c);
            }
            if (!buffer.empty()) { on_field_end_nodes(++field,buffer); on_line_end_nodes(); }
            file.close();
            on_window_end_nodes(i);

            // Haplotype counts
            input_file=directory/tools.at(i)/HAPLOTYPES_FILE;
            file.clear(); file.open(input_file);
            if (!file.good() || !file.is_open()) throw runtime_error("ERROR: could not read file: "+input_file.string());
            on_init_haplotypes();
            file.ignore(STREAMSIZE_MAX,LINE_DELIMITER);  // Skipping CSV header
            field=0;
            while (file.get(c)) {
                if (c==LINE_DELIMITER) { on_line_end_haplotypes(); buffer.clear(); }
                else if (c==CSV_DELIMITER) { on_field_end_haplotypes(++field,buffer); buffer.clear(); }
                else buffer.push_back(c);
            }
            if (!buffer.empty()) { on_field_end_haplotypes(++field,buffer); on_line_end_haplotypes(); }
            file.close();
            on_window_end_haplotypes(i);

            // Remaining counts
            input_file=directory/tools.at(i)/OTHER_COUNTS_FILE;
            file.clear(); file.open(input_file);
            if (!file.good() || !file.is_open()) throw runtime_error("ERROR: could not read file: "+input_file.string());
            file.ignore(STREAMSIZE_MAX,LINE_DELIMITER);  // Skipping CSV header
            field=0;
            while (file.get(c)) {
                if (c==LINE_DELIMITER) { on_window_end_other_counts(i); buffer.clear(); }
                else if (c==CSV_DELIMITER) { on_field_end_other_counts(++field,buffer); buffer.clear(); }
                else buffer.push_back(c);
            }
            if (!buffer.empty()) { on_field_end_other_counts(++field,buffer); on_window_end_other_counts(i); }
            file.close();

            // VCF counts
            input_file=directory/tools.at(i)/SUPPORTED_VCF_FILE;
            for (j=0; j<vcf_counts_supported.size(); j++) vcf_counts_supported.at(j)=0;
            VcfReader reader1(input_file);
            reader1.for_record_in_vcf([&](VcfRecord& record) {
                vcf_counts_supported.at(0)++;
                if (record.sv_type==VcfReader::TYPE_DELETION) vcf_counts_supported.at(1)++;
                else if (record.sv_type==VcfReader::TYPE_INSERTION) vcf_counts_supported.at(2)++;
                else if (record.sv_type==VcfReader::TYPE_INVERSION) vcf_counts_supported.at(3)++;
                else if (record.sv_type==VcfReader::TYPE_DUPLICATION) vcf_counts_supported.at(4)++;
            });
            input_file=directory/tools.at(i)/UNSUPPORTED_VCF_FILE;
            for (j=0; j<vcf_counts_unsupported.size(); j++) vcf_counts_unsupported.at(j)=0;
            VcfReader reader2(input_file);
            reader2.for_record_in_vcf([&](VcfRecord& record) {
                vcf_counts_unsupported.at(0)++;
                if (record.sv_type==VcfReader::TYPE_DELETION) vcf_counts_unsupported.at(1)++;
                else if (record.sv_type==VcfReader::TYPE_INSERTION) vcf_counts_unsupported.at(2)++;
                else if (record.sv_type==VcfReader::TYPE_INVERSION) vcf_counts_unsupported.at(3)++;
                else if (record.sv_type==VcfReader::TYPE_DUPLICATION) vcf_counts_unsupported.at(4)++;
            });
            on_window_end_vcf_counts(i,vcf_counts_supported,vcf_counts_unsupported);
        }
    }


private:
    const string NODES_FILE = "nodes.csv";
    const string HAPLOTYPES_FILE = "haplotypes.csv";
    const string OTHER_COUNTS_FILE = "counts.txt";
    const string SUPPORTED_VCF_FILE = "supported.vcf";
    const string UNSUPPORTED_VCF_FILE = "unsupported.vcf";
    const char CSV_DELIMITER = ',';
    const char LINE_DELIMITER = '\n';
    const int32_t STREAMSIZE_MAX = numeric_limits<streamsize>::max();

    const vector<string>& tools;
    const size_t N_TOOLS;
    const double coverage_threshold;
    const size_t n_haplotypes, n_nonref_haplotypes, n_haplotype_clusters;
    const vector<size_t>& haplotype_cluster_size;
    const unordered_map<string,size_t>& haplotype2cluster;

    size_t n_windows;
    vector<vector<size_t>> n_alignments;

    /**
     * Haplotype averages. Every haplotype contributes equally.
     */
    vector<vector<double>> haplotype_coverage_avg, alignment_identity_avg;

    /**
     * Haplotype averages. Every haplotype contributes equally. Haplotypes marked as reference do not contribute.
     */
    vector<vector<double>> nonref_haplotype_coverage_avg, nonref_alignment_identity_avg;

    /**
     * Haplotype cluster measures. Assume that every haplotype is assigned to a cluster. The following measures are
     * computed by averaging over all the haplotypes in the same cluster, and then by computing the average over all
     * clusters. Every cluster contributes equally to the measure, regardless of the number of haplotypes it contains.
     */
    vector<vector<double>> cluster_coverage_avg, cluster_alignment_identity_avg;

    /**
     * Graph measures
     */
    vector<vector<size_t>> nonref_nodes, nonref_nodes_fully_covered, nonref_nodes_partially_covered;
    vector<vector<size_t>> nonref_nodes_bps, nonref_nodes_bps_covered;
    vector<vector<size_t>> nonref_edges, nonref_edges_covered;

    /**
     * VCF measures
     */
    vector<vector<size_t>> supported_vcf_records, unsupported_vcf_records;
    vector<vector<size_t>> supported_del, unsupported_del, supported_ins, unsupported_ins, supported_inv, unsupported_inv, supported_dup, unsupported_dup;

    /**
     * Reused temporary space, storing the current line of a CSV file.
     */
    string name;
    int32_t length;
    bool is_ref;
    double coverage, identity;
    size_t n_alignments_in_window, n_nonref_edges, n_nonref_edges_covered;

    /**
     * Reused temporary space, for cumulating counts over the current file.
     */
    vector<size_t> node_counts;
    vector<double> haplotype_counts, nonref_haplotype_counts;
    vector<vector<double>> cluster_counts;
    vector<size_t> vcf_counts_supported, vcf_counts_unsupported;


    void on_init_nodes() {
        const size_t size = node_counts.size();
        for (size_t i=0; i<size; i++) node_counts.at(i)=0;
    }


    void on_field_end_nodes(size_t field, const string& buffer) {
        switch (field) {
            case 1: name=buffer; break;
            case 2: length=stoi(buffer); break;
            case 3: is_ref=stoi(buffer)==1; break;
            case 4: coverage=stod(buffer); break;
            case 5: identity=stod(buffer); break;
            default: return;
        }
    }


    /**
     * @param counts nonref_nodes, nonref_nodes_fully_covered, nonref_nodes_partially_covered, nonref_nodes_bps,
     * nonref_nodes_bps_covered;
     */
    void on_line_end_nodes() {
        if (is_ref) return;
        node_counts.at(0)++;
        if (coverage>=coverage_threshold) node_counts.at(1)++;
        else node_counts.at(2)++;
        node_counts.at(3)+=length;
        node_counts.at(4)+=(size_t)(coverage*length);
    }


    /**
     * @param counts nonref_nodes, nonref_nodes_fully_covered, nonref_nodes_partially_covered, nonref_nodes_bps,
     * nonref_nodes_bps_covered;
     */
    void on_window_end_nodes(size_t tool_id) {
        nonref_nodes.at(tool_id).emplace_back(node_counts.at(0));
        nonref_nodes_fully_covered.at(tool_id).emplace_back(node_counts.at(1));
        nonref_nodes_partially_covered.at(tool_id).emplace_back(node_counts.at(2));
        nonref_nodes_bps.at(tool_id).emplace_back(node_counts.at(3));
        nonref_nodes_bps_covered.at(tool_id).emplace_back(node_counts.at(4));
    }


    void on_init_haplotypes() {
        size_t i, j;

        for (i=0; i<haplotype_counts.size(); i++) haplotype_counts.at(i)=0;
        for (i=0; i<nonref_haplotype_counts.size(); i++) nonref_haplotype_counts.at(i)=0;
        for (i=0; i<n_haplotype_clusters; i++) {
            for (j=0; j<cluster_counts.at(j).size(); j++) cluster_counts.at(i).at(j)=0;
        }
    }


    void on_field_end_haplotypes(size_t field, const string& buffer) {
        switch (field) {
            case 1: name=buffer; break;
            case 2: length=stoi(buffer); break;
            case 3: is_ref=stoi(buffer)==1; break;
            case 4: coverage=stod(buffer); break;
            case 5: identity=stod(buffer); break;
            default: return;
        }
    }


    /**
     * @param *_counts coverage_sum, alignment_identity_sum.
     */
    void on_line_end_haplotypes() {
        haplotype_counts.at(0)+=coverage;
        haplotype_counts.at(1)+=identity;
        if (!is_ref) {
            nonref_haplotype_counts.at(0)+=coverage;
            nonref_haplotype_counts.at(1)+=identity;
        }
        const size_t cluster_id = haplotype2cluster.at(name);
        cluster_counts.at(cluster_id).at(0)+=coverage;
        cluster_counts.at(cluster_id).at(1)+=identity;
    }


    /**
     * @param *_counts coverage_sum, alignment_identity_sum.
     */
    void on_window_end_haplotypes(size_t tool_id) {
        const size_t cluster_id = haplotype2cluster.at(name);
        const size_t cluster_size = haplotype_cluster_size.at(cluster_id);

        haplotype_coverage_avg.at(tool_id).emplace_back(haplotype_counts.at(0)/((double)n_haplotypes));
        alignment_identity_avg.at(tool_id).emplace_back(haplotype_counts.at(1)/((double)n_haplotypes));
        if (!is_ref) {
            nonref_haplotype_coverage_avg.at(tool_id).emplace_back(nonref_haplotype_counts.at(0)/((double)n_haplotypes));
            nonref_alignment_identity_avg.at(tool_id).emplace_back(nonref_haplotype_counts.at(1)/((double)n_haplotypes));
        }
        cluster_coverage_avg.at(tool_id).emplace_back(cluster_counts.at(cluster_id).at(0)/((double)cluster_size));
        cluster_alignment_identity_avg.at(tool_id).emplace_back(cluster_counts.at(cluster_id).at(1)/((double)cluster_size));
    }


    void on_field_end_other_counts(size_t field, const string& buffer) {
        switch (field) {
            case 1: n_alignments_in_window=stoi(buffer); break;
            case 2: n_nonref_edges=stoi(buffer); break;
            case 3: n_nonref_edges_covered=stoi(buffer)==1; break;
            default: return;
        }
    }


    void on_window_end_other_counts(size_t tool_id) {
        n_alignments.at(tool_id).emplace_back(n_alignments_in_window);
        nonref_edges.at(tool_id).emplace_back(n_nonref_edges);
        nonref_edges_covered.at(tool_id).emplace_back(n_nonref_edges_covered);
    }


    /**
     * @param vcf_counts_* total_records, del, ins, inv, dup.
     */
    void on_window_end_vcf_counts(size_t tool_id, const vector<size_t>& vcf_counts_supported, const vector<size_t>& vcf_counts_unsupported) {
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
    }
};



class DirectoryName {
public:
    int32_t chromosome, first, last;

    DirectoryName(int32_t chromosome, int32_t first, int32_t last):
        chromosome(chromosome),
        first(first),
        last(last)
    { }

    bool operator==(const DirectoryName& other) const {

    }

    bool operator<(const DirectoryName& other) const {

    }


};






int main (int argc, char* argv[]) {
    const path ROOT_DIR;
    vector<string>& tools;
    vector<double> coverage_threshold;
    const size_t n_haplotypes;
    const size_t n_nonref_haplotypes;
    const size_t n_haplotype_clusters;
    const vector<size_t> haplotype_cluster_size;
    const unordered_map<string,size_t> haplotype2cluster;

    char c;
    size_t i;
    size_t length;
    int32_t first, last;
    string chromosome, buffer, current_dir;

    const char SEPARATOR_1 = '_';
    const char SEPARATOR_2 = '-';


    Counts counts_all(tools,coverage_threshold.at(0),n_haplotypes,n_nonref_haplotypes,n_haplotype_clusters,haplotype_cluster_size,haplotype2cluster);

    vector<DirectoryName> directories;

    // Retrieving and sorting all directory names
    for (const auto& entry: directory_iterator(ROOT_DIR)) {
        if (!entry.is_directory()) continue;
        current_dir.clear(); current_dir.append(entry.path().stem().string());
        if (!current_dir.starts_with("chr")) continue;
        buffer=""; length=current_dir.length();
        for (i=0; i<length; i++) {
            c=current_dir.at(i);
            if (c==SEPARATOR_1) { chromosome.clear(); chromosome.append(buffer); buffer.clear(); }
            else if (c==SEPARATOR_2) { first=stoi(buffer); buffer.clear(); }
            else buffer.push_back(c);
        }
        last=stoi(buffer);
        directories.emplace_back(chromosome,first,last);
    }
    if (directories.size()>1) sort(directories.begin(),directories.end());



}