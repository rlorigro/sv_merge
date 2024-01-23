#include "Filesystem.hpp"
#include "CLI11.hpp"
#include "misc.hpp"

using ghc::filesystem::path;

#include <stdexcept>

using std::string;
using std::vector;
using std::ifstream;
using std::runtime_error;


using namespace sv_merge;


class Counts {
    const string COUNTS_FILE = "counts.txt";
    const char COUNTS_DELIMITER = '\n';

    vector<string>& tools;
    size_t N_TOOLS;

    vector<vector<int32_t>> n_alignments;

    /**
     * Haplotype averages. Every haplotype contributes equally.
     *
     * A haplotype is classified as reference iff it has just one alignment to the graph containing all calls, and such
     * an alignment uses only reference nodes.
     */
    vector<vector<double>> haplotype_coverage_avg;
    vector<vector<double>> nonref_haplotype_coverage_avg;  // All haps excluding the reference hap
    vector<vector<double>> alignment_identity_avg;
    vector<vector<double>> nonref_alignment_identity_avg;  // All haps excluding the reference hap

    /**
     * Haplotype cluster measures. Assume that every haplotype is assigned to a cluster. The following measures are
     * computed by averaging over all the haplotypes in the same cluster, and then by computing the average over all
     * clusters. Every cluster contributes equally to the measure, regardless of the number of haplotypes it contains.
     *
     * Remark: the cluster ID of a haplotype could be its path in the graph that contains all calls.
     */
    vector<vector<double>> cluster_coverage_avg, cluster_alignment_identity_avg;

    /**
     * Graph measures
     */
    vector<vector<int32_t>> nonref_nodes, nonref_nodes_fully_covered, nonref_nodes_partially_covered;
    vector<vector<int32_t>> nonref_nodes_bps, nonref_nodes_bps_covered;
    vector<vector<int32_t>> nonref_edges, nonref_edges_covered;

    Counts(vector<string>& tools):
        tools(tools)
    {
        size_t i;
        N_TOOLS=tools.size();

        n_alignments.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) n_alignments.emplace_back();
        haplotype_coverage_avg.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) haplotype_coverage_avg.emplace_back();
        nonref_haplotype_coverage_avg.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) nonref_haplotype_coverage_avg.emplace_back();
        alignment_identity_avg.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) alignment_identity_avg.emplace_back();
        nonref_alignment_identity_avg.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) nonref_alignment_identity_avg.emplace_back();

        cluster_coverage_avg.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) cluster_coverage_avg.emplace_back();
        cluster_alignment_identity_avg.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) cluster_alignment_identity_avg.emplace_back();

        nonref_nodes.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) nonref_nodes.emplace_back();
        nonref_nodes_fully_covered.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) nonref_nodes_fully_covered.emplace_back();
        nonref_nodes_partially_covered.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) nonref_nodes_partially_covered.emplace_back();
        nonref_nodes_bps.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) nonref_nodes_bps.emplace_back();
        nonref_nodes_bps_covered.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) nonref_nodes_bps_covered.emplace_back();
        nonref_edges.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) nonref_edges.emplace_back();
        nonref_edges_covered.reserve(N_TOOLS);
        for (i=0; i<N_TOOLS; i++) nonref_edges_covered.emplace_back();
    }


    void load_impl(size_t field, size_t tool_id, string& buffer) {
        switch (field) {
            case 1: n_alignments.at(tool_id).push_back(stoi(buffer)); break;
            case 2: haplotype_coverage_avg.at(tool_id).push_back(stod(buffer)); break;
            case 3: nonref_haplotype_coverage_avg.at(tool_id).push_back(stod(buffer)); break;
            case 4: alignment_identity_avg.at(tool_id).push_back(stod(buffer)); break;
            case 5: nonref_alignment_identity_avg.at(tool_id).push_back(stod(buffer)); break;
            case 6: cluster_coverage_avg.at(tool_id).push_back(stod(buffer)); break;
            case 7: cluster_alignment_identity_avg.at(tool_id).push_back(stod(buffer)); break;
            case 8: nonref_nodes.at(tool_id).push_back(stoi(buffer)); break;
            case 9: nonref_nodes_fully_covered.at(tool_id).push_back(stoi(buffer)); break;
            case 10: nonref_nodes_partially_covered.at(tool_id).push_back(stoi(buffer)); break;
            case 11: nonref_nodes_bps.at(tool_id).push_back(stoi(buffer)); break;
            case 12: nonref_nodes_bps_covered.at(tool_id).push_back(stoi(buffer)); break;
            case 13: nonref_edges.at(tool_id).push_back(stoi(buffer)); break;
            case 14: nonref_edges_covered.at(tool_id).push_back(stoi(buffer)); break;
        }
    }


    /**
     * Adds to the lists of measures all the values in the directory of a window.
     *
     * @param directory contains a subdirectory per tool; each subdirectory contains a file named `COUNTS_FILE`.
     */
    void load(const path& directory) {
        char c;
        size_t field;
        string buffer;

        for (size_t i=0; i<N_TOOLS; i++) {
            ifstream file(directory/COUNTS_FILE);
            if (!file.good() || !file.is_open()) throw runtime_error("ERROR: could not read file: " + (directory/COUNTS_FILE).string());
            field=0;
            while (file.get(c)) {
                if (c!=COUNTS_DELIMITER) { buffer.push_back(c); continue; }
                load_impl(++field,i,buffer);
                buffer.clear();
            }
            if (!buffer.empty()) load_impl(++field,i,buffer);
            file.close();
        }
    }


};










int main (int argc, char* argv[]) {

}