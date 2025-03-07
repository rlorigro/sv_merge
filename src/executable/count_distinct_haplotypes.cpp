#include "VcfReader.hpp"
#include "CLI11.hpp"

using sv_merge::VcfRecord;
using sv_merge::VcfReader;

#include <sstream>
#include <iostream>
#include <filesystem>

using std::filesystem::path;
using std::filesystem::exists;
using std::numeric_limits;
using std::streamsize;
using std::string;
using std::vector;
using std::ifstream;
using std::sort;
using std::to_string;
using std::ofstream;
using std::pair;
using std::cout;
using std::stringstream;


using namespace sv_merge;


/**
 * @param window assumed to contain only phased GTs;
 * @return the number of distinct haplotypes in `window`.
 */
int64_t count_distinct_haplotypes(vector<VcfRecord>& window) {
    const int64_t N_SAMPLES = window.at(0).n_samples;
    const size_t N_RECORDS = window.size();
    int32_t i, j;
    int32_t n_gts;
    pair<int8_t,int8_t> tmp_pair;
    vector<string> haplotypes;

    for (j=0; j<N_SAMPLES; j++) {
        string hap1, hap2;
        for (i=0; i<N_RECORDS; i++) {
            if (!VcfRecord::is_phased(window.at(i).genotypes.at(j))) {
                window.at(i).print(cerr);
                throw runtime_error("ERROR: the "+to_string(i)+"-th sample in the record above has an unphased GT: "+window.at(i).genotypes.at(j));
            }
            n_gts=window.at(i).get_gt(j,tmp_pair);
            if (n_gts!=2) {
                window.at(i).print(cerr);
                throw runtime_error("ERROR: the "+to_string(i)+"-th sample in the record above does not have two values in its GT: "+window.at(i).genotypes.at(j));
            }
            hap1+=to_string(tmp_pair.first); hap2+=to_string(tmp_pair.second);
        }
        haplotypes.push_back(hap1); haplotypes.push_back(hap2);
    }
    sort(haplotypes.begin(),haplotypes.end());
    const auto iterator = unique(haplotypes.begin(),haplotypes.end());
    return distance(haplotypes.begin(),iterator);
}


int main (int argc, char* argv[]) {
    int64_t i;
    int64_t neighborhood_id, nid, n_haps, format_id;
    string buffer, WINDOW_FORMAT_FIELD;
    const char FORMAT_DELIMITER = ':';
    path INPUT_VCF;
    vector<string> format;
    vector<VcfRecord> window;

    // Parsing the input
    CLI::App app{"For every window, prints the number of distinct haplotypes in the cohort."};
    app.add_option("--input_vcf",INPUT_VCF,"Input VCF")->required();
    app.add_option("--window_format_field",WINDOW_FORMAT_FIELD,"Format field that specifies the window ID. We assume that the same windows are used in every sample, but that they might have different IDs in different samples.")->required();
    app.parse(argc,argv);
    INPUT_VCF=std::filesystem::weakly_canonical(INPUT_VCF);

    // Scanning the VCF
    VcfReader reader(INPUT_VCF);
    reader.progress_n_lines=0;
    neighborhood_id=-1;
    reader.for_record_in_vcf([&](VcfRecord& record) {
        // Finding the format field
        stringstream stream1(record.format);
        format_id=-1; i=-1;
        while (!stream1.eof()) {
            getline(stream1,buffer,FORMAT_DELIMITER);
            i++;
            if (buffer==WINDOW_FORMAT_FIELD) { format_id=i; break; }
        }
        if (format_id==-1) throw runtime_error("Cannot find "+WINDOW_FORMAT_FIELD+" among the format fields");

        // Reading the window ID from the first sample
        stringstream stream2(record.genotypes.at(0));
        nid=-1; i=-1;
        while (!stream2.eof()) {
            getline(stream2,buffer,FORMAT_DELIMITER);
            i++;
            if (i<format_id) continue;
            if (buffer==".") {
                // The record does not belong to any window
                if (!window.empty()) {
                    n_haps=count_distinct_haplotypes(window);
                    cout << to_string(n_haps) << '\n';
                    window.clear();
                }
                neighborhood_id=-1;
                return;
            }
            nid=stoi(buffer);
            break;
        }
        if (neighborhood_id==-1) {
            neighborhood_id=nid;
            window.emplace_back(record);
        }
        else if (nid==neighborhood_id) window.emplace_back(record);
        else {
            n_haps=count_distinct_haplotypes(window);
            cout << to_string(n_haps) << '\n';
            neighborhood_id=nid;
            window.clear();
            window.emplace_back(record);
        }
    });
    // Last window
    if (!window.empty()) {
        n_haps=count_distinct_haplotypes(window);
        cout << to_string(n_haps) << '\n';
        window.clear();
    }
}
