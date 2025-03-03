#include "cpptrace/from_current.hpp"
#include "Authenticator.hpp"
#include <filesystem>
#include "Alignment.hpp"
#include "Region.hpp"
#include "fasta.hpp"
#include "fetch.hpp"
#include "misc.hpp"
#include "bam.hpp"

using std::filesystem::path;
using namespace sv_merge;

#include <unordered_map>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <cmath>
#include <span>
#include <map>

using std::unordered_map;
using std::runtime_error;
using std::to_string;
using std::ofstream;
using std::string;
using std::cerr;
using std::span;
using std::map;
using std::max;


void get_reads_of_region(vector<Region>& regions, path csv_path, int32_t flank_length){

    Timer t;
    unordered_map<Region, TransMap> region_transmaps;

    fetch_reads_from_clipped_bam(
        t,
        regions,
        csv_path,
        1,
        9999,
        flank_length,
        region_transmaps,
        true,
        true,
        false,
        false,
        false,
        flank_length
    );

    string s;
    for (auto& region: regions){
        auto& transmap = region_transmaps.at(region);

        cerr << '\n';
        cerr << region.to_string() << '\n';

        size_t n_reads = 0;
        transmap.for_each_read([&](const string& name, const BinarySequence<uint64_t>& sequence){
            sequence.to_string(s);
            cerr << std::left << std::setw(60) << name << s << '\n';
            n_reads++;
        });

        if (n_reads < 14){
            cerr << "FAIL: " << n_reads << " of 14 reads fetch for " << region.to_string() << " with flank_length: " << flank_length << '\n';
        }
    }

}


void test_fetch(path data_directory){
    path bam_path = data_directory / "test_alignment_softclip_only_sorted.bam";
    path csv_path = bam_path;
    csv_path.replace_extension(".csv");

    ofstream output_csv_file(csv_path);

    output_csv_file << "TEST_SAMPLE," << bam_path.string() << '\n';

    output_csv_file.close();

    int32_t window_length = 100;
    int32_t flank_length = 20;

    for (int32_t start=2000; start<=2075; start++){
        vector<Region> regions;
        regions.emplace_back("a", start, start+window_length);
        cerr << regions.back().to_string() << '\n';
        get_reads_of_region(regions, csv_path, flank_length);
    }

    window_length = 51;
    flank_length = 25;

    for (int32_t start=2000; start<=2075; start++){
        vector<Region> regions;
        regions.emplace_back("a", start, start+window_length);
        cerr << regions.back().to_string() << '\n';
        get_reads_of_region(regions, csv_path, flank_length);
    }

    window_length = 50;
    flank_length = 25;

    for (int32_t start=2000; start<=2075; start++){
        vector<Region> regions;
        regions.emplace_back("a", start, start+window_length);
        cerr << regions.back().to_string() << '\n';
        get_reads_of_region(regions, csv_path, flank_length);
    }
}


int main(){
    CPPTRACE_TRY {
        HAPESTRY_DEBUG = true;

        path project_directory = path(__FILE__).parent_path().parent_path().parent_path();
        path data_directory = project_directory / "data";

        cerr << data_directory << '\n';

        cerr << "TESTING: fetch\n\n";
        test_fetch(data_directory);
    }
    CPPTRACE_CATCH(const std::exception& e) {
        std::cerr<<"Exception: "<<e.what()<<std::endl;
        cpptrace::from_current_exception().print_with_snippets();
        throw runtime_error("FAIL: exception caught");
    }

    return 0;
}

