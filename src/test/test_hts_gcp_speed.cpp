#include "Authenticator.hpp"
#include "Timer.hpp"
#include "CLI11.hpp"
#include "misc.hpp"
#include "bed.hpp"
#include "bam.hpp"

#include <unordered_map>
#include <exception>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <thread>
#include <memory>
#include <limits>

using std::numeric_limits;
using std::unordered_map;
using std::runtime_error;
using std::exception;
using std::ifstream;
using std::thread;
using std::atomic;
using std::mutex;
using std::cerr;
using std::min;
using std::cref;
using std::ref;
using std::cerr;

using namespace sv_merge;


//
void iterate_regions(path bam_path, vector<char*> regions){
    GoogleAuthenticator authenticator;
    authenticator.update();

    samFile* bam_file;
    bam_hdr_t* bam_header;
    hts_idx_t* bam_index;
    bam1_t* alignment;

    bam_file = nullptr;
    bam_index = nullptr;
    alignment = bam_init1();


    if ((bam_file = hts_open(bam_path.string().c_str(), "r")) == nullptr) {
        throw runtime_error("ERROR: Cannot open bam file: " + bam_path.string());
    }

    // bam index
    if ((bam_index = sam_index_load(bam_file, bam_path.string().c_str())) == nullptr) {
        throw runtime_error("ERROR: Cannot open index for bam file: " + bam_path.string() + "\n");
    }

    // bam header
    if ((bam_header = sam_hdr_read(bam_file)) == nullptr){
        throw runtime_error("ERROR: Cannot open header for bam file: " + bam_path.string() + "\n");
    }

    hts_name2id_f getid;
    int r_count;

    hts_reglist_t* reglist = hts_reglist_create(regions.data(), int(regions.size()), &r_count, bam_header,  getid);

    hts_itr_t *itr = sam_itr_regions(bam_index, bam_header, reglist, r_count);

    while (sam_itr_next(bam_file, itr, alignment) >= 0) {
        if (alignment->core.tid < 0) {
            continue;
        }

        // Placeholder code to demonstrate working properly
        string read_name = bam_get_qname(alignment);
        string ref_name = bam_header->target_name[alignment->core.tid];
        auto pos = alignment->core.pos;
        auto len = bam_cigar2rlen(int(alignment->core.n_cigar), bam_get_cigar(alignment));

        cerr << ref_name << ' ' << pos << ' ' << pos + len << ' ' << read_name << '\n';
    }

    bam_destroy1(alignment);
    hts_itr_destroy(itr);
    bam_hdr_destroy(bam_header);
    hts_close(bam_file);
    hts_idx_destroy(bam_index);

}


void test(){
    path bam_path = "gs://fc-bb12940d-f7ba-4515-8a98-42de84f63c34/HPRC/HG002/HG002.bam";

    // Construct some arbitrary regions that avoid the telomeres and centromeres of chr20
    vector<string> region_strings;
    vector<char*> regions;

    int64_t left_arm_start = 3'000'000;
    int64_t left_arm_stop = 30'000'000;
    int64_t right_arm_start = 40'000'000;
    int64_t right_arm_stop = 64'000'000;

    int64_t interval_length = 1000;
    int64_t n_intervals = 1000;

    int64_t gap_length = (left_arm_stop-left_arm_start)/n_intervals;
    Region r("chr20", left_arm_start, left_arm_start+interval_length);

    // Fill in left arm regions
    for (int i=0; i<n_intervals; i++){
        r.start += gap_length;
        r.stop += gap_length;

        region_strings.push_back(r.to_string());
        regions.push_back(region_strings.back().data());
    }

    gap_length = (right_arm_stop-right_arm_start)/n_intervals;
    r.start = right_arm_start;
    r.stop = right_arm_start+interval_length;

    // Fill in right arm regions
    for (int i=0; i<n_intervals; i++){
        r.start += gap_length;
        r.stop += gap_length;

        region_strings.push_back(r.to_string());
        regions.push_back(region_strings.back().data());
    }

    // Test working
    for (int64_t i=0; i<regions.size(); i++){
        cerr << regions[i] << '\n';
    }

    iterate_regions(bam_path, regions);
}


int main(){
    test();

    return 0;
}
