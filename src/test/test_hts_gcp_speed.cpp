#include "Authenticator.hpp"
#include "Timer.hpp"
#include "CLI11.hpp"
#include "misc.hpp"
#include "bed.hpp"
#include "bam.hpp"
#include "hts.h"

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
using std::ofstream;
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
void iterate_region(path bam_path, string region){
    Authenticator authenticator;
    authenticator.update();

    samFile* bam_file;
    bam_hdr_t* bam_header;
    hts_idx_t* bam_index;
    bam1_t* alignment;

    bam_file = nullptr;
    bam_index = nullptr;
    alignment = bam_init1();

    Timer t;
    cerr << t << "Begin" << '\n';
    if ((bam_file = hts_open(bam_path.string().c_str(), "r")) == nullptr) {
        throw runtime_error("ERROR: Cannot open bam file: " + bam_path.string());
    }
    cerr << t << "hts_open" << '\n';

    // bam index
    if ((bam_index = sam_index_load(bam_file, bam_path.string().c_str())) == nullptr) {
        throw runtime_error("ERROR: Cannot open index for bam file: " + bam_path.string() + "\n");
    }
    cerr << t << "sam_index_load" << '\n';

    // bam header
    if ((bam_header = sam_hdr_read(bam_file)) == nullptr){
        throw runtime_error("ERROR: Cannot open header for bam file: " + bam_path.string() + "\n");
    }
    cerr << t << "sam_hdr_read" << '\n';

    hts_name2id_f getid = (hts_name2id_f)bam_name2id;
    int r_count;

    hts_itr_t *itr = sam_itr_querys(bam_index, bam_header, region.c_str());

    cerr << t << "sam_itr_region" << '\n';

    while (sam_itr_next(bam_file, itr, alignment) >= 0) {
        if (alignment->core.tid < 0) {
            continue;
        }

        // Placeholder code to demonstrate working properly
        string read_name = bam_get_qname(alignment);
        string ref_name = bam_header->target_name[alignment->core.tid];
        auto pos = alignment->core.pos;
        auto len = bam_cigar2rlen(int(alignment->core.n_cigar), bam_get_cigar(alignment));

//        cerr << ref_name << ' ' << pos << ' ' << pos + len << ' ' << read_name << '\n';
    }

    cerr << t << "Done iterating" << '\n';

    bam_destroy1(alignment);
    hts_itr_destroy(itr);
    bam_hdr_destroy(bam_header);
    hts_close(bam_file);
    hts_idx_destroy(bam_index);
}


//
void iterate_regions(path bam_path, vector<char*>& regions){
    Authenticator authenticator;
    authenticator.update();

    samFile* bam_file;
    bam_hdr_t* bam_header;
    hts_idx_t* bam_index;
    bam1_t* alignment;

    bam_file = nullptr;
    bam_index = nullptr;
    alignment = bam_init1();

    Timer t;
    cerr << t << "Begin" << '\n';
    if ((bam_file = hts_open(bam_path.string().c_str(), "r")) == nullptr) {
        throw runtime_error("ERROR: Cannot open bam file: " + bam_path.string());
    }
    cerr << t << "hts_open" << '\n';

    // bam index
    if ((bam_index = sam_index_load(bam_file, bam_path.string().c_str())) == nullptr) {
        throw runtime_error("ERROR: Cannot open index for bam file: " + bam_path.string() + "\n");
    }
    cerr << t << "sam_index_load" << '\n';

    // bam header
    if ((bam_header = sam_hdr_read(bam_file)) == nullptr){
        throw runtime_error("ERROR: Cannot open header for bam file: " + bam_path.string() + "\n");
    }
    cerr << t << "sam_hdr_read" << '\n';

    hts_name2id_f getid = (hts_name2id_f)bam_name2id;
    int r_count;

    hts_reglist_t* reglist = hts_reglist_create(regions.data(), int(regions.size()), &r_count, bam_header,  getid);

    cerr << t << "hts_reglist_create" << '\n';

    hts_itr_t *itr = sam_itr_regions(bam_index, bam_header, reglist, r_count);

    cerr << t << "sam_itr_regions" << '\n';

    while (sam_itr_next(bam_file, itr, alignment) >= 0) {
        if (alignment->core.tid < 0) {
            continue;
        }

        // Placeholder code to demonstrate working properly
        string read_name = bam_get_qname(alignment);
        string ref_name = bam_header->target_name[alignment->core.tid];
        auto pos = alignment->core.pos;
        auto len = bam_cigar2rlen(int(alignment->core.n_cigar), bam_get_cigar(alignment));

//        cerr << ref_name << ' ' << pos << ' ' << pos + len << ' ' << read_name << '\n';
    }

    cerr << t << "Done iterating" << '\n';

    bam_destroy1(alignment);
    hts_itr_destroy(itr);
    bam_hdr_destroy(bam_header);
    hts_close(bam_file);
    hts_idx_destroy(bam_index);
}


//
void iterate_regions_manually(path bam_path, vector<char*>& regions){
    Authenticator authenticator(true);
    authenticator.update();

    samFile* bam_file;
    bam_hdr_t* bam_header;
    hts_idx_t* bam_index;
    bam1_t* alignment;

    bam_file = nullptr;
    bam_index = nullptr;
    alignment = bam_init1();

    Timer t;
    cerr << t << "Begin" << '\n';

    if ((bam_file = hts_open(bam_path.string().c_str(), "r")) == nullptr) {
        throw runtime_error("ERROR: Cannot open bam file: " + bam_path.string());
    }
    cerr << t << "hts_open" << '\n';

    // bam index
    if ((bam_index = sam_index_load(bam_file, bam_path.string().c_str())) == nullptr) {
        throw runtime_error("ERROR: Cannot open index for bam file: " + bam_path.string() + "\n");
    }
    cerr << t << "sam_index_load" << '\n';

    // bam header
    if ((bam_header = sam_hdr_read(bam_file)) == nullptr){
        throw runtime_error("ERROR: Cannot open header for bam file: " + bam_path.string() + "\n");
    }
    cerr << t << "sam_hdr_read" << '\n';
    hts_itr_t* itr;

    for (auto region: regions) {
        itr = sam_itr_querys(bam_index, bam_header, region);

        while (sam_itr_next(bam_file, itr, alignment) >= 0) {
            if (alignment->core.tid < 0) {
                continue;
            }

            string read_name = bam_get_qname(alignment);
            string ref_name = bam_header->target_name[alignment->core.tid];
            auto pos = alignment->core.pos;
            auto len = bam_cigar2rlen(int(alignment->core.n_cigar), bam_get_cigar(alignment));

//            cerr << ref_name << ' ' << pos << ' ' << pos + len << ' ' << read_name << '\n';
        }

        hts_itr_destroy(itr);
    }

    cerr << t << "Done iterating" << '\n';

    bam_destroy1(alignment);
    bam_hdr_destroy(bam_header);
    hts_close(bam_file);
    hts_idx_destroy(bam_index);
}


void test(int type) {
    path bam_path = "gs://fc-bb12940d-f7ba-4515-8a98-42de84f63c34/HPRC/HG002/HG002.bam";

    // Construct some arbitrary regions that avoid the telomeres and centromeres of chr20
    vector<string> region_strings;
    vector<char *> regions;

    int64_t left_arm_start = 3'000'000;
    int64_t left_arm_stop = 30'000'000;
    int64_t right_arm_start = 40'000'000;
    int64_t right_arm_stop = 64'000'000;

    int64_t interval_length = 1000;
    int64_t n_intervals = 1000;

    int64_t gap_length = (left_arm_stop - left_arm_start) / n_intervals;
    Region r("chr20", left_arm_start, left_arm_start + interval_length);

    ofstream bed_file("test_regions_chr20.bed");

    // Fill in left arm regions
    for (int i = 0; i < n_intervals; i++) {
        r.start += gap_length;
        r.stop += gap_length;

        bed_file << r.name << '\t' << r.start << '\t' << r.stop << '\n';

        region_strings.push_back(r.to_string());
        regions.push_back(region_strings.back().data());
    }

    gap_length = (right_arm_stop - right_arm_start) / n_intervals;
    r.start = right_arm_start;
    r.stop = right_arm_start + interval_length;

    // Fill in right arm regions
    for (int i = 0; i < n_intervals; i++) {
        r.start += gap_length;
        r.stop += gap_length;

        bed_file << r.name << '\t' << r.start << '\t' << r.stop << '\n';

        region_strings.push_back(r.to_string());
        regions.push_back(region_strings.back().data());
    }

    if (type == 0) {
        cerr << "Iterating chr20" << '\n';
        iterate_region(bam_path, "chr20");
    }
    else if (type == 1) {
        cerr << "Iterating 2000 regions with hts_iter_regions" << '\n';
        iterate_regions(bam_path, regions);
    }
    else if (type == 2) {
        cerr << "Iterating 2000 regions with hts_iter_querys loop" << '\n';
        iterate_regions_manually(bam_path, regions);
    }
    else{
        throw runtime_error("ERROR: cannot parse user provided type index");
    }
}


int main (int argc, char* argv[]){
    int type;

    CLI::App app{"App description"};

    app.add_option(
            "--type",
            type,
            "Maximum number of threads to use")
            ->required();

    CLI11_PARSE(app, argc, argv);

    test(type);

    return 0;
}

