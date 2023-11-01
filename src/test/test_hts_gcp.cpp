#include <string>
#include <iostream>
#include <cstdlib>
#include <stdexcept>
#include <chrono>
#include <iomanip>

using std::runtime_error;
using std::to_string;
using std::string;
using std::cerr;
using std::cout;

#include "Authenticator.hpp"
#include "Filesystem.hpp"
#include "misc.hpp"

using sv_merge::GoogleAuthenticator;
using sv_merge::get_current_time;
using sv_merge::run_command;
using ghc::filesystem::path;

#include "htslib/include/htslib/hts.h"
#include "htslib/include/htslib/sam.h"



void read_bam_as_system_command(path bam_path){
    string command = "samtools view " + bam_path.string() + " chr20:40000000-40001000";
    string result;

    run_command(command, result);

    cerr << result.size() << '\n';
}


void read_bam_region(path bam_path, string region="chr20:10000000-10000001"){
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

    hts_itr_t *itr = sam_itr_querys(bam_index, bam_header, region.c_str());

    while (sam_itr_next(bam_file, itr, alignment) >= 0) {
        if (alignment->core.tid < 0) {
            continue;
        }

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


void read_bam(path bam_path){
    samFile* bam_file = nullptr;
    bam_hdr_t* bam_header = nullptr;
    hts_idx_t* bam_index = nullptr;
    bam1_t* alignment = bam_init1();

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

    int i = 0;
    while (sam_read1(bam_file, bam_header, alignment) >= 0){

        // Copy contents of C object into strings
        string read_name = bam_get_qname(alignment);
        string ref_name = bam_header->target_name[alignment->core.tid];

        cerr << read_name << ' ' << ref_name << '\n';

        if (++i == 10){
            break;
        }
    }

    hts_close(bam_file);
    bam_hdr_destroy(bam_header);
    bam_destroy1(alignment);
    hts_idx_destroy(bam_index);
}


void fetch_bam_region(){
    path bam_path_requester_pays = "gs://fc-4310e737-a388-4a10-8c9e-babe06aaf0cf/working/HPRC_PLUS/HG002/assemblies/year1_f1_assembly_v2_genbank/alignment/assembly-to-reference/HG002.maternal.CHM13Y_EBV.bam";
    path bam_path_local = "/home/ryan/data/test_hapslap/results/assorted_prc150_min_n_plus_d_x2_cov1/chr20_18689217-18689256/paths/all_paths_vs_ref.bam";
    path bam_path = "gs://genomics-public-data/platinum-genomes/bam/NA12889_S1.bam";

    GoogleAuthenticator auth;
    auth.update();

    hts_set_log_level(HTS_LOG_TRACE);
    read_bam(bam_path);
}


int main(){
    fetch_bam_region();
}
