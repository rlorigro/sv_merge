#include <string>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <array>
#include <stdexcept>
#include <chrono>
#include <iomanip>
#include <cstdlib>

using std::runtime_error;
using std::chrono::system_clock;
using std::chrono::seconds;
using std::chrono::hours;
using std::to_string;
using std::string;
using std::cerr;
using std::cout;

#include "Filesystem.hpp"
#include "misc.hpp"

using hapslap::run_command;
using hapslap::get_current_time;
using ghc::filesystem::path;

#include "htslib/include/htslib/hts.h"
#include "htslib/include/htslib/sam.h"


///
/// Consider using this: https://stackoverflow.com/a/57282480
/// To check for expiration
///
class GoogleAuthenticator{
    /// Attributes
    system_clock::time_point expiration = get_current_time() - hours(128);

    // 60min duration by default
    system_clock::duration token_lifetime = seconds(3600);

public:
    /// Methods
    void update();
};


void GoogleAuthenticator::update() {
    // Check if expires within 30sec
    if (expiration - seconds(30) < get_current_time()){
        cerr << "Updating token" << '\n';
        string command = "gcloud auth print-access-token";
        string token;

        // Will throw error if fails, otherwise returns token
        auto t = get_current_time();
        run_command(command, token);
        expiration = t + token_lifetime;

        cerr << "token: " << token << '\n';

        const char* name = "GCS_OAUTH_TOKEN";
        auto env = getenv(name);

        cerr << "env: " << env << '\n';

        // Update the environment with the token
        const char* value = token.data();
        auto result = setenv("GCS_OAUTH_TOKEN", value, 1);
        cerr << "result: " << result << '\n';

        env = getenv(name);
        cerr << "env: " << env << '\n';
    }
}


void fetch_bam_region(){
    path bam_path = "gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/ea7e459a-5db9-4236-866a-f01f2c524357/minimap2/c867d241-8cc7-4fd9-a89d-f67fbba760f3/call-alignAndSortBAM/cacheCopy/HG002.bam";

    // TODO: figure out why this doesn't work... environment variable appears to be set locally, and yet samtools fails
    GoogleAuthenticator auth;
    auth.update();

    const char* name = "GCS_OAUTH_TOKEN";
    auto env = getenv(name);
    cerr << "env: " << env << '\n';

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

    while (sam_read1(bam_file, bam_header, alignment) >= 0){
        string read_name = bam_get_qname(alignment);
        string ref_name = bam_header->target_name[alignment->core.tid];

        cerr << read_name << ' ' << ref_name << '\n';
    }

    hts_close(bam_file);
    bam_hdr_destroy(bam_header);
    bam_destroy1(alignment);
    hts_idx_destroy(bam_index);
}


int main(){
    fetch_bam_region();
}
