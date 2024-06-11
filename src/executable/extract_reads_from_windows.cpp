#include "interval_tree.hpp"
#include "windows.hpp"
#include "fetch.hpp"
#include "Timer.hpp"
#include "CLI11.hpp"
#include "bed.hpp"

using lib_interval_tree::interval_tree_t;

#include <unordered_map>
#include <exception>
#include <stdexcept>
#include <iostream>
#include <thread>
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
using std::max;
using std::cref;
using std::ref;


using namespace sv_merge;


void extract(
        path output_dir,
        path bam_path,
        path bed_path,
        int32_t flank_length,
        bool require_spanning,
        bool force_forward,
        const vector<string>& tags_to_fetch,
        size_t n_threads
        ){

    if (std::filesystem::exists(output_dir)){
        throw runtime_error("ERROR: output dir exists already: " + output_dir.string());
    }
    else{
        std::filesystem::create_directories(output_dir);
    }

    Timer t;

    vector<Region> regions;

    cerr << t << "Reading BED file" << '\n';
    load_windows_from_bed(bed_path, regions);

    // This is only used while loading VCFs to find where each record belongs
    unordered_map <string, interval_tree_t<int32_t> > contig_interval_trees;

    cerr << t << "Flanking windows and writing BED" << '\n';

    // Add flanks, place the regions in the interval tree, and log the windows
    for (auto& r: regions) {
        r.start = max(1,r.start-flank_length);
        r.stop += flank_length;

        contig_interval_trees[r.name].insert({r.start, r.stop});
    }

    cerr << t << "Fetching reads for all windows" << '\n';

    GoogleAuthenticator authenticator;

    // Intermediate object to store results of multithreaded sample read fetching
    sample_region_flanked_coord_map_t sample_to_region_coords;

    // Intermediate objects
    vector <pair <string, path> > sample_bams;
    TransMap sample_only_transmap;

    cerr << t << "Loading CSV" << '\n';

    string sample_name = "sample";

    sample_region_read_map_t sample_to_region_reads;

    sample_only_transmap.add_sample(sample_name);
    sample_bams.emplace_back(sample_name, bam_path);

    // Initialize every combo of sample,region with an empty vector
    for (const auto& region: regions){
        sample_to_region_reads[sample_name][region] = {};
    }

    cerr << t << "Processing windows" << '\n';

    // Thread-related variables
    atomic<size_t> job_index = 0;
    vector<thread> threads;

    threads.reserve(n_threads);

    bool get_qualities = true;

    // Launch threads
    for (uint64_t n=0; n<n_threads; n++){
        try {
            cerr << "launching: " << n << '\n';
            threads.emplace_back(extract_subsequences_from_sample_thread_fn,
                    std::ref(authenticator),
                    std::ref(sample_to_region_reads),
                    std::cref(sample_bams),
                    std::cref(regions),
                    std::ref(require_spanning),
                    force_forward,
                    get_qualities,
                    std::ref(job_index),
                    std::cref(tags_to_fetch),
                    true
            );
        } catch (const exception& e) {
            throw e;
        }
    }

    // Wait for threads to finish
    for (auto& n: threads){
        n.join();
    }
    cerr << t << "Writing reads to fastq" << '\n';

    for (const auto& sample_reads: sample_to_region_reads){
        for (const auto& [region, reads]: sample_reads.second){
            path fastq_path = output_dir / (region.to_string('_') + ".fastq");

            ofstream fastq_file(fastq_path);

            if ((not fastq_file.is_open()) or (not fastq_file.good())){
                throw runtime_error("ERROR: could not write to file: " + fastq_path.string());
            }

            for (const auto& read: reads){
                fastq_file << "@" << read.name << ' ' << (read.is_reverse ? 'R' : 'F') << (read.tags.empty() ? "" : " ") << read.tags << '\n';
                fastq_file << read.sequence << '\n';
                fastq_file << "+" << '\n';
                for (const auto& q: read.qualities){
                    if (q+33 < 33 or q+33 > 126){
                        throw runtime_error("ERROR: quality score out of range: " + std::to_string(q+33) + " for read: " + read.name + " in region: " + region.to_string() + " at position: " + std::to_string(q));
                    }

                    fastq_file << char(q+33);
                }
                fastq_file << '\n';
            }
        }
    }

    cerr << t << "Done" << '\n';
}


/// Function to parse comma separated string as vector<string>
void parse_comma_separated_string(const string& s, vector<string>& result){
    if (s.empty()){
        return;
    }

    size_t start = 0;
    size_t end = s.find(',');

    while (end != string::npos){
        result.push_back(s.substr(start, end-start));
        start = end + 1;
        end = s.find(',', start);
    }

    result.push_back(s.substr(start, end));
}


int main (int argc, char* argv[]){
    path output_dir;
    path windows_bed;
    path bam_path;
    path bed_path;
    path ref;
    int32_t flank_length = 200;
    size_t n_threads = 1;
    bool require_spanning = false;
    bool force_forward = false;
    vector<string> tags_to_fetch;
    string tags_arg;

    CLI::App app{"App description"};

    app.add_option(
            "--output_dir",
            output_dir,
            "Path to output directory which must not exist yet")
            ->required();

    app.add_option(
            "--bam",
            bam_path,
            "Path to BAM file containing reads to be extracted")
            ->required();

    app.add_option(
            "--bed",
            bed_path,
            "Path to BED file containing windows to extract reads from")
            ->required();

    app.add_option(
            "--flank_length",
            flank_length,
            "How much flanking sequence to use when fetching and aligning reads")
            ->required();

    app.add_option(
            "--n_threads",
            n_threads,
            "Maximum number of threads to use for fetching reads");

    app.add_flag(
            "--require_spanning",
            require_spanning,
            "If this flag is invoked, then only reads that span the entire window will be fetched");

    app.add_flag(
            "--force_forward",
            force_forward,
            "If this flag is invoked, reverse complement any reads that are on the reverse strand");

    app.add_option(
            "--tags",
            tags_arg,
            "A comma separated list of tags to fetch from the BAM file (e.g. NM,PS,HP) and append to the "
            "fastq name as space-separated fields");

    CLI11_PARSE(app, argc, argv);

    parse_comma_separated_string(tags_arg, tags_to_fetch);

    extract(
        output_dir,
        bam_path,
        bed_path,
        flank_length,
        require_spanning,
        force_forward,
        tags_to_fetch,
        n_threads
    );

    return 0;
}
