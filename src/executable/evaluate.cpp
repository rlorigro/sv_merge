#include "TransitiveMap.hpp"
#include "Timer.hpp"
#include "CLI11.hpp"
#include "fetch.hpp"
#include "misc.hpp"
#include "bed.hpp"

#include <unordered_map>
#include <exception>
#include <stdexcept>
#include <iostream>
#include <fstream>
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
using std::max;
using std::min;
using std::cref;
using std::ref;


using namespace sv_merge;


// TODO: replace this with proper implementation when it exists
void generate_gfas_with_script(const path& vcf, const path& ref_fasta, const path& bed, const path& output_dir){
    path vcf_to_gfa_script = "/home/ryan/code/hapslap/scripts/vcf_to_gfa.py";

    string command = "python3 " + vcf_to_gfa_script.string() +
                     " --vcf " + vcf.string() +
                     " --ref " + ref_fasta.string() +
                     " --output_dir " + output_dir.string() +
                     " --bed " + bed.string();

    run_command(command, true);
}


void evaluate(
        vector<path>& vcfs,
        path output_dir,
        path windows_bed,               // Override the interval graph if this is provided
        path bam_csv,
        path ref_fasta,
        int32_t flank_length,
        int32_t n_threads,
        bool debug
){

    if (ghc::filesystem::exists(output_dir)){
        throw runtime_error("ERROR: output dir exists already: " + output_dir.string());
    }
    else{
        ghc::filesystem::create_directories(output_dir);
    }

    Timer t;

    cerr << t << "Initializing" << '\n';

    vector <pair <string,path> > bam_paths;
    vector<Region> regions;

    cerr << t << "Constructing windows" << '\n';

    load_windows_from_bed(windows_bed, regions);

    for (auto& r: regions) {
        r.start -= flank_length;
        r.stop += flank_length;
    }

    // The container to store all fetched reads and their relationships to samples/paths
    unordered_map<Region,TransMap> region_transmaps;

    cerr << t << "Fetching reads for all windows" << '\n';

    fetch_reads_from_clipped_bam(t, regions, bam_csv, n_threads, true, region_transmaps);

    cerr << t << "Processing windows" << '\n';

    // Generate GFAs and folder structure for every VCF * every region
    // - GFAs should be stored for easy access/debugging later
    // - GAFs should also be stored
    // - FASTA files for extracted truth haplotypes should be optional
    // By default, all of these files will be stored in /dev/shm and then copied into the output dir as a final step.
    // TODO: create option to NOT use /dev/shm/ as staging dir
    // Absolutely must delete the /dev/shm/ copy before terminating
    //
    path staging_dir = "/dev/shm/hapestry_staging_area/";
    for (const auto& vcf: vcfs){
        generate_gfas_with_script(vcf, ref_fasta, windows_bed, staging_dir);
    }

    ghc::filesystem::remove_all(staging_dir);

    // Iterate sample reads
    for (const auto& region: regions){
        // Per-region output and logging directory
        path output_subdir = output_dir / region.to_string('_');

        if (not ghc::filesystem::exists(output_subdir)){
            ghc::filesystem::create_directories(output_subdir);
        }

        // Get the sample-read-path transitive mapping for this region
        auto& transmap = region_transmaps.at(region);

        // Iterate samples within this region
        transmap.for_each_sample([&](const string& sample_name, int64_t sample_id) {


            if (debug) {
                path p = output_subdir / (sample_name + "_extracted_reads.fasta");
                ofstream file(p);

                transmap.for_each_read_of_sample(sample_name, [&](const string &name, int64_t id) {
                    file << '>' << name << '\n';
                    file << transmap.get_sequence(id) << '\n';
                });
            }
        });
    }

    cerr << t << "Done" << '\n';
}


int main (int argc, char* argv[]){
    path output_dir;
    path windows_bed;
    path bam_csv;
    path ref_fasta;
    vector<path> vcfs;
    int32_t flank_length;
    int32_t n_threads = 1;
    bool debug = false;

    CLI::App app{"App description"};

    app.add_option(
            "--n_threads",
            n_threads,
            "Maximum number of threads to use");

    app.add_option(
            "--output_dir",
            output_dir,
            "Path to output directory which must not exist yet")
            ->required();

    app.add_option(
            "--windows",
            windows_bed,
            "Path to BED file containing windows to merge (inferred automatically if not provided)");

    app.add_option(
            "--vcfs",
            bam_csv,
            "List of VCFs to evaluate (space-separated)")
            ->required();

    app.add_option(
            "--bam_csv",
            bam_csv,
            "Simple headerless CSV file with the format [sample_name],[hap_name],[bam_path]")
            ->required();

    app.add_option(
            "--ref",
            ref_fasta,
            "Path to reference sequence FASTA file that corresponds to VCF")
            ->required();

    app.add_option(
            "--flank_length",
            flank_length,
            "How much flanking sequence to use when fetching and aligning reads")
            ->required();

    app.add_flag("--debug", debug, "Invoke this to add more logging and output");

    CLI11_PARSE(app, argc, argv);

    evaluate(vcfs, output_dir, windows_bed, bam_csv, ref_fasta, flank_length, n_threads, debug);

    return 0;
}
