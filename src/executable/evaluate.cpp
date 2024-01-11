#include "TransitiveMap.hpp"
#include "Filesystem.hpp"
#include "Timer.hpp"
#include "CLI11.hpp"
#include "fetch.hpp"
#include "misc.hpp"
#include "bed.hpp"

using ghc::filesystem::create_directories;

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


void write_region_subsequences_to_file_thread_fn(
        const unordered_map<Region,TransMap>& region_transmaps,
        const vector<Region>& regions,
        const path& output_dir,
        const path& filename,
        atomic<size_t>& job_index
){
    size_t i = job_index.fetch_add(1);

    while (i < regions.size()){
        const auto& region = regions.at(i);
        const auto& t = region_transmaps.at(region);

        path output_subdir = output_dir / region.to_string('_');

        create_directories(output_subdir);

        path output_fasta = output_subdir / filename;
        ofstream file(output_fasta);

        t.for_each_read([&](const string& name, int64_t id){
            file << '>' << name << '\n';
            file << t.get_sequence(id) << '\n';
        });

        i = job_index.fetch_add(1);
    }
}


void generate_vcf_and_run_graphaligner_thread_fn(
        const vector<Region>& regions,
        const path& vcf,
        const path& output_dir,
        const path& ref_fasta,
        const path& fasta_filename,
        atomic<size_t>& job_index
){
    path vcf_to_gfa_script = "/home/ryan/code/hapslap/scripts/vcf_to_gfa.py";

    size_t i = job_index.fetch_add(1);

    while (i < regions.size()){
        const auto& region = regions.at(i);

        path output_subdir = output_dir / region.to_string('_');

        string command = "python3 " + vcf_to_gfa_script.string() +
                         " --vcf " + vcf.string() +
                         " --ref " + ref_fasta.string() +
                         " --output_dir " + output_subdir.string() +
                         " --region " + region.to_string();

        run_command(command, true);

        path gaf_path = output_subdir / "alignments.gaf";

        string name_prefix = vcf.filename().string();
        auto l = name_prefix.find_first_of('.');
        name_prefix = name_prefix.substr(0,l);

        path gfa_path = output_subdir / (name_prefix + ".gfa");
        path fasta_path = output_subdir / fasta_filename;

        command = "GraphAligner"
                  " -x " "vg"
                  " -t " "1"
                  " -a " + gaf_path.string() +
                  " -g " + gfa_path.string() +
                  " -f " + fasta_path.string();

        run_command(command, true);

        i = job_index.fetch_add(1);
    }
}


// TODO: replace this with proper implementation when it exists
void generate_graph_alignments(
        const unordered_map<Region,TransMap>& region_transmaps,
        const path& vcf,
        const path& ref_fasta,
        const vector<Region>& regions,
        size_t n_threads,
        const path& output_dir){

    path fasta_filename = "haplotypes.fasta";

    {
        // Thread-related variables
        atomic<size_t> job_index = 0;
        vector<thread> threads;

        threads.reserve(n_threads);


        // Launch threads
        for (uint64_t n = 0; n < n_threads; n++) {
            try {
                cerr << "launching: " << n << '\n';
                threads.emplace_back(write_region_subsequences_to_file_thread_fn,
                        std::cref(region_transmaps),
                        std::cref(regions),
                        std::cref(output_dir),
                        std::cref(fasta_filename),
                        std::ref(job_index)
                );
            } catch (const exception &e) {
                throw e;
            }
        }

        // Wait for threads to finish
        for (auto &n: threads) {
            n.join();
        }
    }

    {
        // Thread-related variables
        atomic<size_t> job_index = 0;
        vector<thread> threads;

        threads.reserve(n_threads);

        // Launch threads
        for (uint64_t n = 0; n < n_threads; n++) {
            try {
                cerr << "launching: " << n << '\n';
                threads.emplace_back(generate_vcf_and_run_graphaligner_thread_fn,
                        std::cref(regions),
                        std::cref(vcf),
                        std::cref(output_dir),
                        std::cref(ref_fasta),
                        std::cref(fasta_filename),
                        std::ref(job_index)
                );
            } catch (const exception &e) {
                throw e;
            }
        }

        // Wait for threads to finish
        for (auto &n: threads) {
            n.join();
        }
    }

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

    // Generate GFAs and folder structure for every VCF * every region
    // - GFAs should be stored for easy access/debugging later
    // - GAFs should also be stored
    // - FASTA files for extracted truth haplotypes should be optional
    // By default, all of these files will be stored in /dev/shm and then copied into the output dir as a final step.
    // TODO: create option to use /dev/shm/ as staging dir
    // Absolutely must delete the /dev/shm/ copy or warn the user at termination
    //
    path staging_dir = output_dir;
    for (const auto& vcf: vcfs){
        cerr << "Generating graph alignments for VCF: " << vcf << '\n';
        generate_graph_alignments(
            region_transmaps,
            vcf,
            ref_fasta,
            regions,
            n_threads,
            staging_dir);
    }

    cerr << t << "Done" << '\n';
}


int main (int argc, char* argv[]){
    path output_dir;
    path windows_bed;
    path bam_csv;
    path ref_fasta;
    string vcfs_string;
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
            ->required()
            ->expected(1,-1)
            ->delimiter(',');

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

    app.parse(argc, argv);

    auto vcfs = app.get_option("--vcfs")->as<std::vector<path> >();

    evaluate(vcfs, output_dir, windows_bed, bam_csv, ref_fasta, flank_length, n_threads, debug);

    return 0;
}
