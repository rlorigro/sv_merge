#include "path_optimizer_mathopt.hpp"
#include "TransitiveMap.hpp"
#include "interval_tree.hpp"
#include "VariantGraph.hpp"
#include "VcfReader.hpp"
#include "fasta.hpp"
#include "Timer.hpp"
#include "CLI11.hpp"
#include "fetch.hpp"
#include "misc.hpp"
#include "gaf.hpp"
#include "bed.hpp"

using lib_interval_tree::interval_tree_t;

#include "bdsg/hash_graph.hpp"

using bdsg::HandleGraph;

#include <filesystem>

using std::filesystem::path;
using std::filesystem::create_directories;


#include <unordered_map>
#include <exception>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <thread>
#include <limits>
#include <cmath>
#include <tuple>

using std::this_thread::sleep_for;
using std::numeric_limits;
using std::unordered_map;
using std::runtime_error;
using std::exception;
using std::ifstream;
using std::ofstream;
using std::thread;
using std::atomic;
using std::mutex;
using std::log10;
using std::tuple;
using std::cerr;
using std::max;
using std::min;
using std::cref;
using std::ref;


using namespace sv_merge;


#include "bindings/cpp/WFAligner.hpp"
using namespace wfa;


void write_region_subsequences_to_file(
        const unordered_map<Region,TransMap>& region_transmaps,
        const vector<Region>& regions,
        const path& output_dir,
        vector<path>& result_paths
){
    for (auto& region: regions){
        const auto& t = region_transmaps.at(region);

        t.for_each_sample([&](const string& sample_name, int64_t sample_id){
            path output_fasta = output_dir / sample_name / (sample_name + ".fasta");

            if (not std::filesystem::exists(output_fasta.parent_path())){
                std::filesystem::create_directories(output_fasta.parent_path());
                result_paths.emplace_back(output_fasta);
            }

            ofstream file(output_fasta, std::ios_base::app);

            if (not (file.is_open() and file.good())){
                throw runtime_error("ERROR: could not write file: " + output_fasta.string());
            }

            t.for_each_read_of_sample(sample_id, [&](int64_t read_id){
                auto& s = t.get_sequence(read_id);
                if (s.empty()){
                    return;
                }

                file << '>' << t.get_node(read_id).name << '\n';
                file << t.get_sequence(read_id) << '\n';
            });
        });
    }
}


void align_sequences_to_reference_thread_fn(
        path ref_fasta,
        vector<path> fasta_paths,
        path output_dir,
        atomic<size_t>& job_index
        ){

    size_t i = job_index.fetch_add(1);

    while (i < fasta_paths.size()){
        path fasta_path = fasta_paths[i];
        path sam_path = fasta_path.parent_path() / (fasta_path.stem().string() + ".sam");

        string command = "minimap2 -a -x asm20 --eqx -L -t 1 " + ref_fasta.string() + " " + fasta_path.string();
        run_command(command, sam_path);

        // Sort and output as BAM format
        path bam_path = fasta_path.parent_path() / (fasta_path.stem().string() + ".bam");
        command = "samtools sort -@ 1 -o " + bam_path.string() + " " + sam_path.string();
        run_command(command, bam_path);

        // Index the BAM
        command = "samtools index " + bam_path.string();
        run_command(command);

        i = job_index.fetch_add(1);
    }
}


void make_local_test_set(
        path vcf,
        path output_dir,
        path windows_bed,                // Override the interval graph if this is provided
        path tandem_bed,
        path bam_csv,
        path ref_fasta,
        int32_t flank_length,
        int32_t interval_max_length,
        int32_t n_threads,
        bool debug,
        bool force_unique_reads,
        bool bam_not_hardclipped
){
    Timer t;

    if (std::filesystem::exists(output_dir)){
        throw runtime_error("ERROR: output dir exists already: " + output_dir.string());
    }
    else{
        std::filesystem::create_directories(output_dir);
    }

    output_dir = std::filesystem::weakly_canonical(output_dir);
    tandem_bed = std::filesystem::weakly_canonical(tandem_bed);
    bam_csv = std::filesystem::weakly_canonical(bam_csv);
    ref_fasta = std::filesystem::weakly_canonical(ref_fasta);
    vcf = std::filesystem::weakly_canonical(vcf);

    cerr << t << "Initializing" << '\n';

    vector <pair <string,path> > bam_paths;
    unordered_map<string,string> ref_sequences;
    vector<Region> regions;

    cerr << t << "Loading reference sequences" << '\n';

    // Load all chromosome sequences (in case of BND)
    for_sequence_in_fasta_file(ref_fasta, [&](const Sequence& s){
        ref_sequences[s.name] = s.sequence;
    });

    cerr << t << "Reading BED file" << '\n';
    load_windows_from_bed(windows_bed, regions);

    if (not regions.empty()){
        // This is only used while loading VCFs to find where each record belongs
        unordered_map <string, interval_tree_t<int32_t> > contig_interval_trees;

        // Log which windows were used
        path bed_output_path = output_dir / "windows.bed";
        path bed_flanked_output_path = output_dir / "windows_flanked.bed";
        ofstream output_bed_file(bed_output_path);
        ofstream output_bed_flanked_file(bed_flanked_output_path);

        cerr << t << "Flanking windows and writing BED" << '\n';

        // Add flanks, place the regions in the interval tree, and log the windows
        for (auto& r: regions) {
            output_bed_file << r.to_bed() << '\n';

            r.start -= flank_length;
            r.stop += flank_length;

            contig_interval_trees[r.name].insert({r.start, r.stop});
            output_bed_flanked_file << r.to_bed() << '\n';
        }
        output_bed_file.close();
        output_bed_flanked_file.close();

        cerr << t << "Fetching reads for all windows" << '\n';

        // The container to store all fetched reads and their relationships to samples/paths
        unordered_map<Region,TransMap> region_transmaps;

        auto max_length = size_t(float(interval_max_length) * 2.5);

        if (bam_not_hardclipped){
            cerr << "Fetching from NON-hardclipped BAMs" << '\n';
            fetch_reads(
                    t,
                    regions,
                    bam_csv,
                    n_threads,
                    region_transmaps,
                    true,
                    force_unique_reads,
                    true,
                    false,
                    0
            );
        }
        else{
            cerr << "Fetching from HARDCLIPPED BAMs" << '\n';
            fetch_reads_from_clipped_bam(
                    t,
                    regions,
                    bam_csv,
                    n_threads,
                    max_length,
                    flank_length,
                    region_transmaps,
                    true,
                    false,
                    false,
                    force_unique_reads,
                    true,
                    0
            );
        }

        cerr << t << "Writing sequences to disk" << '\n';

        vector<path> fasta_paths;

        // Dump sequences into each region directory
        write_region_subsequences_to_file(region_transmaps, regions, output_dir, fasta_paths);

        // Launch threads to align the sequences to the reference using minimap2 system calls
        cerr << t << "Aligning sequences to reference" << '\n';

        {
            vector<thread> threads;
            atomic<size_t> job_index = 0;

            threads.reserve(n_threads);

            for (uint64_t n=0; n<n_threads; n++){
                threads.emplace_back(
                        align_sequences_to_reference_thread_fn,
                        ref_fasta,
                        std::cref(fasta_paths),
                        output_dir,
                        std::ref(job_index)
                );
            }

            for (auto& n: threads){
                n.join();
            }
        }

        cerr << t << "Writing VCF records to disk" << '\n';

        // Use bcftools command to subset to the regions of the BED provided
        path vcf_subset_path = output_dir / "subset.vcf";

        string command = "bcftools view -R " + bed_output_path.string() + " " + vcf.string();

        run_command(command, vcf_subset_path);
    }

    cerr << t << "Peak memory usage: " << get_peak_memory_usage() << '\n';
    cerr << t << "Done" << '\n';
}


int main (int argc, char* argv[]){
    path output_dir;
    path windows_bed;
    path tandem_bed;
    string bam_csv;
    path ref_fasta;
    path vcf;
    int32_t flank_length = 150;
    int32_t interval_max_length = 15000;
    int32_t min_sv_length = 20;
    int32_t n_threads = 1;
    bool debug = false;
    bool force_unique_reads = false;
    bool bam_not_hardclipped = false;

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
            "--vcf",
            vcf,
            "Path to VCF file containing variants to annotate (must be biallelic/normed)");

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

    app.add_flag("--force_unique_reads", force_unique_reads, "Invoke this to add append each read name with the sample name so that inter-sample read collisions cannot occur");

    app.add_flag("--bam_not_hardclipped", bam_not_hardclipped, "Invoke this if you expect your BAMs NOT to contain ANY hardclips. Saves time on iterating.");

    try{
        app.parse(argc, argv);
    }
    catch (const CLI::ParseError &e) {
        return app.exit(e);
    }

    make_local_test_set(
            vcf,
            output_dir,
            windows_bed,
            tandem_bed,
            bam_csv,
            ref_fasta,
            flank_length,
            interval_max_length,
            n_threads,
            debug,
            force_unique_reads,
            bam_not_hardclipped
    );

    return 0;
}
