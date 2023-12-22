#include "read_cluster_optimizer.hpp"
#include "TransitiveMap.hpp"
#include "windows.hpp"
#include "Timer.hpp"
#include "CLI11.hpp"
#include "fetch.hpp"
#include "misc.hpp"
#include "bed.hpp"

#include "bindings/cpp/WFAligner.hpp"
using namespace wfa;

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


void cross_align_sample_reads(TransMap& transmap, int64_t score_threshold, const string& sample_name, int64_t flank_length){
    WFAlignerGapAffine aligner(4,6,2,WFAligner::Alignment,WFAligner::MemoryHigh);

    auto sample_id = transmap.get_id(sample_name);
    transmap.for_each_read_of_sample(sample_id, [&](int64_t a){
        transmap.for_each_read_of_sample(sample_id, [&](int64_t b){
            if (transmap.has_edge(a,b)){
                return;
            }

            auto& seq_a = transmap.get_sequence(a);
            auto& seq_b = transmap.get_sequence(b);

            if (llabs(int64_t(seq_a.size()) - int64_t(seq_b.size())) > score_threshold){
//                cerr << "skipping edge: " << a << ',' << b << " l: " << int64_t(seq_a.size()) << ',' << int64_t(seq_b.size()) << '\n';
                return;
            }

            if (a != b) {
                aligner.alignEnd2End(seq_a, seq_b);
                auto score = aligner.getAlignmentScore();

//                // WFA creates negative scores for "distance", we make it positive again and rescale it as a percent
                size_t scaled_score = 100 - (100*size_t(-score)) / min(seq_a.size(), seq_b.size());

//                int64_t scaled_score = score > -100 ? 1 : -1;

                transmap.add_edge(a,b, float(scaled_score));

                cerr << transmap.get_node(a).name << ',' << transmap.get_node(b).name << ' ' << int64_t(seq_a.size()) << ',' << int64_t(seq_b.size()) << ' ' << score << ' ' << scaled_score << '\n';
            }
        });
    });
}


void hapestry(
        path output_dir,
        path windows_bed,               // Override the interval graph if this is provided
        path tandem_bed,
        path bam_csv,
        path vcf,
        path reference,
        int32_t interval_max_length,
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

    // TODO: use percent of min(a,b) where a,b are lengths of seqs?
    int64_t score_threshold = 200;

    cerr << t << "Constructing windows" << '\n';

    if (windows_bed.empty()){
        construct_windows_from_vcf_and_bed(tandem_bed, vcf, flank_length, interval_max_length, regions);
    }
    else{
        load_windows_from_bed(windows_bed, regions);
    }

    for (auto& r: regions) {
        r.start -= flank_length;
        r.stop += flank_length;
    }

    // The container to store all fetched reads and their relationships to samples/paths
    unordered_map<Region,TransMap> region_transmaps;

    cerr << t << "Fetching reads for all windows" << '\n';

    fetch_reads(t, regions, bam_csv, n_threads, false, region_transmaps);

//    throw runtime_error("DEBUG");

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

        // Iterate samples within this region and cluster reads
        transmap.for_each_sample([&](const string& sample_name, int64_t sample_id) {
            if (debug) {
                cerr << '\n';
                cerr << sample_name << '\n';

                path p = output_subdir / (sample_name + "_extracted_reads.fasta");
                ofstream file(p);

                transmap.for_each_read_of_sample(sample_name, [&](const string &name, int64_t id) {
                    file << '>' << name << '\n';
                    file << transmap.get_sequence(id) << '\n';
                });
            }

            // Reads must be aligned all-vs-all to get scoring for cluster scheme
            cross_align_sample_reads(transmap, score_threshold, sample_name, flank_length);

            CpModelBuilder model;
            ReadClusterVariables vars;
            vector <vector<int64_t> > clusters;

            // Cluster reads within samples
            int64_t n_weight = 1;
            int64_t d_weight = 2;
            cluster_reads(transmap, sample_id, d_weight, n_weight, clusters);

            if (debug){
                cerr << sample_name << '\n';

                path p = output_subdir / (sample_name + "_clusters.txt");
                ofstream file(p);

                int c = 0;
                for (const auto& cluster: clusters){
                    cerr << c << '\n';
                    file << c << '\n';
                    for (const auto& id: cluster){
                        cerr << '\t' << id << ' ' << transmap.get_node(id).name << '\n';
                        file << transmap.get_sequence(id) << '\n';
                    }

                    c++;
                }

//                for (auto a: representatives){
//                    cerr << transmap.get_node(a).name << '\n';
//                    file << transmap.get_node(a).name << '\n';
//                    cerr << transmap.get_sequence(a).substr(0, 100) << '\n';
//                    file << transmap.get_sequence(a) << '\n';
//                    transmap.for_each_neighbor_of_type(a, 'R', [&](int64_t b){
//                        cerr << b << ' ' << transmap.get_node(b).name << '\n';
//                        cerr << transmap.get_sequence(b).substr(0, 100) << '\n';
//                        file << transmap.get_sequence(b) << '\n';
//                    });
//                }
            }
        });
    }

    cerr << t << "Done" << '\n';

}


int main (int argc, char* argv[]){
    path output_dir;
    path windows_bed;
    path tandem_bed;
    path bam_csv;
    path vcf;
    path ref;
    int32_t interval_max_length;
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
            "--vcf",
            vcf,
            "Path to VCF file containing variants to be merged");

    app.add_option(
            "--tandems",
            tandem_bed,
            "Path to BED file containing tandem repeat locations (for automated window generation)");

    app.add_option(
            "--windows",
            windows_bed,
            "Path to BED file containing windows to merge (inferred automatically if not provided)");

    app.add_option(
            "--bam_csv",
            bam_csv,
            "Simple headerless CSV file with the format [sample_name],[bam_path]")
            ->required();

    app.add_option(
            "--interval_max_length",
            interval_max_length,
            "Maximum reference window size")
            ->required();

    app.add_option(
            "--flank_length",
            flank_length,
            "How much flanking sequence to use when fetching and aligning reads")
            ->required();

    app.add_flag("--debug", debug, "Invoke this to add more logging and output");

    CLI11_PARSE(app, argc, argv);

    hapestry(
        output_dir,
        windows_bed,
        tandem_bed,
        bam_csv,
        vcf,
        ref,
        interval_max_length,
        flank_length,
        n_threads,
        debug
    );

    return 0;
}
