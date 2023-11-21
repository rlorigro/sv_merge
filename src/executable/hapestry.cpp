#include "read_optimizer.hpp"
#include "TransitiveMap.hpp"
#include "IntervalGraph.hpp"
#include "Authenticator.hpp"
#include "Timer.hpp"
#include "CLI11.hpp"
#include "misc.hpp"
#include "bed.hpp"
#include "bam.hpp"

#include "bindings/cpp/WFAligner.hpp"
using namespace wfa;

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


using namespace sv_merge;


void for_each_sample_bam_path(path bam_csv, const function<void(const string& sample_name, const path& bam_path)>& f){
    if (not (bam_csv.extension() == ".csv")){
        throw runtime_error("ERROR: file does not have compatible csv extension: " + bam_csv.string());
    }

    ifstream file(bam_csv);

    if (not (file.is_open() and file.good())){
        throw runtime_error("ERROR: could not read file: " + bam_csv.string());
    }

    char c;
    string sample_name;
    string bam_path;

    int64_t n_delimiters = 0;
    char delimiter = ',';

    while (file.get(c)){
//        cerr << c << ' ' << sample_name << ' ' << bam_path << '\n';
        if (c == delimiter){
            n_delimiters++;
            continue;
        }
        if (c == '\n'){
            f(sample_name, bam_path);

            sample_name.clear();
            bam_path.clear();
            n_delimiters = 0;
            continue;
        }

        if (n_delimiters == 0){
            sample_name += c;
        }
        else if (n_delimiters == 1){
            bam_path += c;
        }
        else {
            throw runtime_error("ERROR: too many delimiters in bam csv");
        }
    }
}


void construct_windows_from_vcf_and_bed(path tandem_bed, path vcf, vector<Region>& regions){
    vector<labeled_interval_t> intervals;

    // Iterate the VCF file and construct a vector of labeled intervals for the IntervalGraph
    // Need to append `interval_padding` onto intervals and then subtract it afterwards
    // TODO: fill in when vcf reader exists
    IntervalGraph<string> g(intervals);

    g.for_each_connected_component_interval([&](interval_t& interval, unordered_set<string>& values){
        cerr << interval.first << "," << interval.second;
        for (const auto& v: values){
            cerr << ' ' << v;
        }
        cerr << '\n';
    });
}


// This is used when only one half of the ref/query dual iterator `for_alignment_in_bam_region` is desired
void null_fn(const CigarInterval &intersection, const interval_t &interval){}


void extract_subsequences_from_region(
        GoogleAuthenticator& authenticator,
        mutex& authenticator_mutex,
        const string& sample_name,
        const Region& region,
        path bam_path,
        TransMap& transmap,
        mutex& transmap_mutex
){

    CigarInterval placeholder;
    placeholder.query_start = numeric_limits<int64_t>::max();
    placeholder.query_stop = numeric_limits<int64_t>::min();
    placeholder.ref_start = numeric_limits<int64_t>::max();
    placeholder.ref_stop = numeric_limits<int64_t>::min();

    // Keep track of the min and max observed query coordinates that intersect the region of interest
    unordered_map<string,CigarInterval> query_coords;
    unordered_map<string,string> query_sequences;

    // Keep the F (query) oriented sequence of all the reads if they have an alignment that intersects the region
    unordered_map<string,string> alignments;

    // The region of interest is defined in reference coordinate space
    vector<interval_t> ref_intervals = {{region.start, region.stop}};

    // Unused
    vector<interval_t> query_intervals;

    // Make sure the system has the necessary authentication env variable to fetch a remote GS URI
    authenticator_mutex.lock();
    authenticator.update();
    authenticator_mutex.unlock();

    // Iterate each alignment in the ref region
    for_alignment_in_bam_region(bam_path, region.to_string(), [&](Alignment& alignment) {
        if (alignment.is_unmapped() or not alignment.is_primary()){
            return;
        }

        string name;
        alignment.get_query_name(name);

        // Check if this read/query has an existing coord, from a previously iterated supplementary alignment
        auto result = query_coords.find(name);

        // If no previous alignment, initialize with the max placeholder
        if (result == query_coords.end()){
            // Insert a new placeholder cigar interval and keep a reference to the value inserted
            result = query_coords.emplace(name, placeholder).first;
            result->second.is_reverse = alignment.is_reverse();

            // Insert a new empty sequence and keep a reference to the value inserted
            auto result3 = query_sequences.emplace(name,"");
            auto& x = result3.first->second;

            // Fill the value with the sequence
            // TODO: find a way to not extract the whole sequence? Tricky when using the Alignment abstract class, and
            // fetching from remote BAM
            alignment.get_query_sequence(x);
        }

        auto& coord = result->second;

        // Find the widest possible pair of query coordinates which exactly spans the ref region (accounting for DUPs)
        for_cigar_interval_in_alignment(alignment, ref_intervals, query_intervals,
            [&](const CigarInterval& intersection, const interval_t& interval) {
//                cerr << cigar_code_to_char[intersection.code] << ' ' << alignment.is_reverse() << " r: " << intersection.ref_start << ',' << intersection.ref_stop << ' ' << "q: " << intersection.query_start << ',' << intersection.query_stop << '\n';

                // If the alignment touches the START of the ref region, record the query position
                if (intersection.ref_start == region.start){
                    auto [start,stop] = intersection.get_forward_query_interval();

                    if (alignment.is_reverse()){
                        if (stop > coord.query_stop){
                            coord.query_stop = stop;
                        }
                    }
                    else{
                        if (start < coord.query_start){
                            coord.query_start = start;
                        }
                    }
                }

                // If the alignment touches the END of the region, record the query position
                if (intersection.ref_stop == region.stop){
                    auto [start,stop] = intersection.get_forward_query_interval();

                    if (intersection.is_reverse){
                        if (start < coord.query_start){
                            coord.query_start = start;
                        }
                    }
                    else{
                        if (stop > coord.query_stop){
                            coord.query_stop = stop;
                        }
                    }
                }
            },
            null_fn
        );
    });

    for (const auto& [name, coords]: query_coords){
        if (coords.query_start != placeholder.query_start and coords.query_stop != placeholder.query_stop){
            auto i = coords.query_start;
            auto l = coords.query_stop - coords.query_start;

//            cerr << name << ' ' << coords.is_reverse << ' ' << l << ' ' << coords.query_start << ',' << coords.query_stop << '\n';

            if (coords.is_reverse) {
                auto s = query_sequences[name].substr(i, l);
                reverse_complement(s);

                transmap_mutex.lock();
                transmap.add_read(name, s);
                transmap_mutex.unlock();
            }
            else{
                transmap_mutex.lock();
                transmap.add_read(name, query_sequences[name].substr(i, l));
                transmap_mutex.unlock();
            }

            transmap_mutex.lock();
            transmap.add_edge(sample_name, name);
            transmap_mutex.unlock();
        }
    }
}


void extract_subsequences_from_region_thread_fn(
        GoogleAuthenticator& authenticator,
        mutex& authenticator_mutex,
        const vector <pair <string,path> >& bam_paths,
        const Region& region,
        TransMap& transmap,
        mutex& transmap_mutex,
        atomic<size_t>& job_index
){
    size_t i = job_index.fetch_add(1);

    while (i < bam_paths.size()){
        const auto& [sample_name, bam_path] = bam_paths[i];

        Timer t;

        extract_subsequences_from_region(
            authenticator,
            authenticator_mutex,
            sample_name,
            region,
            bam_path,
            transmap,
            transmap_mutex
            );

        cerr << t << "Elapsed for: " << sample_name << ' ' << region.to_string() << '\n';

        i = job_index.fetch_add(1);
    }

}


void cross_align_sample_reads(TransMap& transmap, int64_t score_threshold, const string& sample_name){
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

                // WFA creates negative scores for "distance", we make it positive again
                transmap.add_edge(a,b, float(-score));

                cerr << transmap.get_node(a).name << ',' << transmap.get_node(b).name << ' ' << score << '\n';
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
        int64_t interval_padding,
        int64_t interval_max_length,
        int64_t flank_length,
        int64_t n_threads
        ){
    // Flag to control how much logging/dumping
    bool debug = true;

    if (ghc::filesystem::exists(output_dir)){
        throw runtime_error("ERROR: output dir exists already: " + output_dir.string());
    }
    else{
        ghc::filesystem::create_directories(output_dir);
    }

    Timer t;

    cerr << t << "Initializing" << '\n';

    vector <pair <string,path> > bam_paths;
    GoogleAuthenticator authenticator;
    vector<Region> regions;

    // TODO: reserve the hash tables used in transmap so that they approximately have the space needed for n_samples*n_fold_coverage ?
    TransMap sample_only_transmap;

    // TODO: use percent of min(a,b) where a,b are lengths of seqs?
    int64_t score_threshold = 200;

    cerr << t << "Constructing windows" << '\n';

    if (windows_bed.empty()){
        construct_windows_from_vcf_and_bed(tandem_bed, vcf, regions);
    }
    else{
        for_region_in_bed_file(tandem_bed, [&](const Region &r) {
            regions.push_back(r);
        });
    }

    cerr << t << "Loading CSV" << '\n';

    // Load BAM paths as a map with sample->bam
    for_each_sample_bam_path(bam_csv, [&](const string& sample_name, const path& bam_path){
        sample_only_transmap.add_sample(sample_name);
        bam_paths.emplace_back(sample_name,bam_path);
    });

    cerr << t << "Processing windows" << '\n';

    unordered_map<Region,TransMap> region_transmaps;

    // For each region
    for (auto region: regions){
        path output_subdir = output_dir / region.to_string('_');
        ghc::filesystem::create_directories(output_subdir);

        // Duplicate the base transmap which already has the samples loaded
        auto& transmap = region_transmaps.emplace(region, sample_only_transmap).first->second;

        cerr << t << region.name << ' ' << region.start << ' ' << region.stop << '\n';

        region.start -= flank_length;
        region.stop += flank_length;

        // Thread-related variables
        atomic<size_t> job_index = 0;
        vector<thread> threads;
        mutex transmap_mutex;
        mutex authenticator_mutex;

        threads.reserve(n_threads);

        // Launch threads
        for (uint64_t n=0; n<n_threads; n++){
            try {
                cerr << "launching: " << n << '\n';
                threads.emplace_back(
                    extract_subsequences_from_region_thread_fn,
                    std::ref(authenticator),
                    std::ref(authenticator_mutex),
                    std::cref(bam_paths),
                    std::cref(region),
                    std::ref(transmap),
                    std::ref(transmap_mutex),
                    std::ref(job_index)
                );
            } catch (const exception& e) {
                cerr << e.what() << "\n";
                exit(1);
            }
        }

        // Wait for threads to finish
        for (auto& n: threads){
            n.join();
        }


        for (auto& [sample_name, _]: bam_paths){
            cerr << sample_name << ' ' << region.to_string() << '\n';
            auto sample_id = transmap.get_id(sample_name);

            if (debug){
                path p = output_subdir / (sample_name + "_extracted_reads.fasta");
                ofstream file(p);

                transmap.for_each_read_of_sample(sample_name, [&](const string& name, int64_t id){
                    file << '>' << name << '\n';
                    file << transmap.get_sequence(id) << '\n';
                });
            }

            cross_align_sample_reads(transmap, score_threshold, sample_name);

            CpModelBuilder model;
            ReadVariables vars;
            vector<int64_t> representatives;

            optimize_reads_with_d_and_n(transmap, sample_id, representatives);

            if (debug){
                path p = output_subdir / (sample_name + "_clusters.txt");
                ofstream file(p);

                for (auto a: representatives){
                    cerr << transmap.get_node(a).name << '\n';
                    file << transmap.get_node(a).name << '\n';
                    cerr << transmap.get_sequence(a).substr(0, 100) << '\n';
                    file << transmap.get_sequence(a) << '\n';
                    transmap.for_each_neighbor_of_type(a, 'R', [&](int64_t b){
                        cerr << b << ' ' << transmap.get_node(b).name << '\n';
                        cerr << transmap.get_sequence(b).substr(0, 100) << '\n';
                        file << transmap.get_sequence(b) << '\n';
                    });
                }
            }


        }
    }

    cerr << t << "Done" << '\n';

}


int main (int argc, char* argv[]){
    path output_dir;
    path windows_bed;
    path tandem_bed;
    path bam_csv;
    path vcf = ""; // TODO: add arg for this when vcf reader exists
    path ref;
    int64_t interval_padding;
    int64_t interval_max_length;
    int64_t flank_length;
    int64_t n_threads = 1;

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
            "--tandems",
            windows_bed,
            "Path to BED file containing tandem repeat locations (for automated window generation)")
            ->required();

    app.add_option(
            "--windows",
            tandem_bed,
            "Path to BED file containing windows to merge (inferred automatically if not provided)")
            ->required();

    app.add_option(
            "--bam_csv",
            bam_csv,
            "Simple headerless CSV file with the format [sample_name],[bam_path]")
            ->required();

    app.add_option(
            "--interval_padding",
            interval_padding,
            "How much space to require between adjacent ref windows (if window generation is automated)")
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

    CLI11_PARSE(app, argc, argv);

    hapestry(
        output_dir,
        windows_bed,
        tandem_bed,
        bam_csv,
        vcf,
        ref,
        interval_padding,
        interval_max_length,
        flank_length,
        n_threads
    );

    return 0;
}
