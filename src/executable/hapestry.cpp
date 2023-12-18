#include "read_cluster_optimizer.hpp"
#include "TransitiveMap.hpp"
#include "IntervalGraph.hpp"
#include "Authenticator.hpp"
#include "VcfReader.hpp"
#include "windows.hpp"
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
using std::max;
using std::min;
using std::cref;
using std::ref;


using namespace sv_merge;

using sample_region_read_map_t = unordered_map <string, unordered_map <Region, vector<Sequence> > >;


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


// This is used when only one half of the ref/query dual iterator `for_alignment_in_bam_region` is desired
void null_fn(const CigarInterval &intersection, const interval_t &interval){}


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


void extract_subregions_from_sample(
        GoogleAuthenticator& authenticator,
        mutex& authenticator_mutex,
        sample_region_read_map_t& sample_to_region_reads,
        const string& sample_name,
        const vector<Region>& subregions,
        path bam_path
){
    if (subregions.empty()){
        throw runtime_error("ERROR: subregions empty");
    }

    // Generate a super-region to encompass all subregions, and assume that subregions are sorted, contiguous.
    // If they are not contiguous and sorted, the iterator function will detect that and error out.
    Region super_region;
    super_region.name = subregions[0].name;
    super_region.start = subregions[0].start;
    super_region.stop = subregions.back().stop;

    CigarInterval placeholder;
    placeholder.query_start = numeric_limits<int64_t>::max();
    placeholder.query_stop = numeric_limits<int64_t>::min();
    placeholder.ref_start = numeric_limits<int64_t>::max();
    placeholder.ref_stop = numeric_limits<int64_t>::min();

    // Keep track of the min and max observed query coordinates that intersect the region of interest
    unordered_map <Region, unordered_map<string,CigarInterval> > query_coords_per_region;
    unordered_map<string,string> query_sequences;

    // Keep the F (query) oriented sequence of all the reads if they have an alignment that intersects the region
    unordered_map <Region, unordered_map<string,string> > alignments_per_region;

    // Unused
    vector<interval_t> query_intervals;

    // Make sure the system has the necessary authentication env variable to fetch a remote GS URI
    authenticator_mutex.lock();
    authenticator.update();
    authenticator_mutex.unlock();

    // Iterate each alignment in the ref region
    for_alignment_in_bam_subregions(
            bam_path,
            super_region.to_string(),
            subregions,
            [&](Alignment& alignment, span<const Region>& overlapping_regions){

        if (alignment.is_unmapped() or not alignment.is_primary()){
            return;
        }

        string name;
        alignment.get_query_name(name);

//        cerr << name << ' ' << alignment.get_ref_start() << ' ' << alignment.get_ref_stop() << '\n';
//        for (const auto& item: overlapping_regions){
//            cerr << item.start << ',' << item.stop << '\n';
//        }

        // The region of interest is defined in reference coordinate space
        vector<interval_t> ref_intervals;
        for (auto& r: overlapping_regions){
            ref_intervals.emplace_back(r.start, r.stop);

            // Find/or create coord for this region
            auto& region_coord = query_coords_per_region[r];

            // Check if this read/query has an existing coord, from a previously iterated supplementary alignment
            auto result = region_coord.find(name);

            // If no previous alignment, initialize with the max placeholder
            if (result == region_coord.end()){
                // Insert a new placeholder cigar interval and keep a reference to the value inserted
                result = region_coord.emplace(name, placeholder).first;
                result->second.is_reverse = alignment.is_reverse();

                // Try inserting a new empty sequence and keep a reference to the value inserted
                auto [iter,success] = query_sequences.try_emplace(name,"");

                // If the value existed already, don't do anything, the sequence has already been extracted
                if (not success){
                    continue;
                }

                auto& x = iter->second;

                // Fill the value with the sequence
                alignment.get_query_sequence(x);
            }
        }

        // Find the widest possible pair of query coordinates which exactly spans the ref region (accounting for DUPs)
        for_cigar_interval_in_alignment(alignment, ref_intervals, query_intervals,
            [&](const CigarInterval& intersection, const interval_t& interval) {
//                cerr << cigar_code_to_char[intersection.code] << ' ' << alignment.is_reverse() << " r: " << intersection.ref_start << ',' << intersection.ref_stop << ' ' << "q: " << intersection.query_start << ',' << intersection.query_stop << '\n';

                for (auto& region: overlapping_regions){
                    auto& coord = query_coords_per_region.at(region).at(name);

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
                }
            },
            null_fn
        );
    });

    // Finally trim the sequences and insert the subsequences into a map which has keys pre-filled
    for (const auto& [region, query_coords]: query_coords_per_region){
        for (const auto& [name, coords]: query_coords){
            if (coords.query_start != placeholder.query_start and coords.query_stop != placeholder.query_stop){
                auto i = coords.query_start;
                auto l = coords.query_stop - coords.query_start;

//                cerr << name << ' ' << coords.is_reverse << ' ' << l << ' ' << coords.query_start << ',' << coords.query_stop << '\n';

                if (coords.is_reverse) {
                    auto s = query_sequences[name].substr(i, l);
                    reverse_complement(s);

                    sample_to_region_reads.at(sample_name).at(region).emplace_back(name,s);
                }
                else{
                    sample_to_region_reads.at(sample_name).at(region).emplace_back(name,query_sequences[name].substr(i, l));
                }
            }
        }
    }
}


void extract_subsequences_from_sample_thread_fn(
        GoogleAuthenticator& authenticator,
        mutex& authenticator_mutex,
        sample_region_read_map_t & sample_to_region_reads,
        const vector <pair <string,path> >& sample_bams,
        const vector<Region>& regions,
        atomic<size_t>& job_index
){

    size_t i = job_index.fetch_add(1);

    while (i < sample_bams.size()){
        const auto& [sample_name, bam_path] = sample_bams[i];

        Timer t;

        extract_subregions_from_sample(
            authenticator,
            authenticator_mutex,
            sample_to_region_reads,
            sample_name,
            regions,
            bam_path
        );

        cerr << t << "Elapsed for: " << sample_name << '\n';

        i = job_index.fetch_add(1);
    }
}


void get_reads_for_each_bam_subregion(
        Timer& t,
        vector<Region>& regions,
        GoogleAuthenticator& authenticator,
        sample_region_read_map_t& sample_to_region_reads,
        path bam_csv,
        int64_t n_threads
){
    // Intermediate objects
    vector <pair <string, path> > sample_bams;
    TransMap sample_only_transmap;

    cerr << t << "Loading CSV" << '\n';

    // Load BAM paths as a map with sample->bam
    for_each_sample_bam_path(bam_csv, [&](const string& sample_name, const path& bam_path){
        sample_only_transmap.add_sample(sample_name);
        sample_bams.emplace_back(sample_name, bam_path);

        // Initialize every combo of sample,region with an empty vector
        for (const auto& region: regions){
            sample_to_region_reads[sample_name][region] = {};
        }
    });

    cerr << t << "Processing windows" << '\n';

    // Thread-related variables
    atomic<size_t> job_index = 0;
    vector<thread> threads;
    mutex authenticator_mutex;

    threads.reserve(n_threads);

    // Launch threads
    for (uint64_t n=0; n<n_threads; n++){
        try {
            cerr << "launching: " << n << '\n';
            threads.emplace_back(
                std::ref(extract_subsequences_from_sample_thread_fn),
                std::ref(authenticator),
                std::ref(authenticator_mutex),
                std::ref(sample_to_region_reads),
                std::cref(sample_bams),
                std::cref(regions),
                std::ref(job_index)
            );
        } catch (const exception& e) {
            throw e;
        }
    }

    // Wait for threads to finish
    for (auto& n: threads){
        n.join();
    }

}


void fetch_reads(
        Timer& t,
        vector<Region>& regions,
        path bam_csv,
        int64_t n_threads,
        unordered_map<Region,TransMap>& region_transmaps
){
    GoogleAuthenticator authenticator;
    TransMap template_transmap;

    // Intermediate object to store results of multithreaded sample read fetching
    sample_region_read_map_t sample_to_region_reads;

    // Get a nested map of sample -> region -> vector<Sequence>
    get_reads_for_each_bam_subregion(t, regions, authenticator, sample_to_region_reads, bam_csv, n_threads);

    // Construct template transmap with only samples
    for (const auto& [sample_name,item]: sample_to_region_reads){
        template_transmap.add_sample(sample_name);
    }

    // Compute coverages for regions
    unordered_map<Region, size_t> region_coverage;
    for (const auto& [sample_name,item]: sample_to_region_reads) {
        for (const auto& [region,sequences]: item) {
            region_coverage[region] += sequences.size();
        }
    }

    // Copy the template transmap into every region and reserve approximate space for nodes/edges/sequences
    region_transmaps.reserve(regions.size());
    for (const auto& r: regions){
        auto item = region_transmaps.emplace(r, template_transmap).first->second;
        auto n_reads = region_coverage[r];
        auto n_samples = sample_to_region_reads.size();

        // Number of reads is known exactly
        item.reserve_sequences(n_reads + 2);

        // An approximate upper limit that assumes every sample has 1 unique haplotype
        item.reserve_nodes(n_reads + n_samples*2);

        // An approximate upper limit that assumes every sample has 1 unique haplotype, fully connected in read->hap
        item.reserve_edges(n_reads * n_samples);
    }

    // Move the downloaded data into a transmap, construct edges for sample->read
    for (auto& [sample_name,item]: sample_to_region_reads){
        for (auto& [region,sequences]: item){
            auto& transmap = region_transmaps.at(region);

            auto sample_id = transmap.get_id(sample_name);

            for (auto& s: sequences) {
                transmap.add_read_with_move(s.name, s.sequence);
                transmap.add_edge(sample_id, transmap.get_id(s.name), 0);
            }
        }
    }
}


void load_windows_from_bed(path windows_bed, vector<Region>& regions){
    for_region_in_bed_file(windows_bed, [&](const Region &r) {
        regions.push_back(r);
    });

    // How to sort regions
    auto left_comparator = [](const Region& a, const Region& b){
        if (a.name != b.name) {
            return a.name < b.name;
        }
        else {
            return a.start < b.start;
        }
    };

    sort(regions.begin(), regions.end(), left_comparator);
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

    fetch_reads(t, regions, bam_csv, n_threads, region_transmaps);

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

                throw runtime_error("DEBUG EARLY EXIT");

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
