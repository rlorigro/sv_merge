#include "TransitiveMap.hpp"
#include "Filesystem.hpp"
#include "Timer.hpp"
#include "CLI11.hpp"
#include "fetch.hpp"
#include "misc.hpp"
#include "gaf.hpp"
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


class AlignmentSummary{
public:
    int32_t start;
    int32_t stop;
    int32_t n_match;
    int32_t n_mismatch;
    int32_t n_insert;
    int32_t n_delete;

    void update(const CigarInterval& cigar_interval);
};


void AlignmentSummary::update(const sv_merge::CigarInterval& c) {
    switch (c.code){
        case 7:
            n_match += c.length;    // =
        case 8:
            n_mismatch += c.length; // X
        case 1:
            n_insert += c.length;   // I
        case 2:
            n_delete += c.length;   // D
        default:
            return;
    }
}


class GafSummary{
public:
    unordered_map <string, vector<AlignmentSummary> > ref_summaries;
    unordered_map <string, vector<AlignmentSummary> > query_summaries;

    /**
     * Update alignment stats, for a given ref node
     * @param ref_name : ref node to be assigned this cigar operation
     * @param c : cigar interval obtained from iterating an alignment
     * @param insert : a new alignment should start with this operation (separate alignments will be overlap-resolved)
     */
    void update_ref(const string& ref_name, const CigarInterval& c, bool insert);

    /**
     * Update alignment stats, for a given query node
     * @param query_name : query node to be assigned this cigar operation
     * @param c : cigar interval obtained from iterating an alignment
     * @param insert : a new alignment should start with this operation (separate alignments will be overlap-resolved)
     */
    void update_query(const string& query_name, const CigarInterval& c, bool insert);

    void resolve_all_overlaps();
    void resolve_overlaps(vector<AlignmentSummary>& alignments);
};


void GafSummary::update_ref(const string& ref_name, const CigarInterval& c, bool insert){
    auto& result = ref_summaries.at(ref_name);
    if (insert){
        result.emplace_back();
    }
    // TODO: add ref/query specific rules to update fn for start/stop coord, using reversal flag
    result.back().update(c);
}


void GafSummary::update_query(const string& query_name, const CigarInterval& c, bool insert){
    auto& result = query_summaries.at(query_name);
    if (insert){
        result.emplace_back();
    }
    // TODO: add ref/query specific rules to update fn for start/stop coord, using reversal flag
    result.back().update(c);
}


/**
 * Sweep algorithm to split intervals into intersections
 * @param alignments
 * @param is_ref
 */
void GafSummary::resolve_overlaps(vector<AlignmentSummary>& alignments, bool is_ref) {
    vector <pair <interval_t,unordered_set<size_t> > > labeled_intervals;
    vector <pair <interval_t,unordered_set<size_t> > > intersected_intervals;

    for (size_t i=0; i<alignments.size(); i++){
        const auto& a = alignments[i];
        labeled_intervals.push_back({{a.start,a.stop},{i}});
    }

    // How to sort labeled intervals with structure ((a,b), label) by start (a)
    auto left_comparator = [](const pair <interval_t,unordered_set<size_t> >& a, const pair <interval_t,unordered_set<size_t> >& b){
        return a.first.first < b.first.first;
    };

    // How to sort intervals with structure (a,b) by end (b)
    auto right_comparator = [](const interval_t& a, const interval_t& b){
        return a.second < b.second;
    };

    sort(labeled_intervals.begin(), labeled_intervals.end(), left_comparator);

    // Only need to store/compare the interval (and not the data/label) for this DS
    set<interval_t, decltype(right_comparator)> active_intervals;

    for (const auto& [interval,value] : labeled_intervals){
        if (interval.first > interval.second){
            throw runtime_error("ERROR: interval start is greater than interval stop: " + to_string(interval.first) + ',' + to_string(interval.second));
        }

        vector<interval_t> to_be_removed;

        // TODO: Every time an element is added or removed, update the vector of intersected_intervals:
        //  1. extend the previous interval
        //  2. remove/add an element from the set (check if the start/stop is not identical to prev)

        // Iterate active intervals, which are maintained sorted by interval stop
        for (const auto& other_interval: active_intervals){
            // Flag any expired intervals for deletion
            if (other_interval.second < interval.first){
                to_be_removed.emplace_back(other_interval);
            }
        }

        // Remove the intervals that have been passed already in the sweep
        for (const auto& item: to_be_removed){
            active_intervals.erase(item);
        }

        active_intervals.emplace(interval);
    }
}


void GafSummary::resolve_all_overlaps() {
    for (auto& [name, alignments]: ref_summaries){
        resolve_overlaps(alignments);
    }
    for (auto& [name, alignments]: query_summaries){
        resolve_overlaps(alignments);
    }
}


void compute_summaries_from_gaf(
        const path& gaf_path,
        const unordered_map<string,string>& query_sequences,
        int32_t flank_length,
        GafSummary& gaf_summary
        ){

    // GAF alignments (in practice) always progress in the F orientation of the query. Reverse alignments are possible,
    // but redundant because the "reference" is actually a path, which is divided into individually reversible nodes.
    // Both minigraph and GraphAligner code do not allow for R alignments, and instead they use the path orientation.
    bool unclip_coords = false;

    vector<interval_t> query_intervals = {};

    for_alignment_in_gaf(gaf_path, [&](GafAlignment& alignment){
        string name;
        alignment.get_query_name(name);

        vector<interval_t> ref_intervals;
        unordered_map<interval_t,string> interval_to_node_name;

        // First establish the set of intervals that defines the steps in the path (node lengths)
        int32_t x = 0;
        alignment.for_step_in_path(name, [&](const string& step_name, bool is_reverse){
            auto l = int32_t(query_sequences.at(step_name).size());
            ref_intervals.emplace_back(x, x+l);
            interval_to_node_name.emplace(ref_intervals.back(), step_name);
            x += l;
        });

        cerr << name << '\n';

        int64_t length;
        string prev_node_name;
        CigarInterval intersection;

        for_cigar_interval_in_alignment(unclip_coords, alignment, ref_intervals, query_intervals,
        [&](const CigarInterval& i, const interval_t& interval){
            intersection = i;

            auto node_name = interval_to_node_name.at(interval);

            if (is_ref_move[intersection.code]){
                length = intersection.ref_stop - intersection.ref_start;
            }
            else{
                length = intersection.get_query_length();
            }
            cerr << node_name << " r:" << interval.first << ',' << interval.second << ' ' << cigar_code_to_char[intersection.code] << ',' << intersection.length << ',' << length << " r:" << intersection.ref_start << ',' << intersection.ref_stop << " q:" << intersection.query_start << ',' << intersection.query_stop << '\n';

            // If this is a new alignment for this ref/query, inform the GafSummary to initialize a new block
            bool new_query_alignment = prev_node_name.empty();
            bool new_ref_alignment = node_name != prev_node_name;

            // Update the last alignment block in the GafSummary for ref and query
            gaf_summary.update_ref(node_name, i, new_ref_alignment);
            gaf_summary.update_query(node_name, i, new_query_alignment);

            prev_node_name = node_name;
        },{});

        cerr << '\n';

    });
}


string get_name_prefix_of_vcf(const path& vcf){
    string name_prefix = vcf.filename().string();

    if (name_prefix.ends_with(".gz")){
        name_prefix = name_prefix.substr(0,name_prefix.size() - 3);
    }
    if (name_prefix.ends_with(".vcf")){
        name_prefix = name_prefix.substr(0,name_prefix.size() - 4);
    }

    std::replace(name_prefix.begin(), name_prefix.end(), '.', '_');

    return name_prefix;
}


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

    string vcf_name_prefix = get_name_prefix_of_vcf(vcf);

    while (i < regions.size()){
        const auto& region = regions.at(i);

        path input_subdir = output_dir / region.to_string('_');
        path output_subdir = output_dir / region.to_string('_') / vcf_name_prefix;

        string command = "python3 " + vcf_to_gfa_script.string() +
                         " --vcf " + vcf.string() +
                         " --ref " + ref_fasta.string() +
                         " --output_dir " + output_subdir.string() +
                         " --region " + region.to_string();

        run_command(command, true);

        path gaf_path = output_subdir / "alignments.gaf";

        auto name_prefix = get_name_prefix_of_vcf(vcf);

        path gfa_path = output_subdir / (name_prefix + ".gfa");
        path fasta_path = input_subdir / fasta_filename;

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

    // Dump truth haplotypes into each region directory
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

    // Convert VCFs to graphs and run graph aligner
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

    string vcf_name_prefix = get_name_prefix_of_vcf(vcf);
    // Parse the GAFs to get summary info
    {
        for (const auto& region: regions){
            path output_subdir = output_dir / region.to_string('_') / vcf_name_prefix;



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

    cerr << t << "Aligning haplotypes to variant graphs" << '\n';

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
