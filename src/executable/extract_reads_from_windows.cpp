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

#include "cpptrace/from_current.hpp"


using namespace sv_merge;


void write_region_subsequences_to_fasta(const TransMap& t, const FetchConfig& config, const path& output_fasta) {
    string s;

    ofstream file(output_fasta);

    vector<pair<string,int64_t>> ids;
    ids.reserve(t.get_read_count());

    t.for_each_read([&](const string& name, int64_t id){
        t.get_sequence(id,s);

        if (s.empty()){
            return;
        }

        ids.emplace_back(name,id);
    });

    // Reads will be accessed from an unordered_map<int,...>, but not guaranteed to be shuffled so we enforce shuffling
    std::ranges::shuffle(ids, std::mt19937(1337));

    bool is_reverse;
    string tags;
    for (auto& [name,id]: ids) {
        tags.clear();
        if (not config.tags_to_fetch.empty()) {
            t.get_sequence_tags(id, tags);
        }

        // Sequences are compressed so they cant be fetched as string refs directly
        t.get_sequence(id,s);

        is_reverse = t.get_sequence_reversal(id);

        file << ">" << name << ' ' << (is_reverse ? 'R' : 'F') << (tags.empty() ? "" : " ") << tags << '\n';
        file << s << '\n';
    }
}


void write_region_subsequences_to_fastq(const TransMap& t, const FetchConfig& config, const path& output_fastq) {
    string s;

    ofstream file(output_fastq);

    vector<pair<string,int64_t>> ids;
    ids.reserve(t.get_read_count());

    t.for_each_read([&](const string& name, int64_t id){
        t.get_sequence(id,s);

        if (s.empty()){
            return;
        }

        ids.emplace_back(name,id);
    });

    // Reads will be accessed from an unordered_map<int,...>, but not guaranteed to be shuffled so we enforce shuffling
    std::ranges::shuffle(ids, std::mt19937(1337));

    bool is_reverse;
        string tags;

    for (auto& [name,id]: ids) {
        tags.clear();
        if (not config.tags_to_fetch.empty()) {
            t.get_sequence_tags(id, tags);
        }

        // Sequences are compressed so they cant be fetched as string refs directly
        t.get_sequence(id,s);

        const auto& qualities = t.get_sequence_qualities(id);
        is_reverse = t.get_sequence_reversal(id);

        file << "@" << name << ' ' << (is_reverse ? 'R' : 'F') << (tags.empty() ? "" : " ") << tags << '\n';
        file << s << '\n';
        file << "+" << '\n';
        for (const auto& q: qualities){
            if (q+33 < 33 or q+33 > 126){
                throw runtime_error("ERROR: quality score out of range: " + std::to_string(q+33) + " for read: " + name + " in region: " + output_fastq.string() + " at position: " + std::to_string(q));
            }

            file << char(q+33);
        }
        file << '\n';
    }
}


void write_region_subsequences_to_file_thread_fn(
        const unordered_map<Region,TransMap>& region_transmaps,
        const vector<Region>& regions,
        const path& output_dir,
        const path& filename,
        const FetchConfig& config,
        mutex& err_mutex,
        atomic<size_t>& job_index
){
    size_t i = job_index.fetch_add(1);

    while (i < regions.size()){
        const auto& region = regions.at(i);
        const auto& t = region_transmaps.at(region);

        path output_subdir = output_dir / region.to_unflanked_string('_', config.flank_length);

        create_directories(output_subdir);

        path output_fasta = output_subdir / filename;

        CPPTRACE_TRY {
            if (config.get_qualities) {
                write_region_subsequences_to_fastq(t, config, output_fasta);
            }
            else {
                write_region_subsequences_to_fasta(t, config, output_fasta);
            }
        } CPPTRACE_CATCH(const std::exception& e) {
            err_mutex.lock();

            cerr << "ERROR in thread job " << i << " caught in writing sequence for file: " << filename << '\n';
            cerr << "Exception: " << e.what() << '\n';
            cpptrace::from_current_exception().print_with_snippets();

            err_mutex.unlock();
        }

        i = job_index.fetch_add(1);
    }
}


void extract(
        path output_dir,
        size_t n_threads,
        path bam_csv,
        path bed_path,
        bool bam_not_hardclipped,
        FetchConfig& config
        ){

    if (std::filesystem::exists(output_dir)){
        throw runtime_error("ERROR: output dir exists already: " + output_dir.string());
    }
    else{
        std::filesystem::create_directories(output_dir);
    }

    if (not bam_not_hardclipped and (config.get_qualities or not config.tags_to_fetch.empty())) {
        throw runtime_error("ERROR: tag and quality fetching not implemented for hardclipped BAMs, use a non-"
                            "hardclipped BAM and --bam_not_hardclipped or remove tag and quality parameters."
                            " This feature could be supported in the future. Open an issue if needed.");
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
        r.start = max(1,r.start-config.flank_length);
        r.stop += config.flank_length;

        contig_interval_trees[r.name].insert({r.start, r.stop});
    }

    cerr << t << "Fetching reads for all windows" << '\n';

    Authenticator authenticator;

    // Intermediate object to store results of multithreaded sample read fetching
    sample_region_flanked_coord_map_t sample_to_region_coords;

    // Intermediate objects
    vector <pair <string, path> > sample_bams;

    cerr << t << "Loading CSV" << '\n';

    cerr << t << "Processing windows" << '\n';
    // The container to store all fetched reads and their relationships to samples/paths
    unordered_map<Region,TransMap> region_transmaps;

    config.n_threads = n_threads;

    if (bam_not_hardclipped){
        cerr << "Fetching from NON-hardclipped BAMs" << '\n';

        fetch_reads(
                t,
                regions,
                bam_csv,
                config,
                region_transmaps
        );

    }
    else{
        cerr << "Fetching from HARDCLIPPED BAMs" << '\n';

        config.unclip_coords = true;

        fetch_reads_from_clipped_bam(
                t,
                regions,
                bam_csv,
                config,
                region_transmaps
        );
    }

    cerr << t << "Peak memory usage: " << get_peak_memory_usage() << '\n';
    cerr << t << "Writing sequences to disk" << '\n';
    cerr << t << "Writing sequences" << '\n';

    path output_filename;

    if (config.get_qualities) {
        output_filename = "sequences.fastq";
    }
    else {
        output_filename = "sequences.fasta";
    }

    // Dump sequences into each region directory
    {
        // Thread-related variables
        atomic<size_t> job_index = 0;
        vector<thread> threads;

        threads.reserve(n_threads);

        mutex err_mutex;

        // Launch threads
        for (size_t n=0; n<n_threads; n++) {
            try {
                cerr << "launching: " << n << '\n';
                threads.emplace_back(write_region_subsequences_to_file_thread_fn,
                                     std::cref(region_transmaps),
                                     std::cref(regions),
                                     std::cref(output_dir),
                                     std::cref(output_filename),
                                     std::cref(config),
                                     std::ref(err_mutex),
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
    bool bam_not_hardclipped;
    size_t n_threads;
    path output_dir;
    path windows_bed;
    path bam_csv;
    path bed_path;
    path ref;
    string tags_arg;
    FetchConfig config;

    CLI::App app{"App description"};

    app.add_option(
            "--output_dir",
            output_dir,
            "Path to output directory which must not exist yet")
            ->required();

    app.add_option(
            "--n_threads",
            n_threads,
            "Maximum number of threads to use for fetching. To avoid being throttled by cloud providers.")
            ->required();

    app.add_option(
            "--bam_csv",
            bam_csv,
            "Simple headerless CSV file with the format [sample_name],[bam_path]")
            ->required();

    app.add_option(
            "--windows",
            bed_path,
            "Path to BED file containing windows to extract reads from")
            ->required();

    app.add_option(
            "--flank_length",
            config.flank_length,
            "How much flanking sequence to use when fetching and aligning reads")
            ->required();

    app.add_flag(
            "--require_spanning",
            config.require_spanning,
            "If this flag is invoked, then only reads that span the entire window will be fetched");

    app.add_option(
            "--fetch_max_length",
            config.max_length,
            "How long a sequence within a window can be in bp before it is skipped (important for large contigs with clipping, may be millions bp)")
            ->required();

    app.add_flag(
            "--bam_not_hardclipped",
            bam_not_hardclipped,
            "If this flag is invoked, then only reads that span the entire window will be fetched");

    app.add_flag(
            "--force_forward",
            config.force_forward,
            "If this flag is invoked, reverse complement any reads that are on the reverse strand");

    app.add_flag(
            "--get_qualities",
            config.get_qualities,
            "If this flag is invoked, also fetch the qualities of the reads as a fastq");

    app.add_option(
            "--tags",
            tags_arg,
            "A comma separated list of tags to fetch from the BAM file (e.g. NM,PS,HP) and append to the "
            "fastq name as space-separated fields");

    app.add_flag("--force_unique_reads", config.append_sample_to_read, "Invoke this to add append each read name with the sample name so that inter-sample read collisions cannot occur");

    CLI11_PARSE(app, argc, argv);

    parse_comma_separated_string(tags_arg, config.tags_to_fetch);

    extract(
        output_dir,
        n_threads,
        bam_csv,
        bed_path,
        bam_not_hardclipped,
        config
    );


    return 0;
}
