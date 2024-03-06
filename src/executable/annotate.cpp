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
using std::cerr;
using std::max;
using std::min;
using std::cref;
using std::ref;


using namespace sv_merge;



string get_vcf_name_prefix(const path& vcf){
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
        const int32_t flank_length,
        atomic<size_t>& job_index
){
    size_t i = job_index.fetch_add(1);

    while (i < regions.size()){
        const auto& region = regions.at(i);
        const auto& t = region_transmaps.at(region);

        path output_subdir = output_dir / region.to_unflanked_string('_', flank_length);

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


void get_path_clusters(GafSummary& gaf_summary, const VariantGraph& variant_graph, unordered_map <string,vector<string> >& clusters){
    for (const auto& [name, paths]: gaf_summary.query_paths) {
        // If the query has more than one alignment, dump it into the "unknown" cluster
        if (paths.size() > 1) {
            clusters["unknown"].emplace_back(name);
        }
        // If it has exactly one path, then it can be clustered with all other queries of the same path
        else if (paths.size() == 1) {
            string path_name;
            auto& path = paths.front();

            nid_t id_front = stoll(path.front().first);
            nid_t id_back = stoll(path.back().first);

            // Lord help us
            bool front_is_dangling = variant_graph.is_dangling_node(id_front);
            bool back_is_dangling = variant_graph.is_dangling_node(id_back);

            if (not front_is_dangling or not back_is_dangling){
                path_name = "unknown";
            }
            else {
                // Construct a string identifier for the path (just use GAF convention)
                for (const auto &[node_name, is_reverse]: path) {
                    path_name += (is_reverse ? "<" : ">") + node_name;
                }
            }

            clusters[path_name].emplace_back(name);

        } else {
            throw runtime_error("ERROR: query in gaf summary has no path: " + name);
        }
    }
}


/**
 * Append a log file and write the header if it hasn't been written yet
 * @param output_dir
 * @param vcf_name_prefix
 * @param time_csv result of calling Timer::to_csv() immediately after task exits
 * @param success whether or not the task timed out
 */
void write_time_log(path output_dir, string vcf_name_prefix, string time_csv, bool success){
    // Begin the logging process
    path log_path = output_dir / "log.csv";

    // Check if the log file needs to have a header written to it or not
    bool exists = std::filesystem::exists(log_path);

    ofstream file(log_path, std::ios_base::app);

    // Write the header
    if (not exists){
        file << "name,h,m,s,ms,success" << '\n';
    }
    // Write the results for this region/tool
    file << vcf_name_prefix << ',' << time_csv << ',' << success << '\n';
}


void compute_graph_evaluation_thread_fn(
        unordered_map<Region,vector<VcfRecord> >& region_records,
        const unordered_map<string,vector<interval_t> >& contig_tandems,
        const unordered_map<Region,TransMap>& region_transmaps,
        const unordered_map<string,string>& ref_sequences,
        const vector<Region>& regions,
        const VcfReader& vcf_reader,
        const path& output_dir,
        int32_t flank_length,
        bool cluster,
        atomic<size_t>& job_index
){
    // TODO: finish implementing tandem track as a user input
    unordered_map<string,vector<interval_t>> tandem_track;
    for (const auto& [key,value]: ref_sequences){
        tandem_track[key] = {};
    }

    size_t i = job_index.fetch_add(1);

    path vcf;
    vcf_reader.get_file_path(vcf);
    string vcf_name_prefix = get_vcf_name_prefix(vcf);

    while (i < regions.size()){
        const auto& region = regions.at(i);
        const auto& transmap = region_transmaps.at(region);

        path input_subdir = output_dir / region.to_unflanked_string('_', flank_length);
        path output_subdir = output_dir / region.to_unflanked_string('_', flank_length) / vcf_name_prefix;

        auto records = region_records.at(region);

        create_directories(output_subdir);

        path gfa_path = output_subdir / "graph.gfa";
        path fasta_filename = input_subdir / "haplotypes.fasta";

        VariantGraph variant_graph(ref_sequences, contig_tandems);

        // Check if the region actually contains any usable variants, and use corresponding build() functions
        if (variant_graph.would_graph_be_nontrivial(records)){
            variant_graph.build(records, int32_t(flank_length), numeric_limits<int32_t>::max(), region.start + flank_length, region.stop - flank_length, false);
        }
        else{
            cerr << "TRIVIAL REGION: " + region.to_unflanked_string('_', flank_length) << '\n';
            // VariantGraph assumes that the flank length needs to be added to the region
            variant_graph.build(region.name, region.start + flank_length, region.stop - flank_length, flank_length);
        }

        cerr << "WRITING GFA to file: " << gfa_path << '\n';
        variant_graph.to_gfa(gfa_path);

        path gaf_path = output_subdir / "alignments.gaf";

        auto name_prefix = get_vcf_name_prefix(vcf);

        path fasta_path = input_subdir / fasta_filename;

        string command = "GraphAligner"
                  " -x " "vg"
                  " --multimap-score-fraction " "1"
                  " -t " "1"
                  " -a " + gaf_path.string() +
                  " -g " + gfa_path.string() +
                  " -f " + fasta_path.string();

        // Run GraphAligner and check how long it takes, if it times out
        Timer t;
        bool success = run_command(command, false, 90);
        string time_csv = t.to_csv();

        write_time_log(input_subdir, vcf_name_prefix, time_csv, success);

        // Skip remaining steps for this region/tool if alignment failed and get the next job index for the thread
        if (not success) {
            cerr << "WARNING: Command timed out: " << command << '\n';
            i = job_index.fetch_add(1);
            continue;
        }

        size_t total_coverage = 0;

        // Update variant graph to contain all the paths of the (SPANNING ONLY) alignments
        unordered_map<string,size_t> counter;
        for_alignment_in_gaf(gaf_path, [&](GafAlignment& alignment){
            auto& path = alignment.get_path();

            nid_t id_front = stoll(path.front().first);
            nid_t id_back = stoll(path.back().first);

            // Accumulate total flank coverage (to be averaged later)
            auto l = variant_graph.is_dangling_node(id_front);
            auto r = variant_graph.is_dangling_node(id_back);

            // Skip non-spanning alignments
            if (not l and r){
                return;
            }

            // Fetch (or construct) the count for this read name (default = 0)
            auto& c = counter[alignment.get_query_name()];

            // Create a unique name and add the path to variant graph
            auto name = alignment.get_query_name() + "_" + to_string(c);
            variant_graph.add_gaf_path_to_graph(name, path);

            // Increment counter
            total_coverage++;
            c++;
        });

        path output_path = output_subdir / "annotated.vcf";
        ofstream out_file(output_path);
        if (not out_file.is_open() or not out_file.good()){
            throw runtime_error("ERROR: file could not be written: " + output_path.string());
        }

        out_file << "##INFO=<ID=HAPESTRY_COV,Number=3,Type=.,Description=\"Coverage computed by hapestry of the form [forward_coverage, reverse_coverage, window_coverage]\",Source=\"hapestry\",Version=\"0.0.0.0.0.1\">" << '\n';
        vcf_reader.print_minimal_header(out_file);
        variant_graph.for_each_vcf_record_with_supporting_paths([&](size_t id, const VcfRecord& record, const vector<string>& supporting_paths){
            // Construct a simple summary of read alignments of the form:
            // [forward_coverage, reverse_coverage, window_coverage]
            vector<size_t> coverage = {0,0,total_coverage};
            for (const auto& path_name: supporting_paths){
                auto p = variant_graph.graph.get_path_handle(path_name);
                auto s = variant_graph.graph.path_begin(p);
                auto h = variant_graph.graph.get_handle_of_step(s);

                // Since the read is guaranteed spanning, and the ref is always in F direction,
                // just check the orientation w.r.t. the left flank
                bool is_reverse = variant_graph.graph.get_is_reverse(h);
                coverage[is_reverse]++;
            }

            auto r = record;
            r.info += ";HAPESTRY_COV=" + to_string(coverage[0]) + "," + to_string(coverage[1]) + "," + to_string(coverage[2]);
            r.print(out_file);
            out_file << '\n';
        });

        i = job_index.fetch_add(1);
    }
}


void compute_graph_evaluation(
        const unordered_map <string, interval_tree_t<int32_t> >& contig_interval_trees,
        const unordered_map<string,vector<interval_t> >& contig_tandems,
        const unordered_map<Region,TransMap>& region_transmaps,
        const unordered_map<string,string>& ref_sequences,
        const vector<Region>& regions,
        const path& vcf,
        size_t n_threads,
        int32_t flank_length,
        int32_t interval_max_length,
        bool cluster,
        const path& output_dir
        ){

    unordered_map<Region,vector<VcfRecord> > region_records;
    region_records.reserve(regions.size());

    // Load records for this VCF
    VcfReader vcf_reader(vcf);
    vcf_reader.min_qual = numeric_limits<float>::min();
    vcf_reader.min_sv_length = 1;
    vcf_reader.progress_n_lines = 100'000;
    coord_t record_coord;

    cerr << "Reading VCF... " << '\n';
    vcf_reader.for_record_in_vcf([&](VcfRecord& r){
        // TODO: allow breakends in evaluation
        if (r.sv_type == VcfReader::TYPE_BREAKEND){
            cerr << "WARNING: skipping breakend"  << '\n';
            return;
        }

        r.get_reference_coordinates(false, record_coord);

        // For each overlapping region, put the VcfRecord in that region
        contig_interval_trees.at(r.chrom).overlap_find_all({record_coord.first, record_coord.second}, [&](auto iter){
            coord_t unflanked_window = {iter->low() + flank_length, iter->high() - flank_length};

            // Skip large events in the population
            // TODO: address these as breakpoints in the VariantGraph and avoid constructing windows as intervals
            // for very large events
            if (record_coord.second - record_coord.first > interval_max_length){
                return true;
            }

            // Check if this record exceeds the region
            if (record_coord.first < unflanked_window.first or record_coord.second > unflanked_window.second){
                cerr << "WARNING: skipping record that exceeds the un-flanked window. Record: " << record_coord.first << ',' << record_coord.second << " window: " << unflanked_window.first << ',' << unflanked_window.second << '\n';
                return true;
            }

            Region region(r.chrom, iter->low(), iter->high());
            region_records[region].emplace_back(r);
            return true;
        });
    });

    // Before moving on, make sure every region has at least an empty vector
    for (const auto& r: regions){
        if (region_records.find(r) == region_records.end()){
            region_records[r] = {};
        }
    }

    // Convert VCFs to graphs and run graph aligner
    {
        // Thread-related variables
        atomic<size_t> job_index = 0;
        vector<thread> threads;

        threads.reserve(n_threads);

        // Launch threads
        for (size_t n=0; n<n_threads; n++) {
            try {
                cerr << "launching: " << n << '\n';
                threads.emplace_back(compute_graph_evaluation_thread_fn,
                                     std::ref(region_records),
                                     std::cref(contig_tandems),
                                     std::cref(region_transmaps),
                                     std::cref(ref_sequences),
                                     std::cref(regions),
                                     std::cref(vcf_reader),
                                     std::cref(output_dir),
                                     flank_length,
                                     cluster,
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


void annotate(
        vector<path>& vcfs,
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

    output_dir = std::filesystem::weakly_canonical(output_dir);
    tandem_bed = std::filesystem::weakly_canonical(tandem_bed);
    bam_csv = std::filesystem::weakly_canonical(bam_csv);
    ref_fasta = std::filesystem::weakly_canonical(ref_fasta);

    for (auto& v: vcfs){
        v = std::filesystem::weakly_canonical(v);
    }

    if (std::filesystem::exists(output_dir)){
        throw runtime_error("ERROR: output dir exists already: " + output_dir.string());
    }
    else{
        std::filesystem::create_directories(output_dir);
    }

    cerr << t << "Initializing" << '\n';

    vector <pair <string,path> > bam_paths;
    unordered_map<string,string> ref_sequences;
    vector<Region> regions;

    cerr << t << "Loading reference sequences" << '\n';

    // Load all chromosome sequences (in case of BND)
    for_sequence_in_fasta_file(ref_fasta, [&](const Sequence& s){
        ref_sequences[s.name] = s.sequence;
    });

    cerr << "Reading tandem BED" << '\n';

    unordered_map<string,vector<interval_t> > contig_tandems;
    interval_t interval;
    for_region_in_bed_file(tandem_bed, [&](const Region& r){
        interval.first = r.start;
        interval.second = r.stop;
        contig_tandems[r.name].emplace_back(interval);
    });

    if (windows_bed.empty()){
        cerr << t << "Constructing windows from VCFs and tandem BED" << '\n';
        path bed_log_path = output_dir / "windows_omitted.bed";
        construct_windows_from_vcf_and_bed(ref_sequences, contig_tandems, vcfs, flank_length, interval_max_length, regions, bed_log_path);
    }
    else {
        cerr << t << "Reading BED file" << '\n';
        load_windows_from_bed(windows_bed, regions);
    }

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
                false
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
                force_unique_reads,
                false
        );
    }

    cerr << t << "Aligning sequences to variant graphs" << '\n';

    path fasta_filename = "haplotypes.fasta";
    path staging_dir = output_dir;

    // Dump truth haplotypes into each region directory
    {
        // Thread-related variables
        atomic<size_t> job_index = 0;
        vector<thread> threads;

        threads.reserve(n_threads);

        // Launch threads
        for (size_t n=0; n<n_threads; n++) {
            try {
                cerr << "launching: " << n << '\n';
                threads.emplace_back(write_region_subsequences_to_file_thread_fn,
                                     std::cref(region_transmaps),
                                     std::cref(regions),
                                     std::cref(staging_dir),
                                     std::cref(fasta_filename),
                                     flank_length,
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

    if (debug){
        sleep_for(seconds(60));
        throw runtime_error("DEBUG EARLY EXIT");
    }

    // Generate GFAs/GAFs/CSVs and folder structure for every VCF * every region
    // By default, all of these files will be stored in /dev/shm and then copied into the output dir as a final step.
    // TODO: create option to use /dev/shm/ as staging dir
    // Absolutely must delete the /dev/shm/ copy or warn the user at termination
    //
    for (const auto& vcf: vcfs){
        cerr << "Generating graph alignments for VCF: " << vcf << '\n';

        compute_graph_evaluation(
                contig_interval_trees,
                contig_tandems,
                region_transmaps,
                ref_sequences,
                regions,
                vcf,
                n_threads,
                flank_length,
                interval_max_length,
                false,
                staging_dir
        );
    }

    cerr << t << "Done" << '\n';
}


int main (int argc, char* argv[]){
    path output_dir;
    path windows_bed;
    path tandem_bed;
    string bam_csv;
    path ref_fasta;
    string vcfs_string;
    int32_t flank_length = 150;
    int32_t interval_max_length = 15000;
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
            "--tandems",
            tandem_bed,
            "Path to BED file containing tandem track which will inform how to aggregate variants in windows");

    app.add_option(
            "--vcf",
            vcfs_string,
            "List of VCFs to evaluate (space-separated)")
            ->required()
            ->expected(1,1)
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

    app.add_option(
            "--interval_max_length",
            interval_max_length,
            "How much flanking sequence to use when fetching and aligning reads")
            ->required();

    app.add_flag("--debug", debug, "Invoke this to add more logging and output");

    app.add_flag("--force_unique_reads", force_unique_reads, "Invoke this to add append each read name with the sample name so that inter-sample read collisions cannot occur");

    app.add_flag("--bam_not_hardclipped", bam_not_hardclipped, "Invoke this if you expect your BAMs NOT to contain ANY hardclips. Saves time on iterating.");

    app.parse(argc, argv);

    auto vcfs = app.get_option("--vcf")->as<std::vector<path> >();

    annotate(
        vcfs,
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
