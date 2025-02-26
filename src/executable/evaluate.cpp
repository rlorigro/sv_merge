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
        const int32_t flank_length,
        atomic<size_t>& job_index
){
    size_t i = job_index.fetch_add(1);
    string s;

    while (i < regions.size()){
        const auto& region = regions.at(i);
        const auto& t = region_transmaps.at(region);

        path output_subdir = output_dir / region.to_unflanked_string('_', flank_length);

        create_directories(output_subdir);

        path output_fasta = output_subdir / filename;
        ofstream file(output_fasta);

        t.for_each_read([&](const string& name, int64_t id){
            file << '>' << name << '\n';

            t.get_sequence(id,s);

            file << s << '\n';
        });

        i = job_index.fetch_add(1);
    }
}


/**
 * @param output_dir
 * @param gaf_summary summary object which has been updated during the iteration of alignments
 * @param variant_graph cannot be const, but is effectively const in this context. This object is used to determine
 * if alignments cover variants in the graph
 */
void write_summary(path output_dir, const GafSummary& gaf_summary, VariantGraph& variant_graph, const VcfReader& vcf_reader){
    path nodes_output_path = output_dir / "nodes.csv";
    path edges_output_path = output_dir / "edges.csv";
    path haps_output_path = output_dir / "haps.csv";
    path supported_output_path = output_dir / "supported.vcf";
    path unsupported_output_path = output_dir / "unsupported.vcf";

    ofstream nodes_file(nodes_output_path);
    if (not nodes_file.is_open() or not nodes_file.good()){
        throw runtime_error("ERROR: file could not be written: " + nodes_output_path.string());
    }

    ofstream haps_file(haps_output_path);
    if (not haps_file.is_open() or not haps_file.good()){
        throw runtime_error("ERROR: file could not be written: " + haps_output_path.string());
    }

    // Start by writing the headers
    nodes_file << "name,length,is_ref,is_flank,coverage,identity,color" << '\n';
    haps_file << "name,length,is_ref,is_flank,coverage,identity" << '\n';

    string ref_color = "#7D8FE1";
    string non_ref_color = "#FFC07E";
    string flank_color = "#bebebe";   // bebe be bebe
    string color;

    // GAF ref/target nodes
    // Write a CSV file with the format:
    // name,length,is_ref,coverage,identity,color
    // [string],[int],[bool],[float],[float],[string]
    gaf_summary.for_each_ref_summary([&](const string& name, int32_t length, float identity, float coverage){
        auto id = stoll(name);

        bool is_ref = variant_graph.is_reference_node(id);
        bool is_flank = variant_graph.is_flanking_node(id);

        color = non_ref_color;
        if (is_ref){
            color = ref_color;
        }
        if (is_flank){
            color = flank_color;
        }

        nodes_file << name << ',' << length << ',' << is_ref << ',' << is_flank << ',' << coverage << ',' << identity << ',' << color << '\n';
    });

    // GAF queries (aligned haplotypes)
    // Write a CSV file with the format:
    // name,length,is_ref,coverage,identity,color
    // [string],[int],[bool],[float],[float],[string]
    gaf_summary.for_each_query_summary([&](const string& name, int32_t length, float identity, float coverage){
        haps_file << name << ',' << length << ',' << 0 << ',' << 0 << ',' << coverage << ',' << identity << '\n';
    });

    // Iterate the edges covered by the alignments and compile some stats regarding edge coverage
    size_t n_edges = 0;
    size_t n_non_ref_edges = 0;
    size_t n_non_ref_edges_covered = 0;

    unordered_set <pair <handle_t, handle_t> > covered_edges;

    // Walk through all the paths that each query had
    for (const auto& [name, paths]: gaf_summary.query_paths){
        // Iterate the paths (this query may have multiple alignments)
        for (size_t j=0; j<paths.size(); j++){
            const auto& path = paths[j];
            auto alignment_name = name + "_" + to_string(j);

            // Create a new path in the variant graph
            auto p = variant_graph.graph.create_path_handle(alignment_name);

            // Iterate the path steps and append each step to the prev step
            handle_t h_prev;
            for (size_t i=0; i<path.size(); i++){
                const auto& [step_name, is_reverse] = path[i];

                // Convert GAF name string back into nid, and construct handle of correct orientation
                nid_t id = stoll(step_name);
                auto h = variant_graph.graph.get_handle(id,is_reverse);
                variant_graph.graph.append_step(p,h);

                if (i > 0){
                    auto canonical_edge = variant_graph.graph.edge_handle(h_prev,h);
                    covered_edges.emplace(canonical_edge);
                }
                h_prev = h;
            }
        }
    }

    ofstream edges_file(edges_output_path);
    if (not edges_file.is_open() or not edges_file.good()){
        throw runtime_error("ERROR: file could not be written: " + edges_output_path.string());
    }

    variant_graph.graph.for_each_edge([&](const edge_t e){
        // We want to skip ref-ref edges in the evaluation
        auto is_ref = variant_graph.is_reference_edge(e);

        if (not is_ref) {
            n_non_ref_edges += 1;

            if (covered_edges.find(e) != covered_edges.end()){
                n_non_ref_edges_covered += 1;
            }
        }

        n_edges += 1;
    });

    size_t n_alignments = 0;
    for (const auto& [name, paths]: gaf_summary.query_paths){
        n_alignments += paths.size();
    }

    // Write a CSV of the format
    // n_alignments,n_edges,n_edges_covered,n_nonref_edges,n_nonref_edges_covered
    // [int],[int],[int],[int],[int]
    edges_file << "n_alignments" << ',' << "n_edges" << ',' << "n_edges_covered" << ',' << "n_non_ref_edges" << ',' << "n_non_ref_edges_covered" << '\n';
    edges_file << n_alignments << ',' << n_edges << ',' << covered_edges.size() << ',' << n_non_ref_edges << ',' << n_non_ref_edges_covered << '\n';

    // Update the VariantGraph paths to contain all the alignments
    ofstream supported_file(supported_output_path);
    if (not supported_file.is_open() or not supported_file.good()){
        throw runtime_error("ERROR: file could not be written: " + supported_output_path.string());
    }

    ofstream unsupported_file(unsupported_output_path);
    if (not unsupported_file.is_open() or not unsupported_file.good()){
        throw runtime_error("ERROR: file could not be written: " + unsupported_output_path.string());
    }

    vcf_reader.print_minimal_header(supported_file);
    vcf_reader.print_minimal_header(unsupported_file);
    variant_graph.print_supported_vcf_records(supported_file, unsupported_file, false);
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


void compute_graph_evaluation_thread_fn(
        unordered_map<Region,vector<VcfRecord> >& region_records,
        const unordered_map<string,vector<interval_t> >& contig_tandems,
        const unordered_map<Region,TransMap>& region_transmaps,
        const unordered_map<string,string>& ref_sequences,
        const vector<Region>& regions,
        const VcfReader& vcf_reader,
        const path& output_dir,
        int32_t flank_length,
        size_t graphaligner_timeout,
        bool cluster,
        atomic<size_t>& job_index
){
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
        bool success = run_command(command, false, graphaligner_timeout);
        string time_csv = t.to_csv();

        write_time_log(input_subdir, vcf_name_prefix, time_csv, success);

        // Skip remaining steps for this region/tool if alignment failed and get the next job index for the thread
        if (not success) {
            cerr << "WARNING: Command timed out: " << command << '\n';
            i = job_index.fetch_add(1);
            continue;
        }

        // Do the bulk of the Gaf parsing work here
        GafSummary gaf_summary(variant_graph, transmap, true);

        // Allow only this many bases to be parsed within the flanking region to accumulate cigar ops that are
        // "squished" out of the window
        gaf_summary.flank_buffer = 40;
        gaf_summary.compute(gaf_path);

        // Write out all the alignment dependent results for this region
        write_summary(output_subdir, gaf_summary, variant_graph, vcf_reader);

        // Write some additional info about the "cluster" that each haplotype belongs to in case it is desired and
        // this vcf happens to be the one that was assigned as the one to cluster by.
        // In this case each "cluster" is just a non-redundant path observed in the GAF,
        if (cluster){
            unordered_map <string,vector<string> > clusters;
            get_path_clusters(gaf_summary, variant_graph, clusters);

            path clusters_path = input_subdir / "clusters.csv";
            ofstream file(clusters_path);

            if (not file.is_open() or not file.good()){
                throw runtime_error("ERROR: could not write file: " + clusters_path.string());
            }

            file << "cluster_name,query_names" << '\n';
            for (const auto& [cluster_name, query_names]: clusters){
                file << cluster_name << ',';
                for (size_t q=0; q<query_names.size(); q++){
                    file << query_names[q] << ((q < query_names.size() - 1) ? ' ' : '\n');
                }
            }
        }

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
        int32_t min_sv_length,
        size_t graphaligner_timeout,
        bool cluster,
        const path& output_dir
        ){

    unordered_map<Region,vector<VcfRecord> > region_records;
    region_records.reserve(regions.size());

    // Load records for this VCF
    VcfReader vcf_reader(vcf);
    vcf_reader.min_qual = numeric_limits<float>::min();
    vcf_reader.min_sv_length = min_sv_length;
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

        auto result = contig_interval_trees.find(r.chrom);

        // First make sure there are actually some windows in this contig (might fail with small satellite contigs)
        if (result == contig_interval_trees.end()){
            return;
        }

        // For each overlapping region, put the VcfRecord in that region
        result->second.overlap_find_all({record_coord.first, record_coord.second}, [&](auto iter){
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
                                     graphaligner_timeout,
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


void evaluate(
        vector<path>& vcfs,
        path cluster_by,
        path output_dir,
        path windows_bed,                // Override the interval graph if this is provided
        path tandem_bed,
        path bam_csv,
        path ref_fasta,
        int32_t flank_length,
        int32_t interval_max_length,
        int32_t min_sv_length,
        int32_t n_threads,
        size_t graphaligner_timeout,
        bool debug,
        bool force_unique_reads
){
    Timer t;

    output_dir = std::filesystem::weakly_canonical(output_dir);
    tandem_bed = std::filesystem::weakly_canonical(tandem_bed);
    bam_csv = std::filesystem::weakly_canonical(bam_csv);
    ref_fasta = std::filesystem::weakly_canonical(ref_fasta);
    cluster_by = std::filesystem::weakly_canonical(cluster_by);

    for (auto& v: vcfs){
        v = std::filesystem::weakly_canonical(v);
    }

    if (std::filesystem::exists(output_dir)){
        throw runtime_error("ERROR: output dir exists already: " + output_dir.string());
    }
    else{
        std::filesystem::create_directories(output_dir);
    }

    if (std::find(vcfs.begin(), vcfs.end(), cluster_by) == vcfs.end()){
        for (const auto& v: vcfs){
            cerr << v << '\n';
        }
        throw runtime_error("ERROR: --cluster_by parameter must match one of the paths provided by --vcfs");
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
        construct_windows_from_vcf_and_bed(ref_sequences, contig_tandems, vcfs, flank_length, interval_max_length, min_sv_length, regions, bed_log_path);
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

    fetch_reads_from_clipped_bam(
            t,
            regions,
            bam_csv,
            n_threads,
            max_length,
            flank_length,
            region_transmaps,
            true,
            true,
            true,
            force_unique_reads,
            true,
            flank_length
    );

    cerr << t << "Aligning haplotypes to variant graphs" << '\n';

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

    // Generate GFAs/GAFs/CSVs and folder structure for every VCF * every region
    // By default, all of these files will be stored in /dev/shm and then copied into the output dir as a final step.
    // TODO: create option to use /dev/shm/ as staging dir
    // Absolutely must delete the /dev/shm/ copy or warn the user at termination
    //
    for (const auto& vcf: vcfs){
        cerr << "Generating graph alignments for VCF: " << vcf << '\n';

        bool cluster = (cluster_by == vcf);
        if (cluster){
            cerr << "Is cluster VCF" << '\n';
        }

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
                min_sv_length,
                graphaligner_timeout,
                cluster,
                staging_dir
        );
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
    path cluster_by;
    string vcfs_string;
    int32_t flank_length = 250;
    int32_t interval_max_length = 15000;
    int32_t min_sv_length = 1;
    int32_t n_threads = 1;
    size_t graphaligner_timeout = 120;
    bool debug = false;
    bool force_unique_reads = false;

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
            "--vcfs",
            vcfs_string,
            "List of VCFs to evaluate (space-separated)")
            ->required()
            ->expected(1,-1)
            ->delimiter(',');

    app.add_option(
            "--cluster_by",
            cluster_by,
            "Simple headerless CSV file with the format [sample_name],[hap_name],[bam_path]")
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

    app.add_option(
            "--interval_max_length",
            interval_max_length,
            "How much flanking sequence to use when fetching and aligning reads")
            ->required();

    app.add_option(
            "--min_sv_length",
            min_sv_length,
            "Skip all variants less than this length (bp)");

    app.add_option(
            "--graphaligner_timeout",
            graphaligner_timeout,
            "Maximum time to spend on GraphAligner (seconds)");

    app.add_flag("--debug", debug, "Invoke this to add more logging and output");

    app.add_flag("--force_unique_reads", force_unique_reads, "Invoke this to add append each read name with the sample name so that inter-sample read collisions cannot occur");

    app.parse(argc, argv);

    auto vcfs = app.get_option("--vcfs")->as<std::vector<path> >();

    if (debug) {
        HAPESTRY_DEBUG = true;
    }

    evaluate(
        vcfs,
        cluster_by,
        output_dir,
        windows_bed,
        tandem_bed,
        bam_csv,
        ref_fasta,
        flank_length,
        interval_max_length,
        min_sv_length,
        n_threads,
        graphaligner_timeout,
        debug,
        force_unique_reads
    );

    return 0;
}
