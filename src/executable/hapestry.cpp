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


void cross_align_sample_reads(TransMap& transmap, int64_t score_threshold, const string& sample_name){
    WFAlignerGapAffine aligner(4,6,2,WFAligner::Alignment,WFAligner::MemoryHigh);

    vector <tuple <int64_t,int64_t,float> > edges_to_add;

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

                // Avoid adding while iterating
                edges_to_add.emplace_back(a,b,scaled_score);

                cerr << transmap.get_node(a).name << ',' << transmap.get_node(b).name << ' ' << int64_t(seq_a.size()) << ',' << int64_t(seq_b.size()) << ' ' << score << ' ' << scaled_score << '\n';
            }
        });
    });

    for (const auto& [a,b,score]: edges_to_add){
        transmap.add_edge(a,b,score);
    }
}


void align_reads_vs_paths(TransMap& transmap, const VariantGraph& variant_graph, float min_read_hap_identity){
    // TODO: test these params?? why is gap extension greater than mismatch cost?
    WFAlignerGapAffine aligner(4,6,2,WFAligner::Alignment,WFAligner::MemoryHigh);

    // First extract the sequences of the paths
    unordered_map<string,string> path_sequences;

    transmap.for_each_path([&](const string& name, int64_t id){
        vector <pair <string,bool> > path;
        GafAlignment::parse_string_as_path(name, path);

        string& sequence = path_sequences[name];

        // Get the corresponding bdsg handles from the bdsg handlegraph for each step in the path and fetch + concatenate their sequences
        for (const auto& [node_name, reversal]: path){
            nid_t n = stoll(node_name);
            auto h = variant_graph.graph.get_handle(n, reversal);

            sequence += variant_graph.graph.get_sequence(h);
        }
    });

    vector <tuple <int64_t,int64_t,float> > edges_to_add;

    transmap.for_each_read([&](const string& name, int64_t id){
        auto& read_sequence = transmap.get_sequence(id);

        for (const auto& [path_name, path_sequence]: path_sequences) {
            if (read_sequence.empty() or path_sequence.empty()){
                continue;
            }

            // TODO: Add length cutoff before aligning
            aligner.alignEnd2End(read_sequence, path_sequence);

            // Extract indel edit distance from cigar
            string cigar = aligner.getCIGAR(false);
            string length_token;

            // TODO: switch to using mismatches if we ever put small vars in the graph
            int64_t n_indels = 0;

            for (auto c: cigar){
                if (isalpha(c) or c == '='){
                    if (c == 'I' or c == 'D'){
                        n_indels += stoll(length_token);
                    }

                    length_token.clear();
                }
                else{
                    length_token += c;
                }
            }

            // Can be >1 because indels can be greater than the length of the smaller sequence
            float scaled_score = float(n_indels) / float(min(read_sequence.size(), path_sequence.size()));

//            cerr << "\tname: " << name << "\tpath_name: "  << path_name << "\tscaled_score: "  << scaled_score << "\tn_indels: "  << n_indels << "\tread_sequence: "  << read_sequence.size() << "\tpath_sequence: "  << path_sequence.size() << '\n';

            if ((1.0f - scaled_score) > min_read_hap_identity) {
                // Avoid adding while iterating
                edges_to_add.emplace_back(id, transmap.get_id(path_name), n_indels);
            }
        }
    });

    for (const auto& [a,b,score]: edges_to_add){
        transmap.add_edge(a,b,score);
    }
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
            auto& s = t.get_sequence(id);
            if (s.empty()){
                return;
            }

            file << '>' << name << '\n';
            file << t.get_sequence(id) << '\n';
        });

        i = job_index.fetch_add(1);
    }
}


void get_path_coverages(path gaf_path, const VariantGraph& variant_graph, unordered_map<string,int64_t>& path_coverage){
    // Accumulate coverage of spanning paths in the GAF
    for_alignment_in_gaf(gaf_path, [&](GafAlignment& alignment){
        auto& path = alignment.get_path();
        if (path.size() < 2){
            return;
        }

        nid_t id_front = stoll(path.front().first);
        nid_t id_back = stoll(path.back().first);

        // Check if spanning
        auto l = variant_graph.is_dangling_node(id_front);
        auto r = variant_graph.is_dangling_node(id_back);
        bool is_spanning = l and r;

        if (not is_spanning){
            return;
        }

        auto path_name = alignment.get_path_string();

        path_coverage[path_name]++;
    });
}


void write_solution_to_vcf(
        VariantGraph& variant_graph,
        const TransMap& transmap,
        const vector<string>& sample_names,
        const path& output_path
        ){

    vector <pair <string,int64_t> > samples;

    // Collect all the samples in the window and sort them
    for (const auto& name: sample_names){
        auto sample_id = transmap.get_id(name);
        samples.emplace_back(name, sample_id);
    }

    // TODO: consider just directly overwriting the vector<string> in the VcfRecords stored by VariantGraph
    // Generate a vector of genotypes for each sample, where the vector indexes correspond to the VariantGraph indexes
    unordered_map <string, vector <array<int8_t,2> > > sample_genotypes;

    cerr << output_path << '\n';

    // TODO: find a way to deal with regions for which the transmap was cleared by the optimizer and no calls were made
    // (for now they get 0|0 genotypes)
    for (const auto& [sample_name, sample_id]: samples){
        // Initialize the vectors with arrays of {0,0}
        sample_genotypes[sample_name] = vector <array<int8_t,2> >(variant_graph.vcf_records.size(), {0,0});

        transmap.for_each_phased_variant_of_sample(sample_id, [&](const string& var_name, int64_t _, bool phase){
            // Reconstruct the variant ID from the name
            size_t v = stoull(var_name.substr(1));

            // If there is a transitive edge from sample->path->variant, set the genotype to 1 (with phase)
            sample_genotypes[sample_name][v][phase] = 1;
        });
    }

    // Open the VCF output file
    ofstream vcf_file(output_path);

    if (not (vcf_file.is_open() and vcf_file.good())){
        throw runtime_error("ERROR: could not write file: " + output_path.string());
    }

    // Write the VCF header
    vcf_file << "##fileformat=VCFv4.2" << '\n';
    vcf_file << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (const auto& [sample_name, sample_id]: samples) {
        vcf_file << '\t' << sample_name;
    }
    vcf_file << '\n';

    for (size_t v=0; v<variant_graph.vcf_records.size(); v++){
        auto& record = variant_graph.vcf_records.at(v);

        record.genotypes.clear();
        record.format = "GT";

        // Iterate all the samples and accumulate GTs for the variant object
        for (const auto& [sample_name, sample_id]: samples) {
            string variant_name = "v" + to_string(v);

            // Get the genotype for this sample
            auto& genotype = sample_genotypes.at(sample_name).at(v);

            // Update the record genotypes
            record.genotypes.emplace_back(to_string(int(genotype[0])) + "|" + to_string(int(genotype[1])));
        }

        // Write the record to the VCF
        record.print(vcf_file);
        vcf_file << '\n';
    }
}


void merge_thread_fn(
        unordered_map<Region,vector<VcfRecord> >& region_records,
        const unordered_map<string,vector<interval_t> >& contig_tandems,
        unordered_map<Region,TransMap>& region_transmaps,
        const unordered_map<string,string>& ref_sequences,
        const vector<Region>& regions,
        const VcfReader& vcf_reader,
        const path& output_dir,
        int32_t flank_length,
        int32_t graphaligner_timeout,
        float min_read_hap_identity,
        bool skip_solve,
        atomic<size_t>& job_index
) {
    // TODO: make this a parameter
    int64_t min_path_coverage = 2;

    unordered_map<string, interval_tree_t<int32_t> > contig_interval_trees;

    for (const auto &[contig, intervals]: contig_tandems) {
        for (const auto &interval: intervals) {
            contig_interval_trees[contig].insert({interval.first, interval.second});
        }
    }

    size_t i = job_index.fetch_add(1);

    path vcf;
    vcf_reader.get_file_path(vcf);
    string vcf_name_prefix = get_vcf_name_prefix(vcf);

    // Thread jobs are regions
    while (i < regions.size()) {
        const auto &region = regions.at(i);
        auto& transmap = region_transmaps.at(region);

        path subdir = output_dir / region.to_unflanked_string('_', flank_length);

        auto records = region_records.at(region);

        create_directories(subdir);

        path gfa_path = subdir / "graph.gfa";
        path fasta_filename = subdir / "sequences.fasta";

        VariantGraph variant_graph(ref_sequences, contig_tandems);

        // Check if the region actually contains any usable variants, and use corresponding build() functions
        if (variant_graph.would_graph_be_nontrivial(records)) {
            variant_graph.build(records, int32_t(flank_length), numeric_limits<int32_t>::max(),
                                region.start + flank_length, region.stop - flank_length, false);
        } else {
            // VariantGraph assumes that the flank length needs to be added to the region
            variant_graph.build(region.name, region.start + flank_length, region.stop - flank_length, flank_length);
        }

        cerr << "WRITING GFA to file: " << gfa_path << '\n';
        variant_graph.to_gfa(gfa_path);

        path gaf_path = subdir / "alignments.gaf";

        auto name_prefix = get_vcf_name_prefix(vcf);

        path fasta_path = subdir / fasta_filename;

        string command = "GraphAligner"
                         " --seeds-mum-count " "-1"
                         " --seeds-mxm-windowsize " "0"
                         " -b " "10"
                         " --max-cluster-extend " "10"
                         " --multimap-score-fraction  " "1"
                         " -t " "1"
                         " -a " + gaf_path.string() +
                         " -g " + gfa_path.string() +
                         " -f " + fasta_path.string();

        // Run GraphAligner and check how long it takes, if it times out
        Timer t;
        bool success = run_command(command, false, float(graphaligner_timeout));
        string time_csv = t.to_csv();

        write_graphaligner_time_log(subdir, vcf_name_prefix, time_csv, success);

        // Skip remaining steps for this region/tool if alignment failed and get the next job index for the thread
        if (not success) {
            cerr << "WARNING: Command timed out: " << command << '\n';
            i = job_index.fetch_add(1);
            continue;
        }

        // Iterate the alignments and accumulate their coverages
        unordered_map<string,int64_t> path_coverage;
        get_path_coverages(gaf_path, variant_graph, path_coverage);

        // Add the paths to the TransMap
        for (const auto& [path_name, coverage]: path_coverage) {
            if (coverage >= min_path_coverage){
                transmap.add_path(path_name);
            }
        }

        // Align all reads to all candidate paths and update transmap
        align_reads_vs_paths(transmap, variant_graph, min_read_hap_identity);

        // Write the full transmap to CSV (in the form of annotated edges)
        path output_csv = subdir / "reads_to_paths.csv";
        transmap.write_edge_info_to_csv(output_csv, variant_graph);

        if (not skip_solve){
            try {
                // Optimize
                SolverType solver_type = SolverType::kGscip;
                optimize_reads_with_d_and_n(transmap, 1, 1, 1, subdir, solver_type);

                // Add all the variant nodes to the transmap using a simple name based on the variantgraph ID which likely does
                // not conflict with existing names
                for (size_t v=0; v<variant_graph.vcf_records.size(); v++){
                    transmap.add_variant("v" + to_string(v));
                }

                // Construct mapping of paths to variants
                transmap.for_each_path([&](const string& path_name, int64_t path_id){
                    // Convert the path to a vector
                    vector <pair <string,bool> > path;

                    GafAlignment::parse_string_as_path(path_name, path);

                    variant_graph.for_each_vcf_record(path, [&](size_t id, const vector<edge_t>& edges_of_the_record, const VcfRecord& record){
                        string var_name = "v" + to_string(id);
                        transmap.add_edge(var_name, path_name);
                    });
                });

                // Write the solution to a VCF
                path output_path = subdir / "solution.vcf";

                write_solution_to_vcf(variant_graph, transmap, vcf_reader.sample_ids, output_path);
            }
            catch (const exception& e) {
                cerr << e.what() << '\n';
                cerr << "ERROR caught at " << region.to_string() << '\n';
            }
        }

        i = job_index.fetch_add(1);
    }
}


void merge_variants(
        const unordered_map <string, interval_tree_t<int32_t> >& contig_interval_trees,
        const unordered_map<string,vector<interval_t> >& contig_tandems,
        unordered_map<Region,TransMap>& region_transmaps,
        const unordered_map<string,string>& ref_sequences,
        const vector<Region>& regions,
        const path& vcf,
        size_t n_threads,
        int32_t flank_length,
        int32_t interval_max_length,
        int32_t min_sv_length,
        int32_t graphaligner_timeout,
        float min_read_hap_identity,
        bool skip_solve,
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
                threads.emplace_back(merge_thread_fn,
                                     std::ref(region_records),
                                     std::cref(contig_tandems),
                                     std::ref(region_transmaps),
                                     std::cref(ref_sequences),
                                     std::cref(regions),
                                     std::cref(vcf_reader),
                                     std::cref(output_dir),
                                     flank_length,
                                     graphaligner_timeout,
                                     min_read_hap_identity,
                                     skip_solve,
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


void hapestry(
        path vcf,
        path output_dir,
        path windows_bed,                // Override the interval graph if this is provided
        path tandem_bed,
        path bam_csv,
        path ref_fasta,
        int32_t flank_length,
        int32_t interval_max_length,
        int32_t min_sv_length,
        int32_t n_threads,
        int32_t graphaligner_timeout,
        float min_read_hap_identity,
        bool skip_solve,
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
        construct_windows_from_vcf_and_bed(ref_sequences, contig_tandems, {vcf}, flank_length, interval_max_length, min_sv_length, regions, bed_log_path, false);
    }
    else {
        cerr << t << "Reading BED file" << '\n';
        load_windows_from_bed(windows_bed, regions);
    }

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
                    flank_length
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
                    flank_length
            );
        }

        cerr << t << "Writing sequences to disk" << '\n';

        path fasta_filename = "sequences.fasta";
        path staging_dir = output_dir;

        // Dump sequences into each region directory
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
        cerr << "Generating graph alignments for VCF: " << vcf << '\n';

        merge_variants(
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
                min_read_hap_identity,
                skip_solve,
                staging_dir
        );

        auto vcf_prefix = get_vcf_name_prefix(vcf);
        path out_vcf = output_dir / ("merged.vcf");
        ofstream out_file(out_vcf);

        if (not (out_file.is_open() and out_file.good())){
            throw runtime_error("ERROR: could not write file: " + out_vcf.string());
        }

        ifstream input_vcf(vcf);

        if (not (input_vcf.is_open() and input_vcf.good())){
            throw runtime_error("ERROR: could not write file: " + vcf.string());
        }

        // Just copy over the header lines
        string line;
        while (getline(input_vcf, line)){
            if (line.starts_with("##")){
                out_file << line << '\n';
            }
            else{
                input_vcf.close();
                break;
            }
        }

        // Copy over the mutable parts of the header and then the main contents of the filtered VCF
        // Will be empty if no regions were processed
        for (size_t i=0; i<regions.size(); i++){
            const auto& region = regions[i];
            path sub_vcf = output_dir / region.to_unflanked_string('_', flank_length) / "solution.vcf";

            ifstream file(sub_vcf);
            while (getline(file, line)){
                if (line.starts_with('#')){
                    if (not line.starts_with("##fileformat")){
                        if (i == 0){
                            out_file << line << '\n';
                        }
                    }
                }
                else{
                    out_file << line << '\n';
                }
            }
        }
    }
    else{
        cerr << "WARNING: No regions to process" << '\n';
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
    int32_t graphaligner_timeout = 90;
    float min_read_hap_identity = 0.85;
    bool skip_solve = false;
    bool debug = false;
    bool force_unique_reads = false;
    bool bam_not_hardclipped = false;

    CLI::App app{"App description"};

    // TODO: add skip_solve arg
    // TODO: add graphaligner timeout arg
    // TODO: break out transmap edge filtering and add edge weight parameter

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

    app.add_option(
            "--interval_max_length",
            interval_max_length,
            "How long a window can be in bp before it is skipped")
            ->required();

    app.add_option(
            "--min_sv_length",
            min_sv_length,
            "Skip all variants less than this length (bp)")
            ->required();

    app.add_option(
            "--graphaligner_timeout",
            graphaligner_timeout,
            "Abort graphaligner after this many seconds, and do not compute the remaining steps for that window")
            ->required();

    app.add_option(
            "--min_read_hap_identity",
            min_read_hap_identity,
            "Minimum alignment identity to consider a read-hap assignment in the optimization step")
            ->required();

    app.add_flag("--skip_solve", skip_solve, "Invoke this to skip the optimization step. CSVs for each optimization input will still be written.");

    app.add_flag("--debug", debug, "Invoke this to add more logging and output");

    app.add_flag("--force_unique_reads", force_unique_reads, "Invoke this to add append each read name with the sample name so that inter-sample read collisions cannot occur");

    app.add_flag("--bam_not_hardclipped", bam_not_hardclipped, "Invoke this if you expect your BAMs NOT to contain ANY hardclips. Saves time on iterating.");

    try{
        app.parse(argc, argv);
    }
    catch (const CLI::ParseError &e) {
        return app.exit(e);
    }

    hapestry(
            vcf,
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
            min_read_hap_identity,
            skip_solve,
            debug,
            force_unique_reads,
            bam_not_hardclipped
    );

    return 0;
}
