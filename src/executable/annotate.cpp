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


class VariantSupport{
public:
    // is_spanning -> is_reverse -> Q
    array <array <vector<float>, 2>, 2> identity;
    int32_t length_of_evaluated_region;
    bool is_tandem;

    VariantSupport();
    [[nodiscard]] size_t get_coverage(bool is_spanning, bool is_reverse) const;
    void get_identity_distribution(bool is_spanning, bool is_reverse, array<size_t,6>& result) const;
    void get_support_string(string& s) const;
};


VariantSupport::VariantSupport():
        length_of_evaluated_region(0),
        is_tandem(false)
{}


void VariantSupport::get_support_string(string& s) const{
    s.clear();

    array<size_t,6> d;

//    s += to_string(get_coverage(false)) + ',';
//    s += to_string(get_coverage(true)) + ',';

    get_identity_distribution(false, false, d);
    s += to_string(d[0]) + ',';
    s += to_string(d[1]) + ',';
    s += to_string(d[2]) + ',';
    s += to_string(d[3]) + ',';
    s += to_string(d[4]) + ',';
    s += to_string(d[5]) + ',';

    get_identity_distribution(false, true, d);
    s += to_string(d[0]) + ',';
    s += to_string(d[1]) + ',';
    s += to_string(d[2]) + ',';
    s += to_string(d[3]) + ',';
    s += to_string(d[4]) + ',';
    s += to_string(d[5]) + ',';

    get_identity_distribution(true, false, d);
    s += to_string(d[0]) + ',';
    s += to_string(d[1]) + ',';
    s += to_string(d[2]) + ',';
    s += to_string(d[3]) + ',';
    s += to_string(d[4]) + ',';
    s += to_string(d[5]) + ',';

    get_identity_distribution(true, true, d);
    s += to_string(d[0]) + ',';
    s += to_string(d[1]) + ',';
    s += to_string(d[2]) + ',';
    s += to_string(d[3]) + ',';
    s += to_string(d[4]) + ',';
    s += to_string(d[5]);
}


size_t VariantSupport::get_coverage(bool is_spanning, bool is_reverse) const{
    return identity[is_reverse].size();
}


void VariantSupport::get_identity_distribution(bool is_spanning, bool is_reverse, array<size_t,6>& result) const{
    result.fill(0);
    for (auto i: identity[is_spanning][is_reverse]){
        auto index = int64_t(result.size() - 1);

        if (i < 1.0){
            // q = 0, p(correct) = 0.0  Merged with above
            // q = 1, p(correct) = 0.5  Merged with above
            // q = 2, p(correct) = 0.75
            // q = 3, p(correct) = 0.875
            // q = 4, p(correct) = 0.9375
            // q = 5, p(correct) = 0.96875
            // q = 6, p(correct) = 0.984375
            // q = 7, p(correct) = 0.9921875
            auto q = -log2(1-i);
            // Compute q and make the bottom bin open ended
            q = max(0.0,q-2.0f);

            // Bin the q value into one of the distribution indexes, with the top bin open ended
            index = min(int64_t(result.size()-1), int64_t(floor(q)));

//            cerr << is_spanning << ',' << i << ',' << q << ',' << index << '\n';
        }
        else if (i > 1.0){
            throw runtime_error("ERROR: identity greater than 1.0 encountered: " + to_string(i));
        }

        result[index]++;
    }
}


void compute_graph_evaluation_thread_fn(
        unordered_map<Region,vector<VcfRecord> >& region_records,
        const unordered_map<string,vector<interval_t> >& contig_tandems,
        const unordered_map<Region,TransMap>& region_transmaps,
        const unordered_map<string,string>& ref_sequences,
        const vector<Region>& regions,
        const VcfReader& vcf_reader,
        const string& label,
        const path& output_dir,
        int32_t flank_length,
        atomic<size_t>& job_index
){
    unordered_map <string, interval_tree_t<int32_t> > contig_interval_trees;

    for (const auto& [contig,intervals]: contig_tandems) {
        for (const auto& interval: intervals) {
            contig_interval_trees[contig].insert({interval.first, interval.second});
        }
    }

    size_t i = job_index.fetch_add(1);

    path vcf;
    vcf_reader.get_file_path(vcf);
    string vcf_name_prefix = get_vcf_name_prefix(vcf);
    string buffer;

    while (i < regions.size()){
        const auto& region = regions.at(i);
        path subdir = output_dir / region.to_unflanked_string('_', flank_length);
        auto records = region_records.at(region);
        create_directories(subdir);
        path gfa_path = subdir / "graph.gfa";
        path fasta_filename = subdir / "haplotypes.fasta";
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
        bool success = run_command(command, false, 90);
        string time_csv = t.to_csv();

        write_graphaligner_time_log(subdir, vcf_name_prefix, time_csv, success);

        // Skip remaining steps for this region/tool if alignment failed and get the next job index for the thread
        if (not success) {
            cerr << "WARNING: Command timed out: " << command << '\n';
            i = job_index.fetch_add(1);
            continue;
        }

        const int32_t variant_flank_length = 50;
        size_t l_coverage = 0;
        size_t r_coverage = 0;

        vector<VariantSupport> variant_supports(variant_graph.vcf_records.size());
        vector<float> max_observed_identities(variant_graph.vcf_records.size());
        // First annotate all the variants as tandem or not
        int32_t bnd_pos;
        string bnd_chromosome;
        for (size_t v=0; v<variant_supports.size(); v++){
            auto& record = variant_graph.vcf_records[v];
            coord_t ref_coord;
            record.get_reference_coordinates(false, ref_coord);
            bool is_tandem = false;
            if (contig_interval_trees.contains(region.name)) {
                contig_interval_trees.at(region.name).overlap_find_all({ref_coord.first, ref_coord.second}, [&](auto iter) {
                    is_tandem=true;
                    return false;
                });
            }
            if (!is_tandem && record.sv_type==VcfReader::TYPE_BREAKEND && record.is_breakend_single()==2) {
                record.get_breakend_chromosome(bnd_chromosome);
                bnd_pos=record.get_breakend_pos();
                if (contig_interval_trees.contains(bnd_chromosome)) {
                    contig_interval_trees.at(bnd_chromosome).overlap_find_all({bnd_pos,bnd_pos}, [&](auto iter) {
                        is_tandem=true;
                        return false;
                    });
                }
            }
            variant_supports[v].is_tandem = is_tandem;
        }

        // Iterate the alignments and accumulate their stats w.r.t. variants.
        // For tandem-contained variants, stats pertain to entire tandem.
        for_alignment_in_gaf(gaf_path, [&](GafAlignment& alignment){
            auto& path = alignment.get_path();
            if (path.size()<2) return;
            const bool path_is_reverse = path.front().second;
            const auto l = variant_graph.is_dangling_node(stoll(path.front().first));
            const auto r = variant_graph.is_dangling_node(stoll(path.back().first));
            const bool is_spanning = l and r;
            if (l) l_coverage++;
            if (r) r_coverage++;

            int32_t max_length;
            float max_identity;
            vector<interval_t> ref_intervals, singleton_interval, no_interval;
            variant_graph.for_each_vcf_record(path, [&](size_t id, const vector<edge_t>& edges_of_the_record, const VcfRecord& _record){
                auto record = _record;
                variant_graph.vcf_record_to_path_intervals(path,edges_of_the_record,variant_flank_length,ref_intervals);
                max_identity=0; max_length=0;
                for (const auto& interval: ref_intervals) {
                    AlignmentSummary summary;
                    singleton_interval.clear(); singleton_interval.push_back(interval);
                    for_cigar_interval_in_alignment(false, alignment, singleton_interval, no_interval, [&](const CigarInterval& i, const interval_t& interval){
                        int32_t cigar_length = i.length;
                        if (i.code==1) cigar_length=abs(i.query_stop-i.query_start);  // I operation
                        if (cigar_length>0) summary.update(i,true);
                    },{});

                    // If we are in a DEL, arbitrarily augment the identity with as many matches as the length of the DEL
                    if (record.sv_type==VcfReader::TYPE_DELETION) {
                        CigarInterval placeholder;
                        coord_t variant_coord;
                        record.get_reference_coordinates(false, variant_coord);
                        placeholder.length = variant_coord.second - variant_coord.first;
                        placeholder.code = 7;   // '=' operation
                        summary.update(placeholder,true);
                    }
                    max_identity=max(max_identity,summary.compute_identity());
                    max_length=max(max_length,interval.second-interval.first);
                }

                // Update the max observed identity for this VCF record
                max_observed_identities[id]=max(max_observed_identities[id],max_identity);

                // Update the region length for this VCF record to the max length over all intervals (arbitrary).
                if (variant_supports[id].length_of_evaluated_region==0 || is_spanning) variant_supports[id].length_of_evaluated_region=max(variant_supports[id].length_of_evaluated_region,max_length);

                // Update the histogram for this VCF record to the max identity over all intervals (arbitrary).
                variant_supports[id].identity[is_spanning][path_is_reverse].emplace_back(max_identity);
            });
        });
        const double total_coverage = float(l_coverage+r_coverage)/2.0;

        path output_path = subdir / "annotated.vcf";
        ofstream out_file(output_path);
        if (not out_file.is_open() or not out_file.good()) throw runtime_error("ERROR: file could not be written: " + output_path.string());
        out_file << "##INFO=<ID=" << label << ",Number=27,Type=Float,Description=\"Coverage computed by hapestry stratified by log2 Q values (6 bins), read is spanning (2 bins, boolean), and is reverse (2 bins, boolean), with bin edges q = 2,3,4,5,6,7 open ended both ends, log2 meaning q1 corresponding to identity 50% q2=75% etc. first value in vector is avg window coverage, last 2 values are: is_tandem (0/1) and length_of_evaluated_region (bp)\",Source=\"hapestry\",Version=\"0.0.0.0.0.1\">\n";
        out_file << "##INFO=<ID=" << label << "_MAX,Number=1,Type=Float,Description=\"The maximum observed identity for any alignment that spanned this variant (DELs are given pseudo-matches equal to their length)\",Source=\"hapestry\",Version=\"0.0.0.0.0.1\">\n";
        out_file << "##INFO=<ID=" << label << "_NEIGHBORS,Number=1,Type=Float,Description=\"The number of variants that shared the window/region (including this one)\",Source=\"hapestry\",Version=\"0.0.0.0.0.1\">\n";
        vcf_reader.print_minimal_header(out_file);
        string s;
        for (size_t v=0; v<variant_supports.size(); v++){
            VcfRecord& record = variant_graph.vcf_records[v];
            variant_supports.at(v).get_support_string(s);
            record.info.append(";");
            record.info.append(label);
            record.info.append("=");
            record.info.append(to_string(total_coverage));
            record.info.append(",");
            record.info.append(s);
            record.info.append(",");
            record.info.append(to_string(int(variant_supports.at(v).is_tandem)));
            record.info.append(",");
            record.info.append(to_string(variant_supports.at(v).length_of_evaluated_region));
            record.info.append(";");
            record.info.append(label);
            record.info.append("_MAX");
            record.info.append("=");
            record.info.append(to_string(max_observed_identities.at(v)));
            record.info.append(";");
            record.info.append(label);
            record.info.append("_NEIGHBORS");
            record.info.append("=");
            record.info.append(to_string(variant_supports.size()));
            record.print(out_file);
            out_file << '\n';
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
        const string& label,
        const path& vcf,
        size_t n_threads,
        int32_t flank_length,
        int32_t interval_max_length,
        int32_t min_sv_length,
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
                                     std::cref(label),
                                     std::cref(output_dir),
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
}


void annotate(
        path vcf,
        path output_dir,
        path windows_bed,                // Override the interval graph if this is provided
        path tandem_bed,
        path bam_csv,
        path ref_fasta,
        const string& label,
        int32_t flank_length,
        int32_t interval_max_length,
        int32_t min_sv_length,
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
                false,
                force_unique_reads,
                false,
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
                false,
                false,
                false,
                force_unique_reads,
                false,
                flank_length
        );
    }

    cerr << t << "Writing sequences to disk" << '\n';

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
        cerr << t << "Peak memory usage: " << get_peak_memory_usage() << '\n';
        sleep_for(seconds(60));
        throw runtime_error("DEBUG EARLY EXIT");
    }

    // Generate GFAs/GAFs/CSVs and folder structure for every VCF * every region
    // By default, all of these files will be stored in /dev/shm and then copied into the output dir as a final step.
    // TODO: create option to use /dev/shm/ as staging dir
    // Absolutely must delete the /dev/shm/ copy or warn the user at termination
    //
    cerr << "Generating graph alignments for VCF: " << vcf << '\n';

    compute_graph_evaluation(
            contig_interval_trees,
            contig_tandems,
            region_transmaps,
            ref_sequences,
            regions,
            label,
            vcf,
            n_threads,
            flank_length,
            interval_max_length,
            min_sv_length,
            staging_dir
    );

    auto vcf_prefix = get_vcf_name_prefix(vcf);
    path out_vcf = output_dir / ("annotated.vcf");
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
    for (size_t i=0; i<regions.size(); i++){
        const auto& region = regions[i];
        path sub_vcf = output_dir / region.to_unflanked_string('_', flank_length) / "annotated.vcf";

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

    cerr << t << "Peak memory usage: " << get_peak_memory_usage() << '\n';
    cerr << t << "Done" << '\n';
}


int main (int argc, char* argv[]){
    path output_dir;
    path windows_bed;
    path tandem_bed;
    string label = "HAPESTRY_SUPPORT";
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
            "--label",
            label,
            "Label to give this annotation in the INFO field of the VCF");

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

    app.add_flag("--debug", debug, "Invoke this to add more logging and output");

    app.add_flag("--force_unique_reads", force_unique_reads, "Invoke this to add append each read name with the sample name so that inter-sample read collisions cannot occur");

    app.add_flag("--bam_not_hardclipped", bam_not_hardclipped, "Invoke this if you expect your BAMs NOT to contain ANY hardclips. Saves time on iterating.");

    try{
        app.parse(argc, argv);
    }
    catch (const CLI::ParseError &e) {
        return app.exit(e);
    }

    annotate(
        vcf,
        output_dir,
        windows_bed,
        tandem_bed,
        bam_csv,
        ref_fasta,
        label,
        flank_length,
        interval_max_length,
        min_sv_length,
        n_threads,
        debug,
        force_unique_reads,
        bam_not_hardclipped
    );

    return 0;
}
