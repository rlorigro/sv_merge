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


class VariantSupport{
public:
    // is_spanning -> is_reverse -> Q
    array <array <vector<float>, 2>, 2> identity;
    int32_t length_of_evaluated_region = -1;
    bool is_tandem = false;

    [[nodiscard]] size_t get_coverage(bool is_spanning, bool is_reverse) const;
    void get_identity_distribution(bool is_spanning, bool is_reverse, array<size_t,6>& result) const;
    void get_support_string(string& s) const;
};


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
            q = max(0.0f,q-2.0f);

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

        write_time_log(subdir, vcf_name_prefix, time_csv, success);

        // Skip remaining steps for this region/tool if alignment failed and get the next job index for the thread
        if (not success) {
            cerr << "WARNING: Command timed out: " << command << '\n';
            i = job_index.fetch_add(1);
            continue;
        }

        size_t l_coverage = 0;
        size_t r_coverage = 0;
        int32_t variant_flank_length = 50;

        vector<VariantSupport> variant_supports(variant_graph.vcf_records.size());

        // Iterate the alignments and accumulate their stats w.r.t. variants
        // For tandem-contained variants, stats pertain to entire tandem
        unordered_map<string,size_t> counter;
        for_alignment_in_gaf(gaf_path, [&](GafAlignment& alignment){
            auto& path = alignment.get_path();
            if (path.size() < 2){
                return;
            }

            bool path_is_reverse = path.front().second;

            nid_t id_front = stoll(path.front().first);
            nid_t id_back = stoll(path.back().first);

            // Accumulate total flank coverage (to be averaged later)
            auto l = variant_graph.is_dangling_node(id_front);
            auto r = variant_graph.is_dangling_node(id_back);
            bool is_spanning = l and r;

            // Fetch (or construct) the count for this read name (default = 0)
            auto& c = counter[alignment.get_query_name()];

            // Create a unique name and add the path to variant graph
            auto name = alignment.get_query_name() + "_" + to_string(c);

            variant_graph.for_each_vcf_record(path, [&](size_t id, const vector<edge_t>& edges_of_the_record, const VcfRecord& _record){
//                cerr << "TESTING: " << _record.id << '\n';

                auto record = _record;

                if (record.sv_type == VcfReader::TYPE_BREAKEND){
                    throw runtime_error("ERROR: tandem-aware annotation not implemented for BNDs");
                }

                unordered_set<edge_t> vcf_edges(edges_of_the_record.begin(), edges_of_the_record.end());

                int32_t path_length = 0;
                vector<int32_t> lengths;
                vector<handle_t> handles;
                for (const auto& [n,is_reverse]: path){
                    auto h = variant_graph.graph.get_handle(stoll(n), is_reverse);
                    auto length = int32_t(variant_graph.graph.get_length(h));

                    lengths.emplace_back(length);
                    handles.emplace_back(h);
                    path_length += length;
                }

                coord_t ref_coord;
                record.get_reference_coordinates(false, ref_coord);
                bool is_tandem = false;

                contig_interval_trees.at(region.name).overlap_find_all({ref_coord.first, ref_coord.second}, [&](auto iter) {
                    is_tandem = true;
                    return false;
                });

                int32_t target_min = variant_graph.get_flank_boundary_left(region.name, ref_coord.first, variant_flank_length);
                int32_t target_max = variant_graph.get_flank_boundary_right(region.name, ref_coord.second, variant_flank_length);

                // Coordinates of path sequence to evaluate
                int32_t path_min = path_length + 1;
                int32_t path_max = 0;

                // Temp variable which is incremented during traversal
                int32_t path_coord = 0;
                Region node_region;

                // Observed path coordinates that match target_min/max (there could be more than one, if path is cyclic)
                vector<int32_t> path_target_mins;
                vector<int32_t> path_target_maxes;

//                cerr << id << ',' << record.id << ',' << alignment.get_query_name() << '\n';
                for (size_t p=0; p < path.size(); p++){
                    auto h1 = handles[p];
                    auto length = lengths[p];

                    variant_graph.get_region_of_node(h1, node_region);

//                    cerr << variant_graph.graph.get_id(h1) << ',' << variant_graph.graph.get_is_reverse(h1) << ',' << node_region.to_string() << '\n';
                    // If the node contains one of the target flanks, record its location
                    if (point_is_contained(target_min, node_region, true)){
//                        cerr << "FOUND new target_min: " << target_min << ',' << path_coord << '\n';
                        path_target_mins.emplace_back(path_coord + abs(node_region.start - target_min));
                    }
                    if (point_is_contained(target_max, node_region, true)){
//                        cerr << "FOUND new target_max: " << target_max << ',' << path_coord << '\n';
                        path_target_maxes.emplace_back(path_coord + abs(node_region.start - target_max));
                    }

                    if (p < path.size() - 1) {
                        auto h2 = handles[p+1];

                        // Get forward and reverse edges
                        edge_t e(h1,h2);
                        edge_t e2(variant_graph.graph.flip(h2), variant_graph.graph.flip(h1));

                        // Update the path ref coord so it is now up to the edge junction
                        path_coord += length;

                        // If the edge is one of the relevant VCF edges then track min and max
                        if (vcf_edges.find(e) != vcf_edges.end() or vcf_edges.find(e2) != vcf_edges.end()) {
                            if (path_coord < path_min) {
                                path_min = path_coord;
                            }
                            if (path_coord > path_max) {
                                path_max = path_coord;
                            }
                        }
                    }
                }

                // Find the smallest path interval that intersects the target region
                coord_t target_interval = {0,path_length};

                for (auto a: path_target_mins){
                    for (auto b: path_target_maxes){
//                        cerr << "Testing: " << a << '<' << path_min << '<' << path_max << '<' << b << '\n';
                        if (point_is_contained(path_min, {a, b}, true) and point_is_contained(path_max, {a, b}, true)){
                            auto x = b-a;
                            auto y = target_interval.second - target_interval.first;
                            if (x < y){
                                target_interval = {min(a,b), max(a,b)};
                            }
                        }
                    }
                }

                // Make sure we aren't exceeding 0 or path_length
                target_interval.first = max(0, target_interval.first);
                target_interval.second = min(path_length, target_interval.second);

                vector<interval_t> ref_intervals = {target_interval};
                vector<interval_t> query_intervals = {};

                AlignmentSummary summary;
                for_cigar_interval_in_alignment(false, alignment, ref_intervals, query_intervals, [&](const CigarInterval& i, const interval_t& interval){
                    summary.update(i,true);
                },{});

                variant_supports[id].identity[is_spanning][path_is_reverse].emplace_back(summary.compute_identity());
                variant_supports[id].is_tandem = is_tandem;
                variant_supports[id].length_of_evaluated_region = target_interval.second - target_interval.first;

//                cerr << region.to_string() << ',' << target_min << ',' << target_max << '\n';
//                cerr << id << ',' << record.id << ',' << alignment.get_query_name() << ',' << path_length << ',' << target_interval.first << ',' << target_interval.second << ',' << summary.compute_identity() << ',' << is_tandem << '\n';
            });

            // Increment coverage if this alignment covers the left flanking node
            if (l) {
                l_coverage++;
            }
            if (r) {
                r_coverage++;
            }

            c++;
        });

        auto total_coverage = float(l_coverage + r_coverage) / 2.0;

        path output_path = subdir / "annotated.vcf";
        ofstream out_file(output_path);
        if (not out_file.is_open() or not out_file.good()){
            throw runtime_error("ERROR: file could not be written: " + output_path.string());
        }

        out_file << "##INFO=<ID=" + label + ",Number=27,Type=Float,Description=\"Coverage computed by hapestry stratified by log2 Q values (6 bins), read is spanning (2 bins, boolean), and is reverse (2 bins, boolean), with bin edges q = 2,3,4,5,6,7 open ended both ends, log2 meaning q1 corresponding to identity 50% q2=75% etc. first value in vector is avg window coverage, last 2 values are: is_tandem (0/1) and length_of_evaluated_region (bp)\",Source=\"hapestry\",Version=\"0.0.0.0.0.1\">" << '\n';
        vcf_reader.print_minimal_header(out_file);
        string s;

        for (size_t v=0; v<variant_supports.size(); v++){
            auto& record = variant_graph.vcf_records[v];

            variant_supports[v].get_support_string(s);
            record.info += ";" + label + "=" + to_string(total_coverage) + ',' + s + ',' + to_string(int(variant_supports[v].is_tandem) + ',' + variant_supports[v].length_of_evaluated_region);
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
        construct_windows_from_vcf_and_bed(ref_sequences, contig_tandems, {vcf}, flank_length, interval_max_length, regions, bed_log_path);
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
                false,
                true,
                force_unique_reads,
                false
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
            staging_dir
    );

    auto vcf_prefix = get_vcf_name_prefix(vcf);
    path out_vcf = output_dir / (vcf_prefix + "_annotated.vcf");
    ofstream out_file(out_vcf);

    if (not (out_file.is_open() and out_file.good())){
        throw runtime_error("ERROR: could not write BED log file: " + out_vcf.string());
    }

    ifstream input_vcf(vcf);

    if (not (input_vcf.is_open() and input_vcf.good())){
        throw runtime_error("ERROR: could not write BED log file: " + vcf.string());
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
    // TODO: fix duplicated version line
    for (size_t i=0; i<regions.size(); i++){
        const auto& region = regions[i];
        path sub_vcf = output_dir / region.to_unflanked_string('_', flank_length) / "annotated.vcf";

        ifstream file(sub_vcf);
        while (getline(file, line)){
            if (line.starts_with('#')){
                if (i == 0){
                    out_file << line << '\n';
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
            "How much flanking sequence to use when fetching and aligning reads")
            ->required();

    app.add_flag("--debug", debug, "Invoke this to add more logging and output");

    app.add_flag("--force_unique_reads", force_unique_reads, "Invoke this to add append each read name with the sample name so that inter-sample read collisions cannot occur");

    app.add_flag("--bam_not_hardclipped", bam_not_hardclipped, "Invoke this if you expect your BAMs NOT to contain ANY hardclips. Saves time on iterating.");

    app.parse(argc, argv);

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
        n_threads,
        debug,
        force_unique_reads,
        bam_not_hardclipped
    );

    return 0;
}
