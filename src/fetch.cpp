#include "fetch.hpp"


namespace sv_merge{

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
    placeholder.query_start = numeric_limits<int32_t>::max();
    placeholder.query_stop = numeric_limits<int32_t>::min();
    placeholder.ref_start = numeric_limits<int32_t>::max();
    placeholder.ref_stop = numeric_limits<int32_t>::min();

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


}
