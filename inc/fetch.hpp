#pragma once

#include "TransitiveMap.hpp"
#include "IntervalGraph.hpp"
#include "Authenticator.hpp"
#include "VcfReader.hpp"
#include "windows.hpp"
#include "Region.hpp"
#include "Timer.hpp"
#include "CLI11.hpp"
#include "misc.hpp"
#include "bed.hpp"
#include "bam.hpp"

// #include "bindings/cpp/WFAligner.hpp"
// using namespace wfa;

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

namespace sv_merge{

using sample_region_read_map_t = unordered_map <string, unordered_map <Region, vector<StrandedQSequence> > >;
using sample_region_coord_map_t = unordered_map <string, unordered_map <Region, vector<pair<string, CigarInterval> > > >;
using sample_region_flanked_coord_map_t = unordered_map <string, unordered_map <Region, vector<pair<string, pair<CigarInterval,CigarInterval> > > > >;



class FetchConfig {
public:
    // store BAM tags as strings for each read, corresponding to their tag names provided here
    vector<string> tags_to_fetch;

    // how many threads to use (parallelized by BAM, not region)
    int64_t n_threads;

    // skip excessively long reads, greater than this threshold (useful for supplementaries, fragmented alignments)
    int32_t max_length;

    // length of sequence that is considered flanking, and to be tracked in Transmaps as additional data. The query coordinates corresponding to the inner flank bounds are stored.
    int32_t flank_length;

    // if non-zero, attempt to fetch this many clipped bases beyond the end of non-spanning sequences
    int32_t max_clip_fetch;

    // for any read to be fetched it must, among all its alignments, cover the left and right bounds
    bool require_spanning;

    // Retain information about where the query crosses the flank bounds, for any read to be fetched it must, among all its alignments, cover the left and right bounds AND inner flank bounds.
    bool get_flank_query_coords;

    // only consider the first read observed for each window in the BAM (for omitting misassemblies in asm-to-ref BAMs)
    bool first_only;

    // if true, attempt to force unique read names by using their sample name as a suffix
    bool append_sample_to_read;

    // if true, complement reverse sequences so they are given in ref forward orientation
    bool force_forward;

    // if true, re-write clipped subseqences coordinates as though they are in the original unclipped sequence coordinate space
    bool unclip_coords;

    // if true, also fetch qualities of equal length to the sequence fetched
    bool get_qualities;

    // do not throw error if tags are not found in the sequences
    bool allow_unused_tags;

    FetchConfig();
};


void fetch_reads(
        Timer& t,
        vector<Region>& regions,
        path bam_csv,
        const FetchConfig& config,
        unordered_map<Region,TransMap>& region_transmaps
        );


/**
 *
 * @param t for printing fetch time per BAM
 * @param regions subregions which will be extracted, reads will be clipped to fit the bounds, must be sorted and same contig. Regions are expected to be NON OVERLAPPING!
 * @param bam_csv path to headerless CSV containing format sample_name,bam_path
 * @param region_transmaps result object storing the output of this fetch
 */
void fetch_reads_from_clipped_bam(
        Timer& t,
        vector<Region>& regions,
        path bam_csv,
        const FetchConfig& config,
        unordered_map<Region,TransMap>& region_transmaps
);


void extract_subregion_coords_from_sample(
        Authenticator& authenticator,
        sample_region_coord_map_t& sample_to_region_coords,
        const string& sample_name,
        const vector<Region>& subregions,
        const FetchConfig& config,
        path bam_path
);


void extract_flanked_subregion_coords_from_sample(
        Authenticator& authenticator,
        sample_region_coord_map_t& sample_to_region_coords,
        const string& sample_name,
        const vector<Region>& subregions,
        const FetchConfig& config,
        path bam_path
);


void extract_subsequences_from_sample_thread_fn(
        Authenticator& authenticator,
        sample_region_read_map_t& sample_to_region_reads,
        const vector <pair <string,path> >& sample_bams,
        const vector<Region>& regions,
        const FetchConfig& config,
        atomic<size_t>& job_index
);


}
