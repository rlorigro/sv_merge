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

namespace sv_merge{

using sample_region_read_map_t = unordered_map <string, unordered_map <Region, vector<StrandedQSequence> > >;
using sample_region_coord_map_t = unordered_map <string, unordered_map <Region, vector<pair<string, CigarInterval> > > >;
using sample_region_flanked_coord_map_t = unordered_map <string, unordered_map <Region, vector<pair<string, pair<CigarInterval,CigarInterval> > > > >;


void fetch_reads(
        Timer& t,
        vector<Region>& regions,
        path bam_csv,
        int64_t n_threads,
        unordered_map<Region,TransMap>& region_transmaps,
        bool require_spanning,
        bool append_sample_to_read = false,
        bool force_forward = false,
        bool get_qualities = false
);


/**
 *
 * @param t for printing fetch time per BAM
 * @param regions subregions which will be extracted, reads will be clipped to fit the bounds, must be sorted and same contig
 * @param bam_csv
 * @param n_threads
 * @param max_length skip excessively long reads, greater than this threshold (useful for supplementaries, fragmented alignments)
 * @param flank_length length of sequence that is considered flanking, and to be tracked in Transmaps as additional
 * data. The query coordinates corresponding to the inner flank bounds are stored.
 * @param region_transmaps
 * @param require_spanning for any read to be fetched it must, among all its alignments, cover the left and right bounds
 * @param first_only only consider the first read observed for each window in the BAM (for omitting misassemblies in asm-to-ref BAMs)
 * @param append_sample_to_read if true, attempt to force unique read names by using their sample name as a suffix
 * @param force_forward if true, complement reverse sequences so they are given in ref forward orientation
 */
void fetch_reads_from_clipped_bam(
        Timer& t,
        vector<Region>& regions,
        path bam_csv,
        int64_t n_threads,
        int32_t max_length,
        int32_t flank_length,
        unordered_map<Region,TransMap>& region_transmaps,
        bool require_spanning,
        bool get_flank_query_coords = false,
        bool first_only = false,
        bool append_sample_to_read = false,
        bool force_forward = false
);


void extract_subregion_coords_from_sample(
        GoogleAuthenticator& authenticator,
        sample_region_coord_map_t& sample_to_region_coords,
        const string& sample_name,
        const vector<Region>& subregions,
        bool require_spanning,
        bool unclip_coords,
        path bam_path
);


void extract_flanked_subregion_coords_from_sample(
        GoogleAuthenticator& authenticator,
        sample_region_coord_map_t& sample_to_region_coords,
        const string& sample_name,
        const vector<Region>& subregions,
        bool require_spanning,
        bool unclip_coords,
        path bam_path
);


void get_read_coords_for_each_subregion_in_bam(
        Timer& t,
        vector<Region>& regions,
        GoogleAuthenticator& authenticator,
        sample_region_flanked_coord_map_t& sample_to_region_coords,
        path bam,
        int64_t n_threads,
        int32_t flank_length,
        bool require_spanning,
        bool get_flank_query_coords
);


void extract_subsequences_from_sample_thread_fn(
        GoogleAuthenticator& authenticator,
        sample_region_read_map_t& sample_to_region_reads,
        const vector <pair <string,path> >& sample_bams,
        const vector<Region>& regions,
        bool require_spanning,
        bool force_forward,
        bool get_qualities,
        atomic<size_t>& job_index,
        const vector<string>& tags_to_fetch = {},
        bool allow_unused_tags = false
);


}
