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

using sample_region_read_map_t = unordered_map <string, unordered_map <Region, vector<Sequence> > >;
using sample_region_coord_map_t = unordered_map <string, unordered_map <Region, vector<pair<string, CigarInterval> > > >;

void fetch_reads(
        Timer& t,
        vector<Region>& regions,
        path bam_csv,
        int64_t n_threads,
        bool require_spanning,
        unordered_map<Region,TransMap>& region_transmaps
);


void fetch_reads_from_clipped_bam(
        Timer& t,
        vector<Region>& regions,
        path bam_csv,
        int64_t n_threads,
        int32_t interval_max_length,
        bool require_spanning,
        unordered_map<Region,TransMap>& region_transmaps,
        bool append_sample_to_read = false
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

}
