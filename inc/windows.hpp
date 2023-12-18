#pragma once

#include "IntervalGraph.hpp"
#include "VcfReader.hpp"
#include "bed.hpp"


namespace sv_merge{

void construct_windows_from_vcf_and_bed(path tandem_bed, path vcf, int32_t flank_length, int32_t interval_max_length, vector<Region>& regions);

void get_overlapping_components(vector <pair <interval_t, bool> >& labeled_intervals, vector <interval_t>& result);

}
