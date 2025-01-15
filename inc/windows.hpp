#pragma once

#include "IntervalGraph.hpp"
#include "VcfReader.hpp"
#include "bed.hpp"


namespace sv_merge{

void construct_windows_from_vcf_and_bed(
        const unordered_map<string,string>& ref_sequences,
        const unordered_map<string,vector<interval_t> >& contig_tandems,
        const vector<path>& vcfs,
        int32_t flank_length,
        int32_t interval_max_length,
        int32_t min_sv_length,
        vector<Region>& regions,
        const path& bed_log_path="",
        bool use_confidence_intervals=false
        );

void construct_windows_from_vcf_and_bed(
        const unordered_map<string,vector<interval_t> >& contig_tandems,
        const vector<path>& vcfs,
        int32_t flank_length,
        int32_t interval_max_length,
        int32_t min_sv_length,
        vector<Region>& regions,
        bool use_confidence_intervals=false
        );

void construct_windows_from_vcf_and_bed(
        const unordered_map<string,vector<interval_t> >& contig_tandems,
        path vcf,
        int32_t flank_length,
        int32_t interval_max_length,
        int32_t min_sv_length,
        vector<Region>& regions,
        bool use_confidence_intervals=false
        );

void get_overlapping_components(int32_t min_gap_length, vector <pair <interval_t, bool> >& labeled_intervals, vector <interval_t>& result);

}
