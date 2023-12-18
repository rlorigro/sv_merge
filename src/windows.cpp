#pragma once

#include "IntervalGraph.hpp"
#include "VcfReader.hpp"
#include "bed.hpp"


namespace sv_merge{

void construct_windows_from_vcf_and_bed(path tandem_bed, path vcf, int32_t flank_length, int32_t interval_max_length, vector<Region>& regions){
    VcfReader vcf_reader(vcf);
    vcf_reader.min_sv_length = 0;
    vcf_reader.progress_n_lines = 100'000;
    vcf_reader.min_sv_length = 20;

    coord_t coord;
    labeled_interval_t interval;
    unordered_set<uint32_t> sample_ids;
    unordered_set<string> sample_names;

    unordered_map<string,vector<labeled_interval_t> > contig_intervals;

    cerr << "Reading BED... " << '\n';
    // Add every tandem region with a generic name as the label
    interval.second.emplace("_tandem_");
    for_region_in_bed_file(tandem_bed, [&](const Region& r){
        interval.first.first = r.start;
        interval.first.second = r.stop + flank_length;
        contig_intervals[r.name].emplace_back(interval);
    });

    cerr << "Reading VCF... " << '\n';
    // Add every VCF allele interval with the sample name as the label
    vcf_reader.for_record_in_vcf([&](VcfRecord& r){
        r.get_samples_with_alt(sample_ids);

        if (sample_ids.empty()){
            return;
        }

        r.get_reference_coordinates(true, coord);

        sample_names.clear();
        for (auto id: sample_ids){
            sample_names.emplace(vcf_reader.sample_ids[id]);
        }

        coord.second += flank_length;

        contig_intervals[r.chrom].emplace_back(coord, sample_names);
    });

    cerr << "Computing intervals... " << '\n';

    // For each contig in reference, compute intervals
    for (auto& [contig, intervals]: contig_intervals){
        cerr << "\tStarting: " << contig << '\n';

        // Iterate the VCF file and construct a vector of labeled intervals for the IntervalGraph
        // Need to append `interval_padding` onto intervals and then subtract it afterwards (if desired)
        IntervalGraph<string> g(intervals);

        g.for_each_connected_component_interval([&](interval_t& interval, unordered_set<string>& values){
            if (interval.second - interval.first > interval_max_length){
                return;
            }

            regions.emplace_back(contig, interval.first, interval.second - flank_length);
        });
    }
}

}
