#include "windows.hpp"
#include "VcfReader.hpp"
#include "bed.hpp"


namespace sv_merge{

void construct_windows_from_vcf_and_bed(path tandem_bed, path vcf, int32_t flank_length, int32_t interval_max_length, vector<Region>& regions){
    VcfReader vcf_reader(vcf);
    vcf_reader.min_qual = numeric_limits<float>::min();
    vcf_reader.min_sv_length = 0;
    vcf_reader.progress_n_lines = 100'000;

    unordered_set<uint32_t> sample_ids;
    unordered_set<string> sample_names;

    coord_t coord;
    interval_t interval;

    unordered_map<string,vector <pair <interval_t, bool> > > contig_intervals;

    cerr << "Reading BED... " << '\n';
    // Add every tandem region with a generic name as the label
    for_region_in_bed_file(tandem_bed, [&](const Region& r){
        interval.first = r.start;
        interval.second = r.stop;
        contig_intervals[r.name].emplace_back(interval, true);
    });

    cerr << "Reading VCF... " << '\n';
    // Add every VCF allele interval with the sample name as the label
    vcf_reader.for_record_in_vcf([&](VcfRecord& r){
        r.get_samples_with_alt(sample_ids);

        if (sample_ids.empty()){
            return;
        }

        r.get_reference_coordinates(false, coord);

        contig_intervals[r.chrom].emplace_back(coord, false);
    });

    cerr << "Computing intervals... " << '\n';

    // For each contig in reference, compute intervals
    for (auto& [contig, intervals]: contig_intervals){
        cerr << "\tStarting: " << contig << '\n';

        if (contig != "chr20"){
            continue;
        }

        vector<interval_t> components;

        get_overlapping_components(flank_length, intervals, components);

        for (const auto& c: components){
            if (c.second - c.first > interval_max_length){
                continue;
            }

            regions.emplace_back(contig, c.first, c.second);
        }
    }
}


void get_overlapping_components(int32_t min_gap_length, vector <pair <interval_t, bool> >& labeled_intervals, vector <interval_t>& result){
    // How to sort labeled intervals
    auto left_comparator = [](const pair <interval_t, bool>& a, const pair <interval_t, bool>& b){
        return a.first.first < b.first.first;
    };

    sort(labeled_intervals.begin(), labeled_intervals.end(), left_comparator);

    result.clear();
    result.emplace_back(labeled_intervals.front().first);

    bool has_non_tandem_interval = false;

    for (auto& [interval,is_tandem] : labeled_intervals){
        if (interval.first > interval.second){
            throw runtime_error("ERROR: interval start is greater than interval stop: " + to_string(interval.first) + ',' + to_string(interval.second));
        }

        interval.second += min_gap_length;

        // If the current interval has no overlap with any of the previous intervals, start new component
        if (interval.first > result.back().second){
            // If the only intervals in this component are tandem intervals, just overwrite it, because it is not needed
            if (not has_non_tandem_interval){
                result.back() = interval;
            }
            else{
                // Otherwise, start a new interval
                result.emplace_back(interval);
            }

            has_non_tandem_interval = not is_tandem;
        }

        // Maintain the maximum bound of this component
        if (interval.second > result.back().second){
            result.back().second = interval.second;
        }

        // Make sure there is at least one non-tandem interval
        if (not is_tandem){
            has_non_tandem_interval = true;
        }

        // Return the stored interval to original state
        interval.second -= min_gap_length;
    }

    for (auto& item: result){
        item.second -= min_gap_length;
    }

}



}
