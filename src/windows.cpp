#include "windows.hpp"
#include "VcfReader.hpp"
#include "bed.hpp"


namespace sv_merge{


void construct_windows_from_vcf_and_bed(
        const unordered_map<string,vector<interval_t> >& contig_tandems,
        const vector<path>& vcfs,
        int32_t flank_length,
        int32_t interval_max_length,
        int32_t min_sv_length,
        vector<Region>& regions,
        bool use_confidence_intervals
        ){

    // Provide empty object
    unordered_map<string,string> ref_sequences;
    construct_windows_from_vcf_and_bed(
        ref_sequences,
        contig_tandems,
        vcfs,
        flank_length,
        interval_max_length,
        min_sv_length,
        regions,
        "",
        use_confidence_intervals);
}


void construct_windows_from_vcf_and_bed(
        const unordered_map<string,string>& ref_sequences,
        const unordered_map<string,vector<interval_t> >& contig_tandems,
        const vector<path>& vcfs,
        int32_t flank_length,
        int32_t interval_max_length,
        int32_t min_sv_length,
        vector<Region>& regions,
        const path& bed_log_path,
        bool use_confidence_intervals
        ){


    ofstream log_file;

    if (not bed_log_path.empty()){
        log_file.open(bed_log_path);

        if (not (log_file.is_open() and log_file.good())){
            throw runtime_error("ERROR: could not write BED log file: " + bed_log_path.string());
        }
    }

    interval_t interval;

    unordered_map<string,vector <pair <interval_t, bool> > > contig_intervals;

    for (const auto& [name, tandems]: contig_tandems){
        for (auto& t: tandems){
            if (t.first == t.second){
                throw runtime_error("ERROR: tandem BED interval start == stop: " + name + ' ' + to_string(interval.first) + ',' + to_string(interval.second));
            }

            contig_intervals[name].emplace_back(t, true);
        }
    }

    vector <pair <string,size_t> > vcf_omissions;

    for (const auto& vcf: vcfs) {
        cerr << "Reading VCF: " << vcf << '\n';

        vcf_omissions.emplace_back(vcf.filename(), 0);

        VcfReader vcf_reader(vcf);
        vcf_reader.min_qual = numeric_limits<float>::min();
        vcf_reader.min_sv_length = min_sv_length;
        vcf_reader.progress_n_lines = 100'000;

        unordered_set<int32_t> sample_ids;
        unordered_set<string> sample_names;

        coord_t coord;

        // Add every VCF allele interval with the sample name as the label
        vcf_reader.for_record_in_vcf([&](VcfRecord &r) {
            r.get_samples_with_alt(sample_ids);

            r.get_reference_coordinates(use_confidence_intervals, coord);

            // Skip large events in the population.
            // A better solution consists in pre-processing the VCF to transform every large event into a set of BNDs:
            // see `clean_bnds.cpp`.
            if (coord.second - coord.first > interval_max_length){
                vcf_omissions.back().second++;
                return;
            }

            contig_intervals[r.chrom].emplace_back(coord, false);
        });
    }

    cerr << "Computing intervals... " << '\n';

    size_t total_bp_omitted = 0;
    size_t total_windows_omitted = 0;
    // For each contig in reference, compute intervals
    for (auto& [contig, intervals]: contig_intervals){
        cerr << "\tStarting: " << contig << '\n';

        vector<interval_t> components;

        get_overlapping_components(flank_length, intervals, components);

        for (const auto& c: components){
            if (c.second - c.first > interval_max_length){
                if (not bed_log_path.empty()) {
                    log_file << contig << '\t' << c.first << '\t' << c.second << '\n';
                }

                total_bp_omitted += c.second - c.first;
                total_windows_omitted++;
                continue;
            }

            // Temporary coord to test bounds of flanks
            auto c_flanked = c;
            c_flanked.first -= flank_length;
            c_flanked.second += flank_length;

            // If ref_sequences are provided, check for consistency with contig lengths
            if (not ref_sequences.empty()){
                auto result = ref_sequences.find(contig);

                if (result == ref_sequences.end()){
                    throw runtime_error("ERROR: VCF region name not found in reference sequences: " + contig);
                }

                auto contig_length = int32_t(result->second.size());

                if (c_flanked.first < 0 or c_flanked.first >= contig_length or c_flanked.second < 0 or c_flanked.second > contig_length){
                    if (not bed_log_path.empty()) {
                        log_file << contig << '\t' << c.first << '\t' << c.second << '\n';
                    }

                    total_bp_omitted += c.second - c.first;
                    total_windows_omitted++;

                    cerr << "WARNING: skipping region for which flanking sequence would exceed bounds: " << contig << ':' << c.first << ',' << c.second << '\n';
                    continue;
                }
            }
            else{
                // Just do one check to fix any negative coords
                if (c_flanked.first < 0 or c_flanked.second < 0){
                    if (not bed_log_path.empty()) {
                        log_file << contig << '\t' << c.first << '\t' << c.second << '\n';
                    }

                    total_bp_omitted += c.second - c.first;
                    total_windows_omitted++;

                    cerr << "WARNING: skipping region for which flanking sequence would be < 0 (NO REF PROVIDED): " << contig << ':' << c.first << ',' << c.second << '\n';
                    continue;
                }
            }

            regions.emplace_back(contig, c.first, c.second);
        }
    }

    cerr << "Input variants skipped because longer than " << interval_max_length << "bp:" << '\n';
    for (const auto& [name,count]: vcf_omissions){
        cerr << name << ": " << count << '\n';
    }

    cerr << "Connected components skipped because too long or would exceed contig end when flanked:" << '\n';
    cerr << "Total windows: " << total_windows_omitted << '\n';
    cerr << "Total bp: " << total_bp_omitted << '\n';

    log_file.close();
}


void construct_windows_from_vcf_and_bed(
        const unordered_map<string,vector<interval_t> >& contig_tandems,
        path vcf,
        int32_t flank_length,
        int32_t interval_max_length,
        int32_t min_sv_length,
        vector<Region>& regions,
        bool use_confidence_intervals
        ){

    vector<path> vcfs = {vcf};
    construct_windows_from_vcf_and_bed(
            contig_tandems,
            vcfs,
            flank_length,
            interval_max_length,
            min_sv_length,
            regions,
            use_confidence_intervals
    );
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

    if (not has_non_tandem_interval){
        result.pop_back();
    }

    for (auto& item: result){
        item.second -= min_gap_length;
    }

}



}
