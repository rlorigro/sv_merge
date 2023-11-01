#include "TransitiveMap.hpp"
#include "IntervalGraph.hpp"
#include "Authenticator.hpp"
#include "CLI11.hpp"
#include "misc.hpp"
#include "bed.hpp"
#include "bam.hpp"

#include <unordered_map>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <limits>

using std::numeric_limits;
using std::unordered_map;
using std::runtime_error;
using std::ifstream;
using std::cerr;


using namespace sv_merge;


void extract_reads_from_region(Region r, const unordered_map<string,path>& bam_paths, TransMap& transmap){
    // Do extraction
}


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


void construct_windows_from_vcf_and_bed(path tandem_bed, path vcf, vector<Region>& regions){
    vector<labeled_interval_t> intervals;

    // Iterate the VCF file and construct a vector of labeled intervals for the IntervalGraph
    // Need to append `interval_padding` onto intervals and then subtract it afterwards
    // TODO: fill in when vcf reader exists
    IntervalGraph<string> g(intervals);

    g.for_each_connected_component_interval([&](interval_t& interval, unordered_set<string>& values){
        cerr << interval.first << "," << interval.second;
        for (const auto& v: values){
            cerr << ' ' << v;
        }
        cerr << '\n';
    });
}


// This is used when only one half of the ref/query dual iterator `for_alignment_in_bam_region` is desired
void null_fn(const CigarInterval &intersection, const interval_t &interval){}


void extract_read_subsequences_from_region(GoogleAuthenticator& authenticator, const string& sample_name, const Region& region, path bam_path, TransMap& transmap){
    CigarInterval placeholder;
    placeholder.query_start = numeric_limits<int64_t>::max();
    placeholder.query_stop = numeric_limits<int64_t>::min();
    placeholder.ref_start = numeric_limits<int64_t>::max();
    placeholder.ref_stop = numeric_limits<int64_t>::min();

    // Keep track of the min and max observed query coordinates that intersect the region of interest
    unordered_map<string,CigarInterval> query_coords;
    unordered_map<string,string> query_sequences;

    // Keep the F (query) oriented sequence of all the reads if they have an alignment that intersects the region
    unordered_map<string,string> alignments;

    // The region of interest is defined in reference coordinate space
    vector<interval_t> ref_intervals = {{region.start, region.stop}};

    // Unused
    vector<interval_t> query_intervals;

    // Make sure the system has the necessary authentication env variables to fetch a remote GS URI
    authenticator.update();

    // Iterate each alignment in the ref region
    for_alignment_in_bam_region(bam_path, region.to_string(), [&](Alignment& alignment) {
        if (alignment.is_unmapped() or not alignment.is_primary()){
            return;
        }

        string name;
        alignment.get_query_name(name);

        cerr << name << '\n';

        // Check if this read/query has an existing coord, from a previously iterated supplementary alignment
        auto result = query_coords.find(name);

        // If no previous alignment, initialize with the max placeholder
        if (result == query_coords.end()){
            // Insert a new placeholder cigar interval and keep a reference to the value inserted
            result = query_coords.emplace(name, placeholder).first;

            // Insert a new empty sequence and keep a reference to the value inserted
            auto result3 = query_sequences.emplace(name,"");
            auto& x = result3.first->second;

            // Fill the value with the sequence
            // TODO: find a way to not extract the whole sequence?
            alignment.get_query_sequence(x);
        }

        auto& coord = result->second;

        for_cigar_interval_in_alignment(alignment, ref_intervals, query_intervals,
            [&](const CigarInterval& intersection, const interval_t& interval) {
                // If the alignment touches the START of the ref region, record the query position
                if (intersection.ref_start == region.start and intersection.query_start < coord.query_start){
                    coord.query_start = intersection.query_start;
                }

                // If the alignment touches the END of the region, record the query position
                if (intersection.ref_stop == region.stop and intersection.query_stop > coord.query_stop){
                    coord.query_stop = intersection.query_stop;
                }
            },
            null_fn
        );
    });

    for (const auto& [name, coords]: query_coords){
        if (coords.query_start != placeholder.query_start and coords.query_stop == placeholder.query_stop){
            auto i = coords.query_start;
            auto l = coords.query_stop - coords.query_start;
            transmap.add_read(name, query_sequences[name].substr(i,l));
            transmap.add_edge(sample_name, name);

            cerr << name << '\n';
            cerr << transmap.get_sequence(name) << '\n';
        }
    }
}


//void call_region(const Region& region, path bam_path, TransMap& transmap){
//    extract_read_subsequences_from_region(region, bam_path, transmap);
//
//}


void hapestry(
        path windows_bed,               // Override the interval graph if this is provided
        path tandem_bed,
        path bam_csv,
        path vcf,
        path ref,
        int64_t interval_padding,
        int64_t interval_max_length,
        int64_t flank_length){

    unordered_map<string,path> bam_paths;
    GoogleAuthenticator authenticator;
    vector<Region> regions;
    TransMap transmap;

    if (windows_bed.empty()){
        construct_windows_from_vcf_and_bed(tandem_bed, vcf, regions);
    }
    else{
        for_region_in_bed_file(tandem_bed, [&](const Region &r) {
            cerr << r.name << ' ' << r.start << ' ' << r.stop << '\n';
            regions.push_back(r);
        });
    }

    // Load BAM paths as a map with sample->bam
    for_each_sample_bam_path(bam_csv, [&](const string& sample_name, const path& bam_path){
        transmap.add_sample(sample_name);
        bam_paths[sample_name] = bam_path;
    });

    // For each region
    for (const auto& region: regions){
        for (const auto& [sample_name, bam]: bam_paths) {
            extract_read_subsequences_from_region(authenticator, sample_name, region, bam, transmap);
        }
    }

}


int main (int argc, char* argv[]){
    path windows_bed;
    path tandem_bed;
    path bam_csv;
    path vcf = ""; // TODO: add arg for this when vcf reader exists
    path ref;
    int64_t interval_padding;
    int64_t interval_max_length;
    int64_t flank_length;

    CLI::App app{"App description"};

    app.add_option(
            "--tandems",
            windows_bed,
            "Path to GFA containing phased non-overlapping segments")
            ->required();

    app.add_option(
            "--windows",
            tandem_bed,
            "Path to GFA containing phased non-overlapping segments")
            ->required();

    app.add_option(
            "--bam_csv",
            bam_csv,
            "Path to GFA containing phased non-overlapping segments")
            ->required();

    app.add_option(
            "--interval_padding",
            interval_padding,
            "Path to GFA containing phased non-overlapping segments")
            ->required();

    app.add_option(
            "--interval_max_length",
            interval_max_length,
            "Path to GFA containing phased non-overlapping segments")
            ->required();

    app.add_option(
            "--flank_length",
            flank_length,
            "Path to GFA containing phased non-overlapping segments")
            ->required();

    CLI11_PARSE(app, argc, argv);

    hapestry(
        windows_bed,
        tandem_bed,
        bam_csv,
        vcf,
        ref,
        interval_padding,
        interval_max_length,
        flank_length
    );

    return 0;
}
