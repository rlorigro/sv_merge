#include "TransitiveMap.hpp"
#include "CLI11.hpp"
#include "bed.hpp"
#include "bam.hpp"

#include <unordered_map>
#include <stdexcept>
#include <iostream>
#include <fstream>

using std::unordered_map;
using std::runtime_error;
using std::ifstream;
using std::cerr;


using namespace sv_merge;


void extract_reads_from_region(Region r, const unordered_map<string,path>& bam_paths, TransMap& transmap){
    // Do extraction
}


void load_bam_paths(path bam_csv, unordered_map<string,path>& bam_paths){
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
            bam_paths[sample_name] = bam_path;

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


void extract_read_subsequences_from_region(const Region& region, path bam_path, Tra){

}


void hapestry(
        path windows_bed,
        path tandem_bed,
        path bam_csv,
        path ref,
        int64_t interval_padding,
        int64_t interval_max_length,
        int64_t flank_length){

    // Iterate the BED file windows and construct a vector of labeled intervals for the IntervalGraph
    // Need to append `interval_padding` onto intervals and then subtract it afterwards
//    for_region_in_bed_file(windows_bed, [&](const Region &r) {
//
//    });

    // Select regions using IntervalGraph

    // Load BAM paths as a map with sample->bam
    unordered_map<string,path> bam_paths;
    load_bam_paths(bam_csv, bam_paths);

    for (auto& item: bam_paths){
        cerr << item.first << ' ' << item.second << '\n';
    }

    // For each region


}


int main (int argc, char* argv[]){
    path windows_bed;
    path tandem_bed;
    path ref;
    path bam_csv;
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
        ref,
        interval_padding,
        interval_max_length,
        flank_length
    );

    return 0;
}
