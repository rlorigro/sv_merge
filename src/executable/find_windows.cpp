#include "read_optimizer.hpp"
#include "TransitiveMap.hpp"
#include "IntervalGraph.hpp"
#include "VcfReader.hpp"
#include "Timer.hpp"
#include "CLI11.hpp"
#include "misc.hpp"
#include "bed.hpp"
#include "bam.hpp"

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
using std::min;
using std::cref;
using std::ref;


using namespace sv_merge;

using sample_region_read_map_t = unordered_map <string, unordered_map <Region, vector<Sequence> > >;


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


void construct_windows_from_vcf_and_bed(path tandem_bed, path vcf, int64_t flank_length, int64_t interval_max_length, vector<Region>& regions){
    VcfReader vcf_reader(vcf);
    vcf_reader.min_sv_length = 0;

    pair<uint64_t, uint64_t> coord;
    unordered_set<uint32_t> sample_ids;
    unordered_set<string> sample_names;

    unordered_map<string,vector<labeled_interval_t> > contig_intervals;

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

        contig_intervals[r.chrom].emplace_back(coord, sample_names);
    });

    for (auto& [contig, intervals]: contig_intervals){
        // Iterate the VCF file and construct a vector of labeled intervals for the IntervalGraph
        // Need to append `interval_padding` onto intervals and then subtract it afterwards (if desired)
        IntervalGraph<string> g(intervals);

        g.for_each_connected_component_interval([&](interval_t& interval, unordered_set<string>& values){
            if (interval.second - interval.first > interval_max_length){
                return;
            }

            regions.emplace_back(contig, interval.first, interval.second);
        });
    }
}


void load_windows_from_bed(path windows_bed, vector<Region>& regions){
    for_region_in_bed_file(windows_bed, [&](const Region &r) {
        regions.push_back(r);
    });

    // How to sort regions
    auto left_comparator = [](const Region& a, const Region& b){
        if (a.name != b.name) {
            return a.name < b.name;
        }
        else {
            return a.start < b.start;
        }
    };

    sort(regions.begin(), regions.end(), left_comparator);
}


void find_windows(
        path output_dir,
        path tandem_bed,
        path vcf,
        int64_t interval_max_length,
        int64_t flank_length,
        int64_t chunk_size,
        bool debug
        ){

    if (ghc::filesystem::exists(output_dir)){
        throw runtime_error("ERROR: output dir exists already: " + output_dir.string());
    }
    else{
        ghc::filesystem::create_directories(output_dir);
    }

    Timer t;

    cerr << t << "Initializing" << '\n';

    vector <pair <string,path> > bam_paths;
    vector<Region> regions;

    cerr << t << "Constructing windows" << '\n';

    construct_windows_from_vcf_and_bed(tandem_bed, vcf, flank_length, interval_max_length, regions);

    int64_t n_chunks = 0;
    int64_t r = 0;
    path output_path;
    ofstream file;

    if (chunk_size == 0){
        chunk_size = numeric_limits<int64_t>::max();
    }

    while (r < regions.size()) {
        if (r % chunk_size == 0){
            output_path = output_dir / ("windows_" + to_string(n_chunks) + ".bed");
            file.close();
            file.open(output_path);
            n_chunks++;
        }

        auto& region = regions[r];
        region.start -= flank_length;
        region.stop += flank_length;

        file << region.name << '\t' << region.start << '\t' << region.stop << '\n';

        r++;
    }

    cerr << t << "Done" << '\n';

}


int main (int argc, char* argv[]){
    path output_dir;
    path windows_bed;
    path tandem_bed;
    path bam_csv;
    path vcf;
    path ref;
    int64_t interval_max_length;
    int64_t flank_length;
    int64_t chunk_size = 0;
    bool debug = false;

    CLI::App app{"App description"};

    app.add_option(
            "--output_dir",
            output_dir,
            "Path to output directory which must not exist yet")
            ->required();

    app.add_option(
            "--chunk_size",
            chunk_size,
            "How many chunks to produce, for scattering of intervals. Setting to 0 indicates no chunking");

    app.add_option(
            "--vcf",
            vcf,
            "Path to VCF file containing variants to be merged")
            ->required();

    app.add_option(
            "--tandems",
            tandem_bed,
            "Path to BED file containing tandem repeat locations (for automated window generation)");

    app.add_option(
            "--interval_max_length",
            interval_max_length,
            "Maximum reference window size")
            ->required();

    app.add_option(
            "--flank_length",
            flank_length,
            "How much flanking sequence to use when fetching and aligning reads")
            ->required();

    app.add_flag("--debug", debug, "Invoke this to add more logging and output");

    CLI11_PARSE(app, argc, argv);

    find_windows(
        output_dir,
        tandem_bed,
        vcf,
        interval_max_length,
        flank_length,
        chunk_size,
        debug
    );

    return 0;
}
