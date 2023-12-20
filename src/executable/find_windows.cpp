#include "windows.hpp"
#include "Timer.hpp"
#include "CLI11.hpp"
#include "misc.hpp"
#include "bed.hpp"
#include "bam.hpp"

#include <unordered_map>
#include <exception>
#include <stdexcept>
#include <iostream>
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
using std::max;
using std::cref;
using std::ref;


using namespace sv_merge;


void find_windows(
        path output_dir,
        path tandem_bed,
        path vcf,
        int32_t interval_max_length,
        int32_t flank_length,
        size_t n_chunks,
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

    size_t r = 0;
    path output_path;
    ofstream file;

    size_t chunk_size = max(1ul,regions.size() / n_chunks);
    size_t n = 0;

    string prev_contig;

    while (r < regions.size()) {
        auto& region = regions[r];

//        cerr << region.start << ',' << region.stop << '\n';

        // If the number of chunks iterated exceeds the desired chunk size, start a new chunk.
        // Also, if a new contig is reached, start a new chunk.
        if (r % chunk_size == 0 or (prev_contig != region.name and not prev_contig.empty())){
            output_path = output_dir / ("windows_" + to_string(n) + ".bed");
            file.close();
            file.open(output_path);
            n++;
        }

        region.start -= flank_length;
        region.stop += flank_length;

        file << region.name << '\t' << region.start << '\t' << region.stop << '\n';

        prev_contig = region.name;
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
    size_t n_chunks = 1;
    bool debug = false;

    CLI::App app{"App description"};

    app.add_option(
            "--output_dir",
            output_dir,
            "Path to output directory which must not exist yet")
            ->required();

    app.add_option(
            "--n_chunks",
            n_chunks,
            "How many chunks to produce, for scattering of intervals.");

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
        n_chunks,
        debug
    );

    return 0;
}
