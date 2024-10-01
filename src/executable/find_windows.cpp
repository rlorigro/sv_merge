#include "windows.hpp"
#include "Timer.hpp"
#include "CLI11.hpp"
#include "fasta.hpp"
#include "bed.hpp"

#include <unordered_map>
#include <exception>
#include <stdexcept>
#include <iostream>
#include <thread>
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
        path windows_bed,
        path ref_fasta,
        int32_t interval_max_length,
        int32_t min_sv_length,
        int32_t flank_length,
        size_t n_chunks
        ){

    output_dir = std::filesystem::weakly_canonical(output_dir);
    tandem_bed = std::filesystem::weakly_canonical(tandem_bed);
    windows_bed = std::filesystem::weakly_canonical(windows_bed);
    ref_fasta = std::filesystem::weakly_canonical(ref_fasta);
    vcf = std::filesystem::weakly_canonical(vcf);

    if (std::filesystem::exists(output_dir)){
        throw runtime_error("ERROR: output dir exists already: " + output_dir.string());
    }
    else{
        std::filesystem::create_directories(output_dir);
    }

    Timer t;

    cerr << t << "Initializing" << '\n';

    vector <pair <string,path> > bam_paths;
    vector<Region> regions;

    cerr << t << "Constructing windows" << '\n';

    unordered_map<string,vector<interval_t> > contig_tandems;
    interval_t interval;
    for_region_in_bed_file(tandem_bed, [&](const Region& r){
        interval.first = r.start;
        interval.second = r.stop;
        contig_tandems[r.name].emplace_back(interval);
    });

    vector<path> vcfs = {vcf};
    path bed_log_path = output_dir / "log.bed";

    unordered_map<string,string> ref_sequences;

    if (not ref_sequences.empty()) {
        cerr << t << "Loading ref sequences from FASTA" << '\n';

        // Load all chromosome sequences (in case of BND)
        for_sequence_in_fasta_file(ref_fasta, [&](const Sequence& s){
            ref_sequences[s.name] = s.sequence;
        });
    }

    if (windows_bed.empty()){
        cerr << t << "Constructing windows from VCFs and tandem BED" << '\n';

        construct_windows_from_vcf_and_bed(
                ref_sequences,
                contig_tandems,
                vcfs,
                flank_length,
                interval_max_length,
                min_sv_length,
                regions,
                bed_log_path,
                false);
    }
    else {
        cerr << t << "Reading BED file" << '\n';
        load_windows_from_bed(windows_bed, regions);
    }

    size_t r = 0;
    path output_path;
    path output_path_unflanked;
    ofstream file;
    ofstream file2;

    size_t chunk_size = max(1ul,regions.size() / n_chunks);
    size_t n = 0;

    while (r < regions.size()) {
        auto& region = regions[r];

//        cerr << region.start << ',' << region.stop << '\n';

        // If the number of chunks iterated exceeds the desired chunk size, start a new chunk.
        // Also, if a new contig is reached, start a new chunk.
        if (r % chunk_size == 0){
            output_path_unflanked = output_dir / ("windows_" + to_string(n) + "_unflanked.bed");
            file2.close();
            file2.open(output_path_unflanked);

            output_path = output_dir / ("windows_" + to_string(n) + ".bed");
            file.close();
            file.open(output_path);

            n++;
        }

        file2 << region.name << '\t' << region.start << '\t' << region.stop << '\n';

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
    path ref_fasta;
    path tandem_bed;
    path bam_csv;
    path vcf;
    path ref;
    int32_t interval_max_length = 50000;
    int32_t min_sv_length = 20;
    int32_t flank_length = 200;
    size_t n_chunks = 1;

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
            "--windows",
            windows_bed,
            "Path to BED file containing windows, which will override the typical behavior of "
            "find_windows, only keeping the reference bounds clipping and the chunking behavior");

    app.add_option(
            "--ref_fasta",
            ref_fasta,
            "Path to FASTA file containing reference sequence, which will be loaded in order to ensure"
            " the flanks dont exceed the bounds of the ref sequences");

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
            "--min_sv_length",
            min_sv_length,
            "Skip all variants less than this length (bp)")
            ->required();

    app.add_option(
            "--flank_length",
            flank_length,
            "How much flanking sequence to use when fetching and aligning reads")
            ->required();

    CLI11_PARSE(app, argc, argv);

    find_windows(
        output_dir,
        tandem_bed,
        vcf,
        windows_bed,
        ref_fasta,
        interval_max_length,
        min_sv_length,
        flank_length,
        n_chunks
    );

    return 0;
}
