#include "bed.hpp"
#include "bam.hpp"
#include "CLI11.hpp"

#include <iostream>
#include <stdexcept>

using std::runtime_error;
using std::cerr;


using namespace sv_merge;


void hapestry(
        path windows_bed,
        path tandem_bed,
        path bam_csv,
        path ref,
        int64_t interval_padding,
        int64_t interval_max_length,
        int64_t flank_length){

    for_region_in_bed_file(windows_bed, [&](const Region &r) {

    });

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



    return 0;
}
