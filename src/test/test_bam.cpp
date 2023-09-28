#include "Bam.hpp"
#include "misc.hpp"
#include "Filesystem.hpp"

using ghc::filesystem::path;
using sv_merge::Sequence;
using sv_merge::Region;
using sv_merge::for_read_in_bam_region;
using sv_merge::run_command;
using sv_merge::files_equal;

#include <stdexcept>
#include <iostream>
#include <fstream>

using std::runtime_error;
using std::ofstream;
using std::string;
using std::cerr;



int main(){
    string region_string = "a:0-5000";

    Region r(region_string);

    cerr << r.name << ' ' << r.start << ' ' << r.stop << '\n';

    path project_directory = path(__FILE__).parent_path().parent_path().parent_path();
    path data_directory = project_directory / "data";

    cerr << data_directory << '\n';

    path bam_path = data_directory / "test_alignment_softclip_only_sorted.bam";
    path output_path = data_directory / "test_bam_region_extract.fasta";

    ofstream file(output_path);

    for_read_in_bam_region(bam_path, region_string, [&](Sequence& s){
        file << '>' << s.name << '\n';
        file << s.sequence << '\n';
    });

    path region_output_bam_path = data_directory / "region.sam";

    string command;
    string result;

    command = "samtools view -h " + bam_path.string() + " " + region_string;
    run_command(command, region_output_bam_path);

    path region_output_fasta_path = data_directory / "region.fasta";

    command = "samtools fasta -0 " + region_output_fasta_path.string() + " " + region_output_bam_path.string();
    run_command(command);

    bool success = files_equal(output_path, region_output_fasta_path);

    cerr << success << '\n';

    return 0;
}

