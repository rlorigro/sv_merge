#include "Bam.hpp"
#include "Filesystem.hpp"

using ghc::filesystem::path;
using sv_merge::Sequence;
using sv_merge::Region;
using sv_merge::for_read_in_bam_region;

#include <iostream>
#include <stdexcept>

using std::runtime_error;
using std::string;
using std::cerr;



int main(){
    string region_string = "a:0-5000";

    Region r(region_string);

    cerr << r.name << ' ' << r.start << ' ' << r.stop << '\n';

    path bam_path_local = "/home/ryan/code/sv_merge/data/test_alignment_softclip_only_sorted.bam";

    for_read_in_bam_region(bam_path_local, region_string, [&](Sequence& s){
        cerr << s.name << ' ' << s.sequence << '\n';
    });

    return 0;
}

