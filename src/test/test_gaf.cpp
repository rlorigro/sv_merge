#include "Filesystem.hpp"
#include "gaf.hpp"

using ghc::filesystem::path;
using sv_merge::for_alignment_in_gaf;
using sv_merge::GafAlignment;

#include <iostream>
#include <stdexcept>

using std::runtime_error;
using std::cerr;


int main(){
    path project_directory = path(__FILE__).parent_path().parent_path().parent_path();
    path data_directory = project_directory / "data";

    path gaf_path = data_directory / "test_alignments.gaf";

    for_alignment_in_gaf(gaf_path, [&](GafAlignment& alignment){
        cerr <<
            alignment.get_path_string() << '\t' <<
            alignment.get_query_name() << '\t' <<
            alignment.get_ref_name() << '\t' <<
            alignment.get_query_length() << '\t' <<
            alignment.get_query_start() << '\t' <<
            alignment.get_query_stop() << '\t' <<
            alignment.get_path_length() << '\t' <<
            alignment.get_path_start() << '\t' <<
            alignment.get_path_stop() << '\t' <<
            alignment.get_reversal() << '\t' <<
            alignment.get_n_match() << '\t' <<
            alignment.get_alignment_length() << '\t' <<
            alignment.get_map_quality() << '\t' <<
            '\n';
    });

    return 0;
}