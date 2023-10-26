#include "Filesystem.hpp"
#include "gaf.hpp"

using ghc::filesystem::path;
using sv_merge::for_alignment_in_gaf;
using sv_merge::cigar_code_to_char;
using sv_merge::GafAlignment;
using sv_merge::CigarTuple;

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
            alignment.is_reverse() << '\t' <<
            alignment.get_n_match() << '\t' <<
            alignment.get_alignment_length() << '\t' <<
            alignment.get_map_quality() << '\t';

            alignment.for_each_cigar_tuple([&](const CigarTuple& cigar){
                cerr << cigar.length << cigar_code_to_char[cigar.code] << ' ';
            });
            cerr << '\n';

    });

    return 0;
}
