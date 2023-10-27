#include "Filesystem.hpp"
#include "gaf.hpp"

using ghc::filesystem::path;
using sv_merge::for_alignment_in_gaf;
using sv_merge::cigar_code_to_char;
using sv_merge::GafAlignment;
using sv_merge::CigarTuple;

#include <unordered_map>
#include <iostream>
#include <stdexcept>

using std::unordered_map;
using std::runtime_error;
using std::cerr;


int main(){
    path project_directory = path(__FILE__).parent_path().parent_path().parent_path();
    path data_directory = project_directory / "data";

    path gaf_path = data_directory / "test_constructed_alignments.gaf";
    unordered_map<string,GafAlignment> expected_results;

    // Col 	Type 	Description
    // 1 	string 	Query sequence name
    // 2 	int 	Query sequence length
    // 3 	int 	Query start (0-based; closed)
    // 4 	int 	Query end (0-based; open)
    // 5 	char 	Strand relative to the path: "+" or "-"
    // 6 	string 	Path matching /([><][^\s><]+(:\d+-\d+)?)+|([^\s><]+)/
    // 7 	int 	Path length
    // 8 	int 	Start position on the path (0-based)
    // 9 	int 	End position on the path (0-based)
    // 10 	int 	Number of residue matches
    // 11 	int 	Alignment block length
    // 12 	int 	Mapping quality (0-255; 255 for missing)
    //
    GafAlignment a;
    a.set_query_name("F_R");
    a.set_query_length(100);
    a.set_query_start(0);
    a.set_query_stop(100);
    a.set_reversal(false);
    a.set_path({{"1",true},{"2",true},{"3",true}});
    a.set_path_length(100);
    a.set_path_start(0);
    a.set_path_stop(100);
    a.set_n_match(99);
    a.set_alignment_length(100);
    a.set_map_quality(60);
    a.set_is_primary(true);

    expected_results[a.get_query_name()] = a;

    a.set_query_name("R_R");
    a.set_query_length(100);
    a.set_query_start(0);
    a.set_query_stop(100);
    a.set_reversal(true);
    a.set_path({{"3",false},{"2",false},{"1",false}});
    a.set_path_length(100);
    a.set_path_start(0);
    a.set_path_stop(100);
    a.set_n_match(99);
    a.set_alignment_length(100);
    a.set_map_quality(60);
    a.set_is_primary(true);

    expected_results[a.get_query_name()] = a;

    a.set_query_name("F_F");
    a.set_query_length(100);
    a.set_query_start(0);
    a.set_query_stop(100);
    a.set_reversal(false);
    a.set_path({{"1",false},{"2",false},{"3",false}});
    a.set_path_length(100);
    a.set_path_start(0);
    a.set_path_stop(100);
    a.set_n_match(99);
    a.set_alignment_length(100);
    a.set_map_quality(60);
    a.set_is_primary(true);

    expected_results[a.get_query_name()] = a;

    a.set_query_name("R_F");
    a.set_query_length(100);
    a.set_query_start(0);
    a.set_query_stop(100);
    a.set_reversal(true);
    a.set_path({{"3",true},{"2",true},{"1",true}});
    a.set_path_length(100);
    a.set_path_start(0);
    a.set_path_stop(100);
    a.set_n_match(99);
    a.set_alignment_length(100);
    a.set_map_quality(60);
    a.set_is_primary(false);

    expected_results[a.get_query_name()] = a;

    a.set_query_name("no_tags");
    a.set_query_length(100);
    a.set_query_start(0);
    a.set_query_stop(100);
    a.set_reversal(true);
    a.set_path({{"3",true},{"2",true},{"1",true}});
    a.set_path_length(100);
    a.set_path_start(0);
    a.set_path_stop(100);
    a.set_n_match(99);
    a.set_alignment_length(100);
    a.set_map_quality(60);
    a.set_is_primary(false);

    expected_results[a.get_query_name()] = a;

    for_alignment_in_gaf(gaf_path, [&](GafAlignment& alignment){
        auto result = expected_results.find(alignment.get_query_name());

        if (result == expected_results.end()){
            throw runtime_error("FAIL: alignment not expected: " + alignment.get_query_name());
        }

        auto& e = result->second;

        if (alignment.get_query_name() != e.get_query_name()){
            throw runtime_error("FAIL: alignment.get_query_name does not match expected result for " + alignment.get_query_name());
        }
        if (alignment.get_ref_name() != e.get_ref_name()){
            throw runtime_error("FAIL: alignment.get_ref_name does not match expected result for " + alignment.get_query_name());
        }
        if (alignment.get_query_length() != e.get_query_length()){
            throw runtime_error("FAIL: alignment.get_query_length does not match expected result for " + alignment.get_query_name());
        }
        if (alignment.get_query_start() != e.get_query_start()){
            throw runtime_error("FAIL: alignment.get_query_start does not match expected result for " + alignment.get_query_name());
        }
        if (alignment.get_query_stop() != e.get_query_stop()){
            throw runtime_error("FAIL: alignment.get_query_stop does not match expected result for " + alignment.get_query_name());
        }
        if (alignment.get_path_length() != e.get_path_length()){
            throw runtime_error("FAIL: alignment.get_path_length does not match expected result for " + alignment.get_query_name());
        }
        if (alignment.get_path_string() != e.get_path_string()){
            cerr << alignment.get_path_string() << '\n';
            cerr << e.get_path_string() << '\n';
            throw runtime_error("FAIL: alignment.get_path_length does not match expected result for " + alignment.get_query_name());
        }
        if (alignment.get_path_start() != e.get_path_start()){
            throw runtime_error("FAIL: alignment.get_path_start does not match expected result for " + alignment.get_query_name());
        }
        if (alignment.get_path_stop() != e.get_path_stop()){
            cerr << alignment.get_path_stop() << '\n';
            cerr << e.get_path_stop() << '\n';
            throw runtime_error("FAIL: alignment.get_path_stop does not match expected result for " + alignment.get_query_name());
        }
        if (alignment.is_reverse() != e.is_reverse()){
            throw runtime_error("FAIL: alignment.is_reverse does not match expected result for " + alignment.get_query_name());
        }
        if (alignment.get_n_match() != e.get_n_match()){
            throw runtime_error("FAIL: alignment.get_n_match does not match expected result for " + alignment.get_query_name());
        }
        if (alignment.get_alignment_length() != e.get_alignment_length()){
            throw runtime_error("FAIL: alignment.get_alignment_length does not match expected result for " + alignment.get_query_name());
        }
        if (alignment.get_map_quality() != e.get_map_quality()){
            throw runtime_error("FAIL: alignment.get_map_quality does not match expected result for " + alignment.get_query_name());
        }
        if (alignment.is_primary() != e.is_primary()){
            throw runtime_error("FAIL: alignment.is_primary does not match expected result for " + alignment.get_query_name());
        }

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
            alignment.get_map_quality() << '\t' <<
            alignment.is_primary() << '\t';

            alignment.for_each_cigar_tuple([&](const CigarTuple& cigar){
                cerr << cigar.length << cigar_code_to_char[cigar.code] << ' ';
            });
            cerr << '\n';

    });

    cerr << "PASS" << '\n';

    return 0;
}
