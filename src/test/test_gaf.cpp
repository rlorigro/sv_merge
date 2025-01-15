#include <filesystem>
using std::filesystem::path;

#include "Sequence.hpp"
#include "fasta.hpp"
#include "bed.hpp"
#include "gaf.hpp"

using sv_merge::for_sequence_in_fasta_file;
using sv_merge::for_alignment_in_gaf;
using sv_merge::reverse_complement;
using sv_merge::cigar_code_to_char;
using sv_merge::CigarInterval;
using sv_merge::is_ref_move;
using sv_merge::is_query_move;
using sv_merge::GafAlignment;
using sv_merge::CigarTuple;
using sv_merge::interval_t;
using sv_merge::Sequence;
using sv_merge::Region;

#include <unordered_map>
#include <iostream>
#include <stdexcept>
#include <map>

using std::unordered_map;
using std::runtime_error;
using std::to_string;
using std::cerr;
using std::map;


void test_windowed_cigar_interval_iterator(path data_directory){
    string region_string = "a:0-5000";

    Region r(region_string);

    cerr << r.name << ' ' << r.start << ' ' << r.stop << '\n';

    path gaf_path = data_directory / "test_alignments_cigar.gaf";
    path ref_path = data_directory / "test_ref.fasta";
    path query_path = data_directory / "test_query.fasta";

    Sequence ref_sequence;

    Sequence ref_sequence_forward;
    for_sequence_in_fasta_file(ref_path, [&](const Sequence& s){
        if (s.name == "a"){
            ref_sequence_forward = s;
        }
    });

    Sequence ref_sequence_reverse = ref_sequence_forward;
    reverse_complement(ref_sequence_reverse.sequence);

    unordered_map<string,string> query_sequences;
    for_sequence_in_fasta_file(query_path, [&](const Sequence& s){
        cerr << s.name << '\n';
        query_sequences[s.name] = s.sequence;
    });

    unordered_map<string,int> counter;

    vector<interval_t> ref_intervals = {
            {1,50},                 // out of range
            {2000,2036},
            {2036,2049},
            {2049,2050},
            {2050,2100},
            {4000-100,4000-46},     // 3900,3956
            {4000-46,4000},         // 3956,4000
            {7000-40,7000}          // out of range
    };

    vector<interval_t> query_intervals = {
            {-10,-5},               // out of range
            {0,2},
            {2,2},                  // zero-length interval
            {2,3},                  // one-length interval
            {3,34},                 // Terminates left of DEL
            {34,35},                // Terminates right of DEL
            {35,44},                // Terminates left of INS
            {44,46},                // Terminates middle of INS
            {44,49},                // Terminates right of INS
            {49,100},
            {200,210},
            {2000-100,2000-46},     // 1900,1956
            {2000-46,2000},         // 1956,2000
            {7000-40,7000}          // out of range
    };

    bool unclip_coords = false;

    string s_ref;
    string s_query;
    string s_crossref;

    map<string,string> r_formatted_query;
    map<string,string> r_formatted_ref;
    map<string,string> r_formatted_crossref;

    map<string,string> q_formatted_query;
    map<string,string> q_formatted_ref;
    map<string,string> q_formatted_crossref;

    for_alignment_in_gaf(gaf_path, [&](GafAlignment& alignment){
        string name;
        alignment.get_query_name(name);

        // In practice, this will need to be fetched according to which interval on the path is active.
        // For the purpose of this test there is only one ref node, as with a linear alignment
        const string& query_sequence = query_sequences.at(name);

        string alignment_name = name + "_" + to_string(counter[name]);
        counter[name]++;

        cerr << alignment_name << '\n';

        int64_t length;

        interval_t prev_ref_interval = {-1,-1};
        interval_t prev_query_interval = {-1,-1};
        CigarInterval intersection;

        for_cigar_interval_in_alignment(unclip_coords, alignment, ref_intervals, query_intervals,
            [&](const CigarInterval& i, const interval_t& interval){
                // BAM stores query sequences in reverse, and always aligns reads in ref F orientation.
                // GAF does not store any query sequences in reverse, and in general, the path involves reversed "ref" sequences
                // TODO: properly address this in the intervals
                intersection = i;
                if (alignment.get_path_step(0).second){
                    ref_sequence = ref_sequence_reverse;
                }
                else{
                    ref_sequence = ref_sequence_forward;
                }

                if (is_ref_move[intersection.code]){
                    length = intersection.ref_stop - intersection.ref_start;
                }
                else{
                    length = intersection.get_query_length();
                }
                cerr << "r:" << interval.first << ',' << interval.second << ' ' << cigar_code_to_char[intersection.code] << ',' << intersection.length << ',' << length << " r:" << intersection.ref_start << ',' << intersection.ref_stop << " q:" << intersection.query_start << ',' << intersection.query_stop << '\n';

                if (not intersection.is_clip()){
                    get_formatted_sequence_of_cigar_interval(
                            intersection,
                            ref_sequence.sequence,
                            query_sequence,
                            s_ref,
                            s_query,
                            s_crossref);

                    if (interval != prev_ref_interval and prev_ref_interval.first >= 0 and prev_ref_interval.second >= 0){
                        r_formatted_ref[alignment_name] += " ";
                        r_formatted_crossref[alignment_name] += " ";
                        r_formatted_query[alignment_name] += " ";
                    }

                    r_formatted_ref[alignment_name] += s_ref;
                    r_formatted_crossref[alignment_name] += s_crossref;
                    r_formatted_query[alignment_name] += s_query;
                }

                prev_ref_interval = interval;
            },
            [&](const CigarInterval& i, const interval_t& interval){
                // BAM stores query sequences in reverse, and always aligns reads in ref F orientation.
                // GAF does not store any query sequences in reverse, and in general, the path involves reversed "ref" sequences
                // TODO: properly address this in the intervals
                intersection = i;
                if (alignment.get_path_step(0).second){
                    ref_sequence = ref_sequence_reverse;
                }
                else{
                    ref_sequence = ref_sequence_forward;
                }

                if (is_query_move[intersection.code]){
                    length = intersection.query_stop - intersection.query_start;
                }
                else{
                    length = intersection.get_ref_length();
                }
                cerr << "q:" << interval.first << ',' << interval.second << ' ' << cigar_code_to_char[intersection.code] << ',' << intersection.length << ',' << length << " r:" << intersection.ref_start << ',' << intersection.ref_stop << " q:" << intersection.query_start << ',' << intersection.query_stop << '\n';

                if (not intersection.is_clip()){
                    get_formatted_sequence_of_cigar_interval(
                            intersection,
                            ref_sequence.sequence,
                            query_sequence,
                            s_ref,
                            s_query,
                            s_crossref);

                    if (interval != prev_query_interval and prev_query_interval.first >= 0 and prev_query_interval.second >= 0){
                        q_formatted_ref[alignment_name] += " ";
                        q_formatted_crossref[alignment_name] += " ";
                        q_formatted_query[alignment_name] += " ";
                    }

                    q_formatted_ref[alignment_name] += s_ref;
                    q_formatted_crossref[alignment_name] += s_crossref;
                    q_formatted_query[alignment_name] += s_query;
                }

                prev_query_interval = interval;
            });

        cerr << '\n';

    });

    cerr << '\n';
    cerr << "TESTING REF INTERVALS" << '\n';
    for (const auto& [name,result]: r_formatted_query){
        cerr << name << '\n';

        cerr << r_formatted_ref[name] << '\n';
        cerr << r_formatted_crossref[name] << '\n';
        cerr << r_formatted_query[name] << '\n';

        cerr << '\n';
    }

    cerr << '\n';
    cerr << "TESTING QUERY INTERVALS" << '\n';
    for (const auto& [name,result]: q_formatted_query){
        cerr << name << '\n';

        cerr << q_formatted_ref[name] << '\n';
        cerr << q_formatted_crossref[name] << '\n';
        cerr << q_formatted_query[name] << '\n';

        cerr << '\n';
    }

}


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
    GafAlignment x;
    x.set_query_name("F_R");
    x.set_query_length(100);
    x.set_query_start(0);
    x.set_query_stop(100);
    x.set_reversal(false);
    x.set_path({{"1", true}, {"2", true}, {"3", true}});
    x.set_path_length(100);
    x.set_path_start(0);
    x.set_path_stop(100);
    x.set_n_match(99);
    x.set_alignment_length(100);
    x.set_map_quality(60);
    x.set_is_primary(true);

    expected_results[x.get_query_name()] = x;

    x.set_query_name("R_R");
    x.set_query_length(100);
    x.set_query_start(0);
    x.set_query_stop(100);
    x.set_reversal(true);
    x.set_path({{"3", false}, {"2", false}, {"1", false}});
    x.set_path_length(100);
    x.set_path_start(0);
    x.set_path_stop(100);
    x.set_n_match(99);
    x.set_alignment_length(100);
    x.set_map_quality(60);
    x.set_is_primary(true);

    expected_results[x.get_query_name()] = x;

    x.set_query_name("F_F");
    x.set_query_length(100);
    x.set_query_start(0);
    x.set_query_stop(100);
    x.set_reversal(false);
    x.set_path({{"1", false}, {"2", false}, {"3", false}});
    x.set_path_length(100);
    x.set_path_start(0);
    x.set_path_stop(100);
    x.set_n_match(99);
    x.set_alignment_length(100);
    x.set_map_quality(60);
    x.set_is_primary(true);

    expected_results[x.get_query_name()] = x;

    x.set_query_name("R_F");
    x.set_query_length(100);
    x.set_query_start(0);
    x.set_query_stop(100);
    x.set_reversal(true);
    x.set_path({{"3", true}, {"2", true}, {"1", true}});
    x.set_path_length(100);
    x.set_path_start(0);
    x.set_path_stop(100);
    x.set_n_match(99);
    x.set_alignment_length(100);
    x.set_map_quality(60);
    x.set_is_primary(false);

    expected_results[x.get_query_name()] = x;

    x.set_query_name("no_tags");
    x.set_query_length(100);
    x.set_query_start(0);
    x.set_query_stop(100);
    x.set_reversal(true);
    x.set_path({{"3", true}, {"2", true}, {"1", true}});
    x.set_path_length(100);
    x.set_path_start(0);
    x.set_path_stop(100);
    x.set_n_match(99);
    x.set_alignment_length(100);
    x.set_map_quality(60);
    x.set_is_primary(false);

    expected_results[x.get_query_name()] = x;

    vector<CigarTuple> expected_cigar_tuples_f = {
            {49,7},
            {1,8},
            {50,7},
    };

    vector<CigarTuple> expected_cigar_tuples_r = {
            {50,7},
            {1,8},
            {49,7},
    };

    for_alignment_in_gaf(gaf_path, [&](GafAlignment& a){
        cerr << "TESTING: " << a.get_query_name() << '\n';

        auto result = expected_results.find(a.get_query_name());

        if (result == expected_results.end()){
            throw runtime_error("FAIL: alignment not expected: " + a.get_query_name());
        }

        auto& e = result->second;

        if (a.get_query_name() != e.get_query_name()){
            throw runtime_error("FAIL: alignment.get_query_name does not match expected result for " + a.get_query_name());
        }
        if (a.get_ref_name() != e.get_ref_name()){
            throw runtime_error("FAIL: alignment.get_ref_name does not match expected result for " + a.get_query_name());
        }
        if (a.get_query_length() != e.get_query_length()){
            throw runtime_error("FAIL: alignment.get_query_length does not match expected result for " + a.get_query_name());
        }
        if (a.get_query_start() != e.get_query_start()){
            throw runtime_error("FAIL: alignment.get_query_start does not match expected result for " + a.get_query_name());
        }
        if (a.get_query_stop() != e.get_query_stop()){
            throw runtime_error("FAIL: alignment.get_query_stop does not match expected result for " + a.get_query_name());
        }
        if (a.get_path_length() != e.get_path_length()){
            throw runtime_error("FAIL: alignment.get_path_length does not match expected result for " + a.get_query_name());
        }
        if (a.get_path_string() != e.get_path_string()){
            cerr << a.get_path_string() << '\n';
            cerr << e.get_path_string() << '\n';
            throw runtime_error("FAIL: alignment.get_path_length does not match expected result for " + a.get_query_name());
        }
        if (a.get_path_start() != e.get_path_start()){
            throw runtime_error("FAIL: alignment.get_path_start does not match expected result for " + a.get_query_name());
        }
        if (a.get_path_stop() != e.get_path_stop()){
            cerr << a.get_path_stop() << '\n';
            cerr << e.get_path_stop() << '\n';
            throw runtime_error("FAIL: alignment.get_path_stop does not match expected result for " + a.get_query_name());
        }
        if (a.is_reverse() != e.is_reverse()){
            throw runtime_error("FAIL: alignment.is_reverse does not match expected result for " + a.get_query_name());
        }
        if (a.get_n_match() != e.get_n_match()){
            throw runtime_error("FAIL: alignment.get_n_match does not match expected result for " + a.get_query_name());
        }
        if (a.get_alignment_length() != e.get_alignment_length()){
            throw runtime_error("FAIL: alignment.get_alignment_length does not match expected result for " + a.get_query_name());
        }
        if (a.get_map_quality() != e.get_map_quality()){
            throw runtime_error("FAIL: alignment.get_map_quality does not match expected result for " + a.get_query_name());
        }
        if (a.is_primary() != e.is_primary()){
            throw runtime_error("FAIL: alignment.is_primary does not match expected result for " + a.get_query_name());
        }

        int i = 0;
        a.for_each_cigar_tuple([&](const CigarTuple& cigar){
            if (not a.is_reverse()){
                if (cigar.code != expected_cigar_tuples_f[i].code or cigar.length != expected_cigar_tuples_f[i].length){
                    cerr << int(cigar.code) << ',' << cigar.length << " != " << int(expected_cigar_tuples_f[i].code) << ',' << expected_cigar_tuples_f[i].length << '\n';
                    throw runtime_error("FAIL: cigar not expected");
                }
            }
            else{
                if (cigar.code != expected_cigar_tuples_r[i].code or cigar.length != expected_cigar_tuples_r[i].length){
                    cerr << int(cigar.code) << ',' << cigar.length << " != " << int(expected_cigar_tuples_r[i].code) << ',' << expected_cigar_tuples_r[i].length << '\n';
                    throw runtime_error("FAIL: cigar not expected");
                }
            }

            i++;
        });

        cerr <<
             a.get_path_string() << '\t' <<
             a.get_query_name() << '\t' <<
             a.get_ref_name() << '\t' <<
             a.get_query_length() << '\t' <<
             a.get_query_start() << '\t' <<
             a.get_query_stop() << '\t' <<
             a.get_path_length() << '\t' <<
             a.get_path_start() << '\t' <<
             a.get_path_stop() << '\t' <<
             a.is_reverse() << '\t' <<
             a.get_n_match() << '\t' <<
             a.get_alignment_length() << '\t' <<
             a.get_map_quality() << '\t' <<
             a.is_primary() << '\t';

            a.for_each_cigar_tuple([&](const CigarTuple& cigar){
                cerr << cigar.length << cigar_code_to_char[cigar.code] << ' ';
            });
            cerr << '\n';

    });

    cerr << "PASS exact tests" << '\n';

    cerr << "now testing output (to be manually inspected) ..." << '\n';

    test_windowed_cigar_interval_iterator(data_directory);

    return 0;
}
