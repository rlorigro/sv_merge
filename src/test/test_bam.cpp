#include "IntervalGraph.hpp"
#include "Filesystem.hpp"
#include "Alignment.hpp"
#include "Region.hpp"
#include "fasta.hpp"
#include "misc.hpp"
#include "bam.hpp"

using ghc::filesystem::path;
using sv_merge::for_cigar_interval_in_alignment;
using sv_merge::for_alignment_in_bam_region;
using sv_merge::for_sequence_in_fasta_file;
using sv_merge::decompress_bam_sequence;
using sv_merge::for_read_in_bam_region;
using sv_merge::reverse_complement;
using sv_merge::cigar_code_to_format_char;
using sv_merge::cigar_code_to_char;
using sv_merge::is_query_move;
using sv_merge::is_ref_move;
using sv_merge::run_command;
using sv_merge::files_equal;
using sv_merge::CigarInterval;
using sv_merge::CigarTuple;
using sv_merge::interval_t;
using sv_merge::Alignment;
using sv_merge::Sequence;
using sv_merge::Region;

#include <unordered_map>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <cmath>
#include <map>

using std::unordered_map;
using std::runtime_error;
using std::to_string;
using std::ofstream;
using std::string;
using std::cerr;
using std::map;
using std::max;


void test_bam_sequence_extraction(path data_directory){
    string region_string = "a:0-5000";

    Region r(region_string);

    cerr << r.name << ' ' << r.start << ' ' << r.stop << '\n';

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

}


void test_cigar_iterator(path data_directory){
    string region_string = "a:0-5000";

    Region r(region_string);

    cerr << r.name << ' ' << r.start << ' ' << r.stop << '\n';

    path bam_path = data_directory / "test_alignment_softclip_only_sorted.bam";

    string name;
    for_alignment_in_bam_region(bam_path, region_string, [&](Alignment& alignment){
        alignment.get_query_name(name);
        cerr << name << '\n';

        alignment.for_each_cigar_tuple([&](const CigarTuple& cigar_tuple) {
            cerr << cigar_code_to_char[cigar_tuple.code] << ' ' << cigar_tuple.length << '\n';
        });
    });
}


void get_formatted_sequence_of_cigar_interval(
        const CigarInterval& cigar,
        const string& ref_sequence,
        const string& query_sequence,
        string& s_ref,
        string& s_query,
        string& s_crossref
){
    s_ref.clear();
    s_query.clear();
    s_crossref.clear();

    if (not cigar.is_clip()){
        const auto& [a_query,b_query] = cigar.get_forward_query_interval();
        const auto& [a_ref,b_ref] = cigar.get_forward_ref_interval();

        auto l_query = b_query - a_query;
        auto l_ref = b_ref - a_ref;
        auto l = max(l_ref,l_query);

//        cerr << l_query << ',' << l_ref << '\n';

        s_query = query_sequence.substr(a_query, l_query);
        s_ref = ref_sequence.substr(a_ref, l_ref);

        if (cigar.is_reverse){
            reverse_complement(s_query);
        }

        // Append filler characters to keep lengths equal, if needed
        if (not is_query_move[cigar.code]){
          s_query = string(l, '-');
        }

        if (not is_ref_move[cigar.code]){
          s_ref = string(l, '-');
        }

        // Show whether there is a match, mismatch, or indel
        s_crossref = string(l, cigar_code_to_format_char[cigar.code]);
    }
}


void test_cigar_interval_iterator(path data_directory){
    string region_string = "a:0-5000";

    Region r(region_string);

    cerr << r.name << ' ' << r.start << ' ' << r.stop << '\n';

    path bam_path = data_directory / "test_alignment_softclip_only_sorted.bam";
    path ref_path = data_directory / "test_ref.fasta";

    Sequence ref_sequence;
    for_sequence_in_fasta_file(ref_path, [&](const Sequence& s){
        if (s.name == "a"){
            ref_sequence = s;
        }
    });

    string s_ref;
    string s_query;
    string s_crossref;

    map<string,string> results;
    map<string,string> formatted_query;
    map<string,string> formatted_ref;
    map<string,string> formatted_crossref;

    unordered_map<string,int> counter;

    for_alignment_in_bam_region(bam_path, region_string, [&](Alignment& alignment){
        string name;
        alignment.get_query_name(name);

        string query_sequence;
        alignment.get_query_sequence(query_sequence);

        string alignment_name = name + "_" + to_string(counter[name]);
        counter[name]++;

        alignment.for_each_cigar_interval([&](const CigarInterval& cigar){
            auto result = string(1,cigar_code_to_char[cigar.code]) + " " + to_string(cigar.length) + \
            '\t' + "r:" + to_string(cigar.ref_start) + "," + to_string(cigar.ref_stop) + \
            '\t' + "q:" + to_string(cigar.query_start) + "," + to_string(cigar.query_stop) + '\n';

            cerr << '\n' << name << '\n';
            cerr << alignment.get_query_length() << '\n';
            cerr << result;

            results[name] += result;

            if (not cigar.is_clip()){
                get_formatted_sequence_of_cigar_interval(
                        cigar,
                        ref_sequence.sequence,
                        query_sequence,
                        s_ref,
                        s_query,
                        s_crossref);

                formatted_ref[alignment_name] += s_ref;
                formatted_crossref[alignment_name] += s_crossref;
                formatted_query[alignment_name] += s_query;

            }
        });
    });

    for (const auto& [name,result]: results){
        cerr << name << '\n';
        cerr << result << '\n';
        cerr << '\n';
    }

    for (const auto& [name,result]: formatted_query){
        cerr << name << '\n';

        cerr << formatted_ref[name] << '\n';
        cerr << formatted_crossref[name] << '\n';
        cerr << formatted_query[name] << '\n';

        cerr << '\n';
    }

}


void test_windowed_cigar_interval_iterator(path data_directory){
    string region_string = "a:0-5000";

    Region r(region_string);

    cerr << r.name << ' ' << r.start << ' ' << r.stop << '\n';

    path bam_path = data_directory / "test_alignment_softclip_only_sorted.bam";
    path ref_path = data_directory / "test_ref.fasta";

    Sequence ref_sequence;
    for_sequence_in_fasta_file(ref_path, [&](const Sequence& s){
        if (s.name == "a"){
            ref_sequence = s;
        }
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

    string s_ref;
    string s_query;
    string s_crossref;

    map<string,string> r_formatted_query;
    map<string,string> r_formatted_ref;
    map<string,string> r_formatted_crossref;

    map<string,string> q_formatted_query;
    map<string,string> q_formatted_ref;
    map<string,string> q_formatted_crossref;

    for_alignment_in_bam_region(bam_path, region_string, [&](Alignment& alignment){
        string name;
        alignment.get_query_name(name);

        string query_sequence;
        alignment.get_query_sequence(query_sequence);

        string alignment_name = name + "_" + to_string(counter[name]);
        counter[name]++;

        cerr << alignment_name << '\n';

        int64_t length;

        interval_t prev_ref_interval = {-1,-1};
        interval_t prev_query_interval = {-1,-1};

        for_cigar_interval_in_alignment(alignment, ref_intervals, query_intervals,
        [&](const CigarInterval& intersection, const interval_t& interval){
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
        [&](const CigarInterval& intersection, const interval_t& interval){
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

    cerr << data_directory << '\n';

    cerr << "TESTING: bam_sequence_extraction\n\n";
    test_bam_sequence_extraction(data_directory);

    cerr << "TESTING: cigar iterator\n\n";
    test_cigar_iterator(data_directory);

    cerr << "TESTING: cigar interval iterator\n\n";
    test_cigar_interval_iterator(data_directory);

    cerr << "TESTING: windowed cigar interval iterator\n\n";
    test_windowed_cigar_interval_iterator(data_directory);


    return 0;
}

