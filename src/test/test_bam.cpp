#include "cpptrace/from_current.hpp"
#include "Authenticator.hpp"
#include <filesystem>
#include "Alignment.hpp"
#include "Region.hpp"
#include "fasta.hpp"
#include "fetch.hpp"
#include "misc.hpp"
#include "bam.hpp"

using std::filesystem::path;
using namespace sv_merge;

#include <unordered_map>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <cmath>
#include <span>
#include <map>

using std::unordered_map;
using std::runtime_error;
using std::to_string;
using std::ofstream;
using std::string;
using std::cerr;
using std::span;
using std::map;
using std::max;


void test_bam_sequence_extraction(path data_directory){
    string region_string = "a:0-5000";

    Region r(region_string);

    cerr << r.name << ' ' << r.start << ' ' << r.stop << '\n';

    path bam_path = data_directory / "test_alignment_softclip_only_sorted.bam";
    path output_path = data_directory / "test_bam_region_extract.fasta";

    ofstream file(output_path);

    // BAM file is iterated (region 0:5000) using our libraries
    for_read_in_bam_region(bam_path, region_string, [&](Sequence& s){
        cerr << "---" << '\n';
        cerr << s.name << ',' << s.sequence.size() << '\n';
        file << '>' << s.name << '\n';
        file << s.sequence << '\n';
    });

    file.close();

    path region_output_bam_path = data_directory / "region.sam";

    string command;
    string result;

    // BAM file is viewed (0:5000) with samtools
    // Now part of data directory on github repo to avoid samtools dependency for testing
    command = "samtools view -h " + bam_path.string() + " " + region_string;
    run_command(command, region_output_bam_path);

    path samtools_output_fasta_path = data_directory / "region.fasta";

    // BAM file is converted to FASTA with samtools
    // Now part of data directory on github repo to avoid samtools dependency for testing
    command = "samtools fasta -0 " + samtools_output_fasta_path.string() + " " + region_output_bam_path.string();
    run_command(command);

    map<string,string> expected_sequences;
    for_sequence_in_fasta_file(samtools_output_fasta_path, [&](const Sequence& s){
        expected_sequences[s.name] = s.sequence;
    });

    map<string,string> result_sequences;
    for_sequence_in_fasta_file(output_path, [&](const Sequence& s){
        cerr << "result: " << output_path << ',' << s.name << ',' << s.sequence.size() << '\n';
        result_sequences[s.name] = s.sequence;
    });

    for (auto& [name, sequence]: expected_sequences){
        auto result2 = result_sequences.find(name);

        if (result2 == result_sequences.end()){
            throw runtime_error("FAIL: expected sequence not found in result sequences");
        }

        if (sequence == result2->second){
            result_sequences.erase(name);
            cerr << "PASS: " << name << '\n';
        }
        else{
            cerr << "expected: " << sequence.size() << '\n';
            cerr << "resulted: " << result2->second.size() << '\n';
            throw runtime_error("FAIL: expected sequence not equivalent to result sequence: "  + name);
        }
    }

    if (not result_sequences.empty()){
        throw runtime_error("FAIL: unexpected sequence in results");
    }

    cerr << '\n';
}


void test_bam_prefetched_subsequence_extraction(path data_directory) {
    string region_string = "a:0-5000";

    Region r(region_string);

    cerr << r.name << ' ' << r.start << ' ' << r.stop << '\n';

    path bam_path = data_directory / "test_alignment_hardclipped_sorted.bam";
    path output_path = data_directory / "test_alignment_hardclipped_sorted.fasta";

    ofstream file(output_path);
    if (not file.good() or not file.is_open()){
        throw runtime_error("ERROR: file could not be written: " + output_path.string());
    }

    cerr << "pre-fetching sequences..." << '\n';
    map<string, string> sequences;
    for_read_in_bam(bam_path, [&](Sequence& sequence){
        file << ">" << sequence.name << '\n';
        file << sequence.sequence << '\n';
        sequences.emplace(sequence.name, sequence.sequence);
    });

    vector<Region> subregions = {
            {"a",2000,2010},
            {"a",3090,4000},
            {"a",4990,5000},
            {"a",60000,6010},
            {"a",0,10}
    };

    // How to sort intervals
    auto left_comparator = [](const Region& a, const Region& b){
        return a.start < b.start;
    };

    sort(subregions.begin(), subregions.end(), left_comparator);

    Authenticator authenticator;
    sample_region_coord_map_t sample_to_region_coords;

    cerr << "computing coordinates..." << '\n';

    string sample_name = "A";

    // Initialize every combo of sample,region with an empty vector
    for (const auto& region: subregions){
        sample_to_region_coords[sample_name][region] = {};
    }

    extract_subregion_coords_from_sample(
            authenticator,
            sample_to_region_coords,
            sample_name,
            subregions,
            true,
            true,
            bam_path
    );

    for (const auto& region: subregions){
        cerr << sample_name << ' ' << region.to_string() << '\n';

        for (const auto& [name,coord]: sample_to_region_coords.at(sample_name).at(region)){
            cerr << '\t' << name << ' ' << coord.query_start << ',' << coord.query_stop << '\n';
            auto i = coord.query_start;
            auto l = coord.query_stop - coord.query_start;

            string s;
            if (coord.is_reverse) {
                s = sequences.at(name).substr(i, l);
                reverse_complement(s);
            }
            else{
                s = sequences[name].substr(i, l);
            }

            cerr << '\t' << s << '\n';
            cerr << '\n';
        }
    }

}


void test_bam_subsequence_extraction(path data_directory){
    string region_string = "a:0-5000";

    Region r(region_string);

    cerr << r.name << ' ' << r.start << ' ' << r.stop << '\n';

    path bam_path = data_directory / "test_alignment_softclip_only_sorted.bam";

    map<string,string> result_sequences;
    for_alignment_in_bam_region(bam_path, region_string, [&](Alignment& alignment){
        if (alignment.is_unmapped() or not alignment.is_primary()) {
            return;
        }

        Sequence s;

        alignment.get_query_name(s.name);
        alignment.get_query_sequence(s.sequence, 10, 110);
        result_sequences[s.name] = s.sequence;
    });

    path samtools_output_bam_path = data_directory / "region.sam";

    string command;
    string result;

    command = "samtools view -h " + bam_path.string() + " " + region_string;
    run_command(command, samtools_output_bam_path);

    path samtools_output_fasta_path = data_directory / "region.fasta";

    command = "samtools fasta -0 " + samtools_output_fasta_path.string() + " " + samtools_output_bam_path.string();
    run_command(command);

    map<string,string> expected_sequences = {
        {"del_5_at_34", "CATAGCCTATGACAGACTACAAGATATGAACGATCAATATGGTTCGCGACAGAATCGCTGTCCTAGTGGTTACGCGCTTACGGAAGTCGGACCACCCTAT"},
        {"del_5_at_34_reverse", "ATAGGGTGGTCCGACTTCCGTAAGCGCGTAACCACTAGGACAGCGATTCTGTCGCGAACCATATTGATCGTTCATATCTTGTAGTCTGTCATAGGCTATG"},
        {"del_5_at_34_ins_5_at_49", "CATAGCCTATGACAGACTACAAGATATGAACGATGGGGGCAATATGGTTCGCGACAGAATCGCTGTCCTAGTGGTTACGCGCTTACGGAAGTCGGACCAC"},
        {"del_5_at_34_ins_5_at_49_reverse", "GTGGTCCGACTTCCGTAAGCGCGTAACCACTAGGACAGCGATTCTGTCGCGAACCATATTGCCCCCATCGTTCATATCTTGTAGTCTGTCATAGGCTATG"},
        {"exact_match", "CATAGCCTATGACAGACTACAAGACCCACTATGAACGATCAATATGGTTCGCGACAGAATCGCTGTCCTAGTGGTTACGCGCTTACGGAAGTCGGACCAC"},
        {"exact_match_reverse", "GTGGTCCGACTTCCGTAAGCGCGTAACCACTAGGACAGCGATTCTGTCGCGAACCATATTGATCGTTCATAGTGGGTCTTGTAGTCTGTCATAGGCTATG"},
        {"ins_5_at_39", "CATAGCCTATGACAGACTACAAGACCCACGGGGGTATGAACGATCAATATGGTTCGCGACAGAATCGCTGTCCTAGTGGTTACGCGCTTACGGAAGTCGG"},
        {"ins_5_at_39_reverse", "CCGACTTCCGTAAGCGCGTAACCACTAGGACAGCGATTCTGTCGCGAACCATATTGATCGTTCATACCCCCGTGGGTCTTGTAGTCTGTCATAGGCTATG"},
        {"mismatch_1_at_49", "CATAGCCTATGACAGACTACAAGACCCACTATGAACGATGAATATGGTTCGCGACAGAATCGCTGTCCTAGTGGTTACGCGCTTACGGAAGTCGGACCAC"},
        {"mismatch_1_at_49_reverse", "GTGGTCCGACTTCCGTAAGCGCGTAACCACTAGGACAGCGATTCTGTCGCGAACCATATTCATCGTTCATAGTGGGTCTTGTAGTCTGTCATAGGCTATG"},
        {"softclip_2000", "TCGTGAGAACCTCGACCGCAGTGTGCTAATTTTTCAGCAACTTGTTCAACGAAGAACTTCCAAGTCTTAAGGTCGTAGGACTGCCGCGGCTGCCCATGTG"},
        {"softclip_2000_reverse", "CACATGGGCAGCCGCGGCAGTCCTACGACCTTAAGACTTGGAAGTTCTTCGTTGAACAAGTTGCTGAAAAATTAGCACACTGCGGTCGAGGTTCTCACGA"},
        {"supplementary_softclip_2000_at_1000", "ATTTAAATGCTTTCTTAGGGTGACAAGCCGACCACTGACCCGTGCTGAGGAATTCTCGCTTACAGTCGTCTGTCTGGAGCTCATCCCGCTGCGACAATTA"},
        {"supplementary_softclip_2000_at_1000_reverse", "TAATTGTCGCAGCGGGATGAGCTCCAGACAGACGACTGTAAGCGAGAATTCCTCAGCACGGGTCAGTGGTCGGCTTGTCACCCTAAGAAAGCATTTAAAT"}
    };

    for (auto& [name, sequence]: expected_sequences){
        cerr << name << '\n';
        auto result2 = result_sequences.find(name);

        if (result2 == result_sequences.end()){
            throw runtime_error("FAIL: expected sequence not found in result sequences");
        }

        if (sequence == result2->second){
            result_sequences.erase(name);
            cerr << "PASS" << '\n';
        }
        else{
            cerr << sequence << '\n';
            cerr << result2->second << '\n';
            throw runtime_error("FAIL: expected sequence not equivalent to result sequence");
        }
    }

    if (not result_sequences.empty()){
        throw runtime_error("FAIL: unexpected sequence in results");
    }

    cerr << '\n';
}


void test_clipped_bam_subsequence_extraction(path data_directory){
    string region_string = "a:0-5000";

    Region r(region_string);

    cerr << r.name << ' ' << r.start << ' ' << r.stop << '\n';

    path bam_path = data_directory / "test_alignment_hardclipped_sorted.bam";

    Timer t;

    vector<Region> regions = {
            {"a",2000,2010},
            {"a",3090,4000},
            {"a",4990,5000},
            {"a",60000,6010},
            {"a",0,10}
    };

    // How to sort intervals
    auto left_comparator = [](const Region& a, const Region& b){
        return a.start < b.start;
    };

    sort(regions.begin(), regions.end(), left_comparator);

    path bam_csv = data_directory / "test_bams.csv";

    ofstream file(bam_csv);

    if (not file.good() or not file.is_open()){
        throw runtime_error("ERROR: file could not be written: " + bam_csv.string());
    }

    file << "A," << bam_path.string() << '\n';
    file.close();

    int64_t n_threads = 4;

    unordered_map<Region,TransMap> region_transmaps;

    fetch_reads_from_clipped_bam(
            t,
            regions,
            bam_csv,
            n_threads,
            100'000,
            0,
            region_transmaps,
            true
            );

    cerr << '\n';

    string sequence;

    // Iterate sample reads
    for (const auto& region: regions){
        // Get the sample-read-path transitive mapping for this region
        auto& transmap = region_transmaps.at(region);

        // Iterate samples within this region and cluster reads
        transmap.for_each_sample([&](const string& sample_name, int64_t sample_id) {
            cerr << '\n';
            cerr << sample_name << ' ' << region.to_string() << '\n';

            transmap.for_each_read_of_sample(sample_name, [&](const string &name, int64_t id) {
                cerr << '>' << name << '\n';

                transmap.get_sequence(id, sequence);
                cerr << sequence << '\n';
            });
        });
    }
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

    bool unclip_coords = false;

    for_alignment_in_bam_region(bam_path, region_string, [&](Alignment& alignment){
        string name;
        alignment.get_query_name(name);

        string query_sequence;
        alignment.get_query_sequence(query_sequence);

        string alignment_name = name + "_" + to_string(counter[name]);
        counter[name]++;

        alignment.for_each_cigar_interval(unclip_coords, [&](const CigarInterval& cigar){
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


void test_for_alignment_in_bam_subregions(path data_directory){
    string region_string = "a:0-5000";

    vector<Region> subregions = {
            {"a",2000,2010},
            {"a",3090,4000},
            {"a",4990,5000},
            {"a",60000,6010},
            {"a",0,10}
    };

    // How to sort intervals
    auto left_comparator = [](const Region& a, const Region& b){
        return a.start < b.start;
    };

    sort(subregions.begin(), subregions.end(), left_comparator);

    Region r(region_string);

    cerr << r.name << ' ' << r.start << ' ' << r.stop << '\n';

    path bam_path = data_directory / "test_alignment_softclip_only_sorted.bam";
    path ref_path = data_directory / "test_ref.fasta";

    for_alignment_in_bam_subregions(
            bam_path,
            region_string,
            subregions,
            [&](Alignment& alignment, span<const Region>& overlapping_regions){

        string name;
        alignment.get_query_name(name);
        cerr << name << ' ' << alignment.get_ref_start() << ' ' << alignment.get_ref_stop() << '\n';
        for (const auto& item: overlapping_regions){
            cerr << item.start << ',' << item.stop << '\n';
        }
    });
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

        bool unclip_coords = false;

        int64_t length;

        interval_t prev_ref_interval = {-1,-1};
        interval_t prev_query_interval = {-1,-1};

        for_cigar_interval_in_alignment(unclip_coords, alignment, ref_intervals, query_intervals,
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
    CPPTRACE_TRY {
        path project_directory = path(__FILE__).parent_path().parent_path().parent_path();
        path data_directory = project_directory / "data";

        cerr << data_directory << '\n';

        cerr << "TESTING: bam_sequence_extraction\n\n";
        test_bam_sequence_extraction(data_directory);

        cerr << "TESTING: test_bam_subsequence_extraction\n\n";
        test_bam_subsequence_extraction(data_directory);

        cerr << "TESTING: test_bam_prefetched_subsequence_extraction\n\n";
        test_bam_prefetched_subsequence_extraction(data_directory);

        cerr << "TESTING: test_clipped_bam_subsequence_extraction\n\n";
        test_clipped_bam_subsequence_extraction(data_directory);

        cerr << "TESTING: cigar iterator\n\n";
        test_cigar_iterator(data_directory);

        cerr << "TESTING: cigar interval iterator\n\n";
        test_cigar_interval_iterator(data_directory);

        cerr << "TESTING: windowed cigar interval iterator\n\n";
        test_windowed_cigar_interval_iterator(data_directory);

        cerr << "TESTING: windowed cigar interval iterator\n\n";
        test_for_alignment_in_bam_subregions(data_directory);

    }
    CPPTRACE_CATCH(const std::exception& e) {
        std::cerr<<"Exception: "<<e.what()<<std::endl;
        cpptrace::from_current_exception().print_with_snippets();
        throw runtime_error("FAIL: exception caught");
    }

    return 0;
}

