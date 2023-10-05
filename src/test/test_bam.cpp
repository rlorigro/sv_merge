#include "Filesystem.hpp"
#include "Region.hpp"
#include "misc.hpp"
#include "bam.hpp"

using ghc::filesystem::path;
using sv_merge::for_cigar_tuple_in_alignment;
using sv_merge::for_cigar_interval_in_alignment;
using sv_merge::for_alignment_in_bam_region;
using sv_merge::decompress_bam_sequence;
using sv_merge::for_read_in_bam_region;
using sv_merge::reverse_complement;
using sv_merge::cigar_code_to_char;
using sv_merge::run_command;
using sv_merge::files_equal;
using sv_merge::CigarInterval;
using sv_merge::CigarTuple;
using sv_merge::Sequence;
using sv_merge::Region;

#include <stdexcept>
#include <iostream>
#include <fstream>
#include <map>

using std::runtime_error;
using std::to_string;
using std::ofstream;
using std::string;
using std::cerr;
using std::map;


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
    path output_path = data_directory / "test_bam_region_extract.fasta";

    ofstream file(output_path);

    for_alignment_in_bam_region(bam_path, region_string, [&](const bam1_t *alignment){
        cerr << bam_get_qname(alignment) << '\n';
        for_cigar_tuple_in_alignment(alignment, [&](CigarTuple& cigar_tuple){
            cerr << cigar_code_to_char[cigar_tuple.code] << ' ' << cigar_tuple.length << '\n';
        });
    });


}


void test_cigar_interval_iterator(path data_directory){
    string region_string = "a:0-5000";

    Region r(region_string);

    cerr << r.name << ' ' << r.start << ' ' << r.stop << '\n';

    path bam_path = data_directory / "test_alignment_softclip_only_sorted.bam";
    path output_path = data_directory / "test_bam_region_extract.fasta";

    ofstream file(output_path);

    map<string,string> results;
    map<string,string> formatted_query;
    map<string,string> formatted_ref;

    for_alignment_in_bam_region(bam_path, region_string, [&](const bam1_t* alignment){
        string name = bam_get_qname(alignment);

        string query_sequence;
        decompress_bam_sequence(alignment, query_sequence);

        for_cigar_interval_in_alignment(alignment, [&](CigarInterval& cigar){
            auto result = string(1,cigar_code_to_char[cigar.code]) + " " + to_string(cigar.length) + \
            '\t' + "r:" + to_string(cigar.ref_start) + "," + to_string(cigar.ref_stop) + \
            '\t' + "q:" + to_string(cigar.query_start) + "," + to_string(cigar.query_stop) + '\n';

            cerr << '\n' << name << '\n';
            cerr << alignment->core.l_qseq << '\n';
            cerr << result;

            auto [a,b] = cigar.get_forward_query_interval();
            auto s = query_sequence.substr(a, b-a + 1);

            if (cigar.is_reverse){
                reverse_complement(s);
            }

            formatted_query[name] += s;

            results[name] += result;
        });
    });

    for (const auto& [name,result]: results){
        cerr << name << '\n';
        cerr << result << '\n';
        cerr << formatted_query[name] << '\n';
        cerr << '\n';
    }

}


int main(){
    path project_directory = path(__FILE__).parent_path().parent_path().parent_path();
    path data_directory = project_directory / "data";

    cerr << data_directory << '\n';

    cerr << "TESTING: bam_sequence_extraction\n\n";
    test_bam_sequence_extraction(data_directory);

    cerr << "TESTING: bam_sequence_extraction\n\n";
    test_cigar_iterator(data_directory);

    cerr << "TESTING: bam_sequence_extraction\n\n";
    test_cigar_interval_iterator(data_directory);


    return 0;
}

