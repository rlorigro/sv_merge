#include <iostream>
#include <stdexcept>

using std::runtime_error;
using std::cerr;

#include <iostream>
#include <stdexcept>

#include "SimpleAlignment.hpp"
#include "bindings/cpp/WFAligner.hpp"

using namespace wfa;

using namespace sv_merge;
using std::runtime_error;
using std::cerr;
using std::cout;

void test_simple_alignment() {
    // Define test sequences
    string ref_seq =   "agactacataggttcaagtctccactccagcatcactagctgattgctcaacatccagctactgagcctgtttcctaggagtacaatgaaaatataaggtgctttgtatagtgcctggtacacagcaagtgctccataagtgtttgctgctcttgggtttgGAGGACAATAATGGCAAACTTGCTCCCTAGCACTCATGTCCAGTGTAGACAAACCTGAAGAGGGTGTGTCTCCTCTGACCCTCTCAAAGTACACAGGTGAAAAGCAGGAGTGAGAATGTTTGGGTTCTGAGCTCACACATGGACCAAGCTTCCCTAGGCTGGCTTGGGGCTATGTCTGACTTCTGTTCTTAAAAAGCCAAttcatgggccgggcgcggtggctcacacctgcaatcccagcactttgggaggccgaggcgggcggatcacgaggtcaggagatagagaccatcctggctaacacggtgaaaccccgtctctactaaaaatacaaaaaattagcttggcgaggtggcgggcgcctgtagtcccagctacttgggaggctgaggcaggagaatggcgtgaaccccagggggcggagcctgcagtgagccgagatcacgccactgcactccagcctgggcgacagcgagactccgtctcaaaaaacaaacaaacaaacaaacaaacaaaaagccaattcatgctgggtgcagtggctcacacctataatcccagcactttgggaggctgaggtgggaggattgctttgagctcaggagagagcagcctgggcaacagggcaaaaacccatctctacaaaaaaaaaatacaaaaaattagccaggtatggtggcgtggtggctcacatctgtggttccagctactcaagaggctgaggcaggggaattacttgaacccaggaggtggaggttgccgtgagccgagattgcaccactgtacttcagcctgggtgacagagtgagaccctgtctcagaaaaaaaaaaaaaaaaaaaaaaaaaaaggccaattcAACAGCCATCTCCTGCAAGTTCTGTGATGCCAGGCcctcaccacaaccctgtgtggtaggtagcactgttcccatttgacagatgaggagactgaggctcagagaggtgacatggccaccttgaggccagaccgcCCCTTCCAGGCCAGGCTCACTAGGCTCCCCTGATGTGGAGCTGCAGGAAAGCACGTGAGGCCCAAG";
    string query_seq = "AGACTACATAGGTTCAAGTCTCCACTCCAGCATCACTAGCTGATTGCTCAACATCCAGCTACTGAGCCTGTTTCCTAGGAGTACAATGAAAATATAAGGTGCTTTGTATAGTGCCTGGTACACAGCAAGTGCTCCATAAGTGTTTGCTGCTCTTGGGTTTGGAGGACAATAATGGCAAACTTGCTCCCTAGCACTCATGTCCAGTGTAGACAAACCTGAAGAGGGTGTGTCTCCTCCTGACCCAGTGTAGACAAAGCTGAAGAGGGTGTGTCTCCTCCTGACCCTCTCAAAGTACACAGGTGAAAAGCAGGAGTGAGAATGTTTGGGTTCTGAGCTCACACATGGACCAAGCTTCCCTAGGCTGGCTTGGGGCTATGTCTGACTTCTGTTCTTAAAAAGCCAATTCATGGGCCGGGCGCGGTGGCTCACACCTGCAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAGATAGAGACCATCCTGGCTAACACGGTGAAACCCCGTCTCTACTAAAAATACAAAAAATTAGCTTGGCGAGGTGGCGGGCGCCTGTAGTCCCAGCTACTTGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCCAGGGGCGGAGCCTGCAGTGAGCCGAGATCACGCCACTGCACTCCAGCCTGGGCGACAGCGAGACTCCGTCTCAAAAAACAAACAAACAAACAAACAAACAAAAAGCCAATTCATGCTGGGTGCAGTGGCTCACACCTATAATCCCAGCACTTTGGGAGGCTGAGGTGGGAGGATTGCTTTGAGCTCAGGAGAGAGCAGCCTGGGCAACAGGGCAAAACCCATCTCTACAAAAAAAAAATACAAAAAATTAGCCAGGTATGGTGGCGTGGTGGCTCACATCTGTGGTTCCAGCTACTCAAGAGGCTGAGGCAGGGGAATTACTTGAACCCAGGAGGTGGAGGTTGCCGTGAGCCGAGATTGCACCACTGTACTTCAGCCTGGGTGACAGAGTGAGACCCTGTCTCAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGCCAATTCAACAGCCATCTCCTGCAAGTTCTGTGATGCCAGGCCCTCACCACAACCCTGTGTGGTAGGTAGCACTGTTCCCATTTGACAGATGAGGAGACTGAGGCTCAGAGAGGTGACATGGCCACCTTGAGGCCAGACCGCCCCTTCCAGGCCAGGCTCACTAGGCTCCCCTGATGTGGAGCTGCAGGAAAGCATGTGAGGCCCAAG";

    std::transform(ref_seq.begin(), ref_seq.end(), ref_seq.begin(), ::toupper);
    std::transform(query_seq.begin(), query_seq.end(), query_seq.begin(), ::toupper);

    // Use WFA2 to generate the CIGAR string
    WFAlignerGapAffine aligner(4,2,2,WFAligner::Alignment,WFAligner::MemoryHigh);

    aligner.alignEnd2End(ref_seq, query_seq);

    // Extract indel edit distance from cigar
    string cigar = aligner.getCIGAR(false);

    cerr << cigar << '\n';

//    cigar = "236M42I343M1D212M1D216M1I210M";
//    cerr << cigar << '\n';

    // Construct SimpleAlignment
    SimpleAlignment alignment(ref_seq, query_seq, cigar);

    // Define ref intervals to be iterated
//    std::vector<interval_t> ref_intervals = {{0,ref_seq.size()}};
    std::vector<interval_t> ref_intervals = {{170,1049}};
    std::vector<interval_t> query_intervals = {};

    string a;
    string a_all;
    string x;
    string x_all;
    string b;
    string b_all;

    // Test this function:
    /**
     * Iterate cigar intervals and return cigar intervals that are clipped to match the bounds of each window provided by
     * the user as a vector of interval_t. Useful for parsing cigars over a region of interest. WARNING: empty intervals
     * such as [2,2) are not ever returned by this iterator.
     * @param unclip_coords - re-interpret hardclips as softclips. Intended to fetch coords for native/unclipped sequence
     * @param alignment - pointer to htslib alignment struct
     * @param ref_intervals - intervals which MUST BE NON-OVERLAPPING and NO CHECK is performed to verify this! Intervals
     * must be half-open in F orientation, e.g.: [[0,2),[2,4)] are two adjacent intervals of length 2.
     * @param query_intervals - intervals which MUST BE NON-OVERLAPPING and NO CHECK is performed to verify this! Intervals
     * must be half-open in F orientation, e.g.: [[0,2),[2,4)] are two adjacent intervals of length 2.
     * @param f_ref lambda function to operate on ref intervals
     * @param f_query lambda function to operate on query intervals
     */
    for_cigar_interval_in_alignment(
            false,
            alignment,
            ref_intervals,
            query_intervals,
            [&](const CigarInterval& intersection, const interval_t& interval) {
                cerr << "R: " << intersection.ref_start << " " << intersection.ref_stop << " " << cigar_code_to_char[intersection.code] << " " << intersection.get_op_length() << '\n';

                // use the cigars to construct a formatted alignment
                get_formatted_sequence_of_cigar_interval(
                        intersection,
                        ref_seq,
                        query_seq,
                        a,
                        b,
                        x
                );

                a_all += a;
                x_all += x;
                b_all += b;

            },
            [&](const CigarInterval& intersection, const interval_t& interval) {
                return;
            }
    );

    cerr << a_all << '\n';
    cerr << x_all << '\n';
    cerr << b_all << '\n';

    string a_all_expected = "AATGGCAAACTTGCTCCCTAGCACTCATGTCCAGTGTAGACAAACCTGAAGAGGGTGTGTCTCCTC------------------------------------------TGACCCTCTCAAAGTACACAGGTGAAAAGCAGGAGTGAGAATGTTTGGGTTCTGAGCTCACACATGGACCAAGCTTCCCTAGGCTGGCTTGGGGCTATGTCTGACTTCTGTTCTTAAAAAGCCAATTCATGGGCCGGGCGCGGTGGCTCACACCTGCAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAGATAGAGACCATCCTGGCTAACACGGTGAAACCCCGTCTCTACTAAAAATACAAAAAATTAGCTTGGCGAGGTGGCGGGCGCCTGTAGTCCCAGCTACTTGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCCAGGGGGCGGAGCCTGCAGTGAGCCGAGATCACGCCACTGCACTCCAGCCTGGGCGACAGCGAGACTCCGTCTCAAAAAACAAACAAACAAACAAACAAACAAAAAGCCAATTCATGCTGGGTGCAGTGGCTCACACCTATAATCCCAGCACTTTGGGAGGCTGAGGTGGGAGGATTGCTTTGAGCTCAGGAGAGAGCAGCCTGGGCAACAGGGCAAAAACCCATCTCTACAAAAAAAAAATACAAAAAATTAGCCAGGTATGGTGGCGTGGTGGCTCACATCTGTGGTTCCAGCTACTCAAGAGGCTGAGGCAGGGGAATTACTTGAACCCAGGAGGTGGAGGTTGCCGTGAGCCGAGATTGCACCACTGTACTTCAGCCTGGGTGACAGAGTGAGACCCTGTCTCAGAAAAAAAAAAAAAAAAAAAAAAAAAAA-GGCCAATTCAACAGCCATCTCCTGCAAGTTCTGTGATGCC";
    string x_all_expected = "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||                                          ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ||||||||||||||||||||||||||||||||||||||||";
    string b_all_expected = "AATGGCAAACTTGCTCCCTAGCACTCATGTCCAGTGTAGACAAACCTGAAGAGGGTGTGTCTCCTCCTGACCCAGTGTAGACAAAGCTGAAGAGGGTGTGTCTCCTCCTGACCCTCTCAAAGTACACAGGTGAAAAGCAGGAGTGAGAATGTTTGGGTTCTGAGCTCACACATGGACCAAGCTTCCCTAGGCTGGCTTGGGGCTATGTCTGACTTCTGTTCTTAAAAAGCCAATTCATGGGCCGGGCGCGGTGGCTCACACCTGCAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAGATAGAGACCATCCTGGCTAACACGGTGAAACCCCGTCTCTACTAAAAATACAAAAAATTAGCTTGGCGAGGTGGCGGGCGCCTGTAGTCCCAGCTACTTGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCCAGGGG-CGGAGCCTGCAGTGAGCCGAGATCACGCCACTGCACTCCAGCCTGGGCGACAGCGAGACTCCGTCTCAAAAAACAAACAAACAAACAAACAAACAAAAAGCCAATTCATGCTGGGTGCAGTGGCTCACACCTATAATCCCAGCACTTTGGGAGGCTGAGGTGGGAGGATTGCTTTGAGCTCAGGAGAGAGCAGCCTGGGCAACAGGGCAAAA-CCCATCTCTACAAAAAAAAAATACAAAAAATTAGCCAGGTATGGTGGCGTGGTGGCTCACATCTGTGGTTCCAGCTACTCAAGAGGCTGAGGCAGGGGAATTACTTGAACCCAGGAGGTGGAGGTTGCCGTGAGCCGAGATTGCACCACTGTACTTCAGCCTGGGTGACAGAGTGAGACCCTGTCTCAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGCCAATTCAACAGCCATCTCCTGCAAGTTCTGTGATGCC";


    if (a_all != a_all_expected) {
        throw runtime_error("ERROR: a_all != a_all_expected");
    }

    if (x_all != x_all_expected) {
        throw runtime_error("ERROR: x_all != x_all_expected");
    }

    if (b_all != b_all_expected) {
        throw runtime_error("ERROR: b_all != b_all_expected");
    }

    cout << "PASS: SimpleAlignment test" << std::endl;
}


int main() {
    try {
        test_simple_alignment();
    } catch (const std::exception& e) {
        cerr << e.what() << std::endl;
        return 1;
    }
    return 0;
}
