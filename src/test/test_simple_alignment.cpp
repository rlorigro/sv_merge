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
    // ref_seq:   ACGTACGTACGT
    // query_seq: ACGT-TG-ACGT
    string ref_seq = "ACGTACGTACGT";
    string query_seq = "ACGTTGACGT";

    // Use WFA2 to generate the CIGAR string
    WFAlignerGapAffine aligner(4,2,2,WFAligner::Alignment,WFAligner::MemoryHigh);

    aligner.alignEnd2End(ref_seq, query_seq);

    // Extract indel edit distance from cigar
    string cigar = aligner.getCIGAR(false);

    // Construct SimpleAlignment
    SimpleAlignment alignment(ref_seq, query_seq, cigar);

    // Define ref intervals to be iterated
    std::vector<interval_t> ref_intervals = {{0, 2}, {3, 8}};
    std::vector<interval_t> query_intervals = {};

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
                cerr << "Ref interval: " << interval.first << " " << interval.second << '\n';
                cerr << "Intersection: " << cigar_code_to_char[intersection.code] << " " << intersection.length << '\n';
            },
            [&](const CigarInterval& intersection, const interval_t& interval) {
                return;
            }
    );

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
