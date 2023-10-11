#include "IntervalGraph.hpp"
#include "bam.hpp"

#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <ranges>
#include <deque>
#include <cmath>

using std::runtime_error;
using std::to_string;
using std::deque;
using std::swap;
using std::cerr;
using std::cout;
using std::min;
using std::max;


#include "htslib/include/htslib/hts.h"


namespace sv_merge {

/**
 * Inefficient complementation, good luck to the optimizer TODO: switch to array based or char ordinal offset math
 * @param c - the character to be complemented, must be ACGTN or acgtn
 * @return The complement of c
 */
char get_complement(char c){
    if (isupper(c)) {
        switch (c){
            case 'A': return 'T';
            case 'C': return 'G';
            case 'G': return 'C';
            case 'T': return 'A';
            case 'N': return 'N';
            default:
                throw runtime_error("ERROR: cannot complement base: " + string(1,c));
        }
    }

    else {
        switch (c){
            case 'a': return 't';
            case 'c': return 'g';
            case 'g': return 'c';
            case 't': return 'a';
            case 'n': return 'n';
            default:
                throw runtime_error("ERROR: cannot complement base: " + string(1,c));
        }
    }
}


void reverse_complement(string& seq){
    string s;

    for (char & iter : std::ranges::reverse_view(seq)){
        s += get_complement(iter);
    }

    seq = s;
}


pair<int64_t,int64_t> CigarInterval::get_forward_ref_interval() const{
    if (ref_start > ref_stop){
        return {ref_stop,ref_start};
    }
    else {
        return {ref_start,ref_stop};
    }
}


void CigarInterval::set_ref_interval_forward(){
    if (ref_start > ref_stop){
        swap(ref_stop,ref_start);
    }
}


void CigarInterval::set_ref_interval_reverse(){
    if (ref_start < ref_stop){
        swap(ref_stop,ref_start);
    }
}


pair<int64_t,int64_t> CigarInterval::get_forward_query_interval() const{
    if (query_start > query_stop){
        return {query_stop,query_start};
    }
    else {
        return {query_start,query_stop};
    }
}


void CigarInterval::set_query_interval_forward(){
    if (query_start > query_stop){
        swap(query_stop,query_start);
    }
}


void CigarInterval::set_query_interval_reverse(){
    if (query_start < query_stop){
        swap(query_stop,query_start);
    }
}


bool CigarInterval::is_softclip() const{
    return code == 4;
}


bool CigarInterval::is_hardclip() const{
    return code == 5;
}


bool CigarInterval::is_clip() const{
    return code == 4 or code == 5;
}


void decompress_cigar_bytes(uint32_t bytes, CigarTuple& cigar){
    cigar.code = int8_t(bytes & bam_cigar_mask);
    cigar.length = int64_t(bytes >> bam_cigar_shift);
}


void decompress_cigar_bytes(uint32_t bytes, CigarInterval& cigar){
    cigar.code = int8_t(bytes & bam_cigar_mask);
    cigar.length = int64_t(bytes >> bam_cigar_shift);
}


void for_cigar_tuple_in_alignment(const bam1_t *alignment, const function<void(CigarTuple& cigar)>& f){
    auto cigar_bytes = bam_get_cigar(alignment);
    for (uint32_t i=0; i<alignment->core.n_cigar; i++){
        CigarTuple c;
        decompress_cigar_bytes(cigar_bytes[i], c);
        f(c);
    }
}


void for_cigar_interval_in_alignment(const bam1_t *alignment, const function<void(CigarInterval& cigar)>& f){
    auto cigar_bytes = bam_get_cigar(alignment);

    CigarInterval c;
    c.query_start = 0;
    c.ref_start = alignment->core.pos;
    c.is_reverse = bam_is_rev(alignment);

    if (c.is_reverse){
        c.query_start = alignment->core.l_qseq;
    }

    for (uint32_t i=0; i<alignment->core.n_cigar; i++){
        decompress_cigar_bytes(cigar_bytes[i], c);

        if (c.is_reverse) {
            c.query_stop = c.query_start - is_query_move[c.code]*c.length;
        }
        else{
            c.query_stop = c.query_start + is_query_move[c.code]*c.length;
        }

        c.ref_stop = c.ref_start + is_ref_move[c.code]*c.length;

        f(c);

        c.ref_start = c.ref_stop;
        c.query_start = c.query_stop;
    }
}


// Cases A and B indicate that more windows need to be consumed before advancing the cigar interval
// FORWARD REF ALIGNMENT
// CASE A:
//         0123456789
// Cigar      [--+--)
// Window  [--+--)
//
// CASE B:
//         0123456
// Cigar   [--+--)
// Window   [-+-)
//
// CASE B2:
//         0123456
// Cigar   [--+--)
// Window  [--+-)
// CASE C:
// Cigar    [-+-)
// Window  [--+--)
//
// CASE C2:
// Cigar   [--+--)
// Window     [--+--)
//
// CASE C3:
// Cigar      [--)
// Window     [--+--)
//
// CASE D:
// Cigar    [-+--)
// Window  [--+--)
//
// CASE D2:
// Cigar   [--+-----)
// Window     [-----)
//
// CASE D3:
// Cigar      [-----)
// Window     [-----)

/**
 * Iterate cigar intervals and return cigar intervals that are clipped to match the bounds of each window provided by
 * the user as a vector of interval_t. Useful for parsing cigars over a region of interest.
 * @param alignment - pointer to htslib alignment struct
 * @param ref_intervals - intervals which MUST BE NON-OVERLAPPING and NO CHECK is performed to verify this! Intervals
 * must be half-open in F orientation, e.g.: [[0,2),[2,4)] are two adjacent intervals of length 2.
 * @param query_intervals - intervals which MUST BE NON-OVERLAPPING and NO CHECK is performed to verify this! Intervals
 * must be half-open in F orientation, e.g.: [[0,2),[2,4)] are two adjacent intervals of length 2.
 * @param f_ref lambda function to operate on ref intervals
 * @param f_query lambda function to operate on query intervals
 */
void for_cigar_interval_in_alignment(
        const bam1_t *alignment,
        vector<interval_t>& ref_intervals,
        vector<interval_t>& query_intervals,
        const function<void(const CigarInterval& intersection, const CigarInterval& original, const interval_t& interval)>& f_ref,
        const function<void(const CigarInterval& intersection, const CigarInterval& original, const interval_t& interval)>& f_query
        ){

    // Comparator to sort intervals [(a,b),...,(a,b)] by start, 'a'
    // Since intervals are expected to be non overlapping, choosing left/right start/end doesn't matter
    auto left_comparator = [](const interval_t& a, const interval_t& b){
        return a.first < b.first;
    };

//    auto left_comparator_reverse = [](const interval_t& a, const interval_t& b){
//        return a.first > b.first;
//    };

    auto cigar_bytes = bam_get_cigar(alignment);

    CigarInterval intersection;
    CigarInterval c;
    c.query_start = 0;
    c.ref_start = alignment->core.pos;
    c.is_reverse = bam_is_rev(alignment);

    // Sort the query intervals
    sort(query_intervals.begin(), query_intervals.end(), left_comparator);

    // Sort the ref intervals
    sort(ref_intervals.begin(), ref_intervals.end(), left_comparator);

    // Setup iterators which walk along the ref/query intervals (if they exist)
    auto ref_iter = ref_intervals.begin();
    auto query_iter = query_intervals.begin();

    if (c.is_reverse){
        c.query_start = alignment->core.l_qseq;
    }

    for (uint32_t i=0; i<alignment->core.n_cigar; i++){
        decompress_cigar_bytes(cigar_bytes[i], c);

        // Update interval bounds for this cigar interval
        if (c.is_reverse) {
            c.query_stop = c.query_start - is_query_move[c.code]*c.length;
            c.set_query_interval_forward();
        }
        else{
            c.query_stop = c.query_start + is_query_move[c.code]*c.length;
        }

        c.ref_stop = c.ref_start + is_ref_move[c.code]*c.length;

        cerr << "-- r:" << c.ref_start << ',' << c.ref_stop << " q:" << c.query_start << ',' << c.query_stop << '\n';

        while (ref_iter != ref_intervals.end()){
            intersection.code = c.code;
            intersection.length = c.length;

            intersection.ref_start = max(ref_iter->first, c.ref_start);
            intersection.ref_stop = min(ref_iter->second, c.ref_stop);

            auto l = intersection.ref_stop - intersection.ref_start;

            // If this is an M/=/X operation, then the fates of the ref/query intervals are tied
            if (is_ref_move[c.code] and is_query_move[c.code]){
                if (c.is_reverse) {
                    intersection.query_start = c.query_stop - l;
                    intersection.query_stop = c.query_stop;
                }
                else{
                    intersection.query_stop = c.query_start + l;
                    intersection.query_start = c.query_start;
                }
            }
            // Otherwise the ref/query advance independently
            else{
                intersection.query_start = c.query_start;
                intersection.query_stop = c.query_stop;
            }

            // In some cases, the window could be completely non overlapping
            if (l >= 0){
                f_ref(intersection, c, *ref_iter);
            }

            // Cases A and B indicate that more windows need to be consumed before advancing the cigar interval.
            // Otherwise, stop iterating
            if (c.ref_stop <= ref_iter->second){
                // If the window is exactly caught up with the cigar interval, then advance to the next window
                if (c.ref_stop == ref_iter->second){
                    ref_iter++;
                }

                // Otherwise, need to keep this window open for the next cigar
                break;
            }

            ref_iter++;
        }

        while (query_iter != query_intervals.end()){
            intersection.code = c.code;
            intersection.length = c.length;

            intersection.query_start = max(query_iter->first, c.query_start);
            intersection.query_stop = min(query_iter->second, c.query_stop);

            auto l = intersection.query_stop - intersection.query_start;

            // If this is an M/=/X operation, then the fates of the ref/query intervals are tied
            if (is_ref_move[c.code] and is_query_move[c.code]){
                intersection.ref_stop = c.ref_start + l;
                intersection.ref_start = c.ref_start;
            }
            // Otherwise the ref/query advance independently
            else{
                intersection.ref_stop = c.ref_stop;
                intersection.ref_start = c.ref_start;
            }

            // In some cases, the window could be completely non overlapping
            if (l >= 0){
                f_query(intersection, c, *query_iter);
            }

            // Cases A and B indicate that more windows need to be consumed before advancing the cigar interval.
            // Otherwise, stop iterating
            if (c.query_stop <= query_iter->second){

                // If the window is exactly caught up with the cigar interval, then advance to the next window
                if (c.query_stop == query_iter->second){
                    query_iter++;
                }

                // Otherwise, need to keep this window open for the next cigar
                break;
            }

            query_iter++;
        }

        if (c.is_reverse){
            c.set_query_interval_reverse();
        }

        c.ref_start = c.ref_stop;
        c.query_start = c.query_stop;
    }
}


void decompress_bam_sequence(const bam1_t* alignment, string& sequence){
    ///
    /// Convert the compressed representation of an aligned sequence into a string.
    /// Does NOT reverse complement the sequence
    ///
    auto length = alignment->core.l_qseq;
    auto compressed_sequence = bam_get_seq(alignment);
    bool is_reverse = bam_is_rev(alignment);

    sequence.clear();

    uint8_t base_code;
    string base;

    // Fetch 4 bit base code from the correct 8 bit integer and convert to a char
    int64_t start = 0;
    int64_t stop = length;
    int64_t increment = 1;

    if (is_reverse){
        start = length - 1;
        stop = -1;
        increment = -1;
    }

    for (int64_t i=start; i!=stop; i+=increment){
        uint64_t index = i/2;

        if (i%2 == 0){
            // Perform bit SHIFT and decode using the standard or complemented base map
            base_code = compressed_sequence[index] >> bam_sequence_shift;
            sequence += bases[is_reverse][base_code];
        }
        else {
            // Perform bit MASK and decode using the standard or complemented base map
            base_code = compressed_sequence[index] & bam_sequence_mask;
            sequence += bases[is_reverse][base_code];
        }
    }
}


void for_read_in_bam_region(path bam_path, string region, const function<void(Sequence& sequence)>& f) {
    for_alignment_in_bam_region(bam_path, region, [&](const bam1_t *alignment){
        if (alignment->core.tid < 0) {
            return;
        }

        Sequence s;
        s.name = bam_get_qname(alignment);

        decompress_bam_sequence(alignment, s.sequence);

        f(s);
    });
}


void for_alignment_in_bam_region(path bam_path, string region, const function<void(const bam1_t *alignment)>& f) {
    samFile *bam_file;
    bam_hdr_t *bam_header;
    hts_idx_t *bam_index;
    bam1_t *alignment;

    bam_file = nullptr;
    bam_index = nullptr;
    alignment = bam_init1();

    if ((bam_file = hts_open(bam_path.string().c_str(), "r")) == nullptr) {
        throw runtime_error("ERROR: Cannot open bam file: " + bam_path.string());
    }

    // bam index
    if ((bam_index = sam_index_load(bam_file, bam_path.string().c_str())) == nullptr) {
        throw runtime_error("ERROR: Cannot open index for bam file: " + bam_path.string() + "\n");
    }

    // bam header
    if ((bam_header = sam_hdr_read(bam_file)) == nullptr) {
        throw runtime_error("ERROR: Cannot open header for bam file: " + bam_path.string() + "\n");
    }

    hts_itr_t *itr = sam_itr_querys(bam_index, bam_header, region.c_str());

    while (sam_itr_next(bam_file, itr, alignment) >= 0) {
        if (alignment->core.tid < 0) {
            continue;
        }

        f(alignment);
    }

    bam_destroy1(alignment);
    hts_itr_destroy(itr);
    bam_hdr_destroy(bam_header);
    hts_close(bam_file);
    hts_idx_destroy(bam_index);
}


}
