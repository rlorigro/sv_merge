#include "Alignment.hpp"

#include <stdexcept>
#include <iostream>

using std::runtime_error;
using std::cerr;
using std::swap;
using std::cerr;
using std::cout;
using std::min;
using std::max;


namespace sv_merge{


CigarTuple::CigarTuple(int64_t length, int8_t code):
    length(length),
    code(code)
{}


/**
 * Regardles of the F/R orientation of the alignment, return an interval in which the second item is greater than the
 * first
 * @return a pair of integers representing the interval
 */
interval_t CigarInterval::get_forward_ref_interval() const{
    if (ref_start > ref_stop){
        return {ref_stop,ref_start};
    }
    else {
        return {ref_start,ref_stop};
    }
}


int32_t CigarInterval::get_ref_length() const{
    return abs(ref_stop-ref_start);
}


int32_t CigarInterval::get_query_length() const{
    return abs(query_stop-query_start);
}


int32_t CigarInterval::get_op_length() const{
    // Return the maximum of the two lengths, since the cigar interval could be a deletion or insertion
    return max(abs(ref_stop-ref_start), abs(query_stop-query_start));
}



/**
 * Similar to the get_forward_..._interval method but this updates the Cigar object as [a,b] where a < b
 */
void CigarInterval::set_ref_interval_forward(){
    if (ref_start > ref_stop){
        swap(ref_stop,ref_start);
    }
}



/**
 * Similar to the get_reverse_..._interval method but this updates the Cigar object as [b,a] where a < b
 */
void CigarInterval::set_ref_interval_reverse(){
    if (ref_start < ref_stop){
        swap(ref_stop,ref_start);
    }
}


/**
 * Regardles of the F/R orientation of the alignment, return an interval in which the second item is greater than the
 * first
 * @return a pair of integers representing the interval
 */
interval_t CigarInterval::get_forward_query_interval() const{
    if (query_start > query_stop){
        return {query_stop,query_start};
    }
    else {
        return {query_start,query_stop};
    }
}


/**
 * Similar to the get_forward_..._interval method but this updates the Cigar object as [a,b] where a < b
 */
void CigarInterval::set_query_interval_forward(){
    if (query_start > query_stop){
        swap(query_stop,query_start);
    }
}


/**
 * Similar to the get_reverse_..._interval method but this updates the Cigar object as [b,a] where a < b
 */
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
void for_cigar_interval_in_alignment(
        bool unclip_coords,
        Alignment& alignment,
        vector<interval_t>& ref_intervals,
        vector<interval_t>& query_intervals,
        const function<void(const CigarInterval& intersection, const interval_t& interval)>& f_ref,
        const function<void(const CigarInterval& intersection, const interval_t& interval)>& f_query
){

    // Comparator to sort intervals [(a,b),...,(a,b)] by start, 'a'
    // Since intervals are expected to be non overlapping, choosing left/right start/end doesn't matter
    auto left_comparator = [](const interval_t& a, const interval_t& b){
        return a.first < b.first;
    };

    auto left_comparator_reverse = [](const interval_t& a, const interval_t& b){
        return a.first > b.first;
    };

    CigarInterval intersection;

    // Sort the query intervals
    sort(query_intervals.begin(), query_intervals.end(), alignment.is_reverse() ? left_comparator_reverse : left_comparator);

    // Sort the ref intervals
    sort(ref_intervals.begin(), ref_intervals.end(), left_comparator);

    // Setup iterators which walk along the ref/query intervals (if they exist)
    auto ref_iter = ref_intervals.begin();
    auto query_iter = query_intervals.begin();

    // Make sure to transfer the reversal status
    intersection.is_reverse = alignment.is_reverse();

    alignment.for_each_cigar_interval(unclip_coords, [&](const sv_merge::CigarInterval& c) {
//        cerr << "(" << cigar_code_to_char[c.code] << "," << c.length << ") " << c.is_reverse << " -- r:" << c.ref_start << ',' << c.ref_stop << " q:" << c.query_start << ',' << c.query_stop << '\n';

        while (ref_iter != ref_intervals.end()) {
            intersection.code = c.code;

            intersection.ref_start = max(ref_iter->first, c.ref_start);
            intersection.ref_stop = min(ref_iter->second, c.ref_stop);

            auto l = intersection.ref_stop - intersection.ref_start;
            intersection.length = l;

//            cerr << "++ r:" << ref_iter->first << ',' << ref_iter->second << " q:" << c.query_start << ',' << c.query_stop << '\n';

            // If this is an M/=/X operation, then the fates of the ref/query intervals are tied
            if (is_ref_move[c.code] and is_query_move[c.code]) {
                if (c.is_reverse) {
                    intersection.query_start = c.query_stop - (intersection.ref_start - c.ref_start);
                    intersection.query_stop = c.query_stop - (intersection.ref_stop - c.ref_start);
                }
                else {
                    intersection.query_start = c.query_start + (intersection.ref_start - c.ref_start);
                    intersection.query_stop = c.query_start + (intersection.ref_stop - c.ref_start);
                }
            }
            // Otherwise the ref/query advance independently
            else {
                intersection.query_start = c.query_start;
                intersection.query_stop = c.query_stop;
            }

            // In some cases, the window could be completely non overlapping
            if (l >= 0) {
                f_ref(intersection, *ref_iter);
            }

            // Cases A and B indicate that more windows need to be consumed before advancing the cigar interval.
            // Otherwise, stop iterating
            if (c.ref_stop <= ref_iter->second) {
                // If the window is exactly caught up with the cigar interval, then advance to the next window
                if (c.ref_stop == ref_iter->second) {
                    ref_iter++;
                }

                // Otherwise, need to keep this window open for the next cigar
                break;
            }

            ref_iter++;
        }


        // del_5_at_34_reverse
        // Cigars in ref coords:
        //        =34       [----)
        //        D5             [--)
        //        =1951             [----------------------)

        // window [0,46)                                (--]
        // cigar  [46,100)                          (---]
        // cigar  [120,140)                   (--]

        // del_5_at_34
        // Cigars in ref coords:
        //        =34       [----)
        //        D5             [--)
        //        =1951             [----------------------)

        // window [0,46)    [---------)
        // cigar  [46,100)            [-----)
        // cigar  [120,140)                    [--)

        while (query_iter != query_intervals.end()) {
            intersection.code = c.code;

            intersection.query_start = max(query_iter->first, c.query_start);
            intersection.query_stop = min(query_iter->second, c.query_stop);

            auto l = intersection.query_stop - intersection.query_start;
            intersection.length = l;

//            cerr << "++ r:" << c.ref_start << ',' << c.ref_stop << " q:" << query_iter->first << ',' << query_iter->second << '\n';

            // If this is an M/=/X operation, then the fates of the ref/query intervals are tied
            if (is_ref_move[c.code] and is_query_move[c.code]) {
                if (c.is_reverse) {
                    intersection.ref_start = c.ref_stop - (intersection.query_start - c.query_start);
                    intersection.ref_stop = c.ref_stop - (intersection.query_stop - c.query_start);
                }
                else {
                    intersection.ref_start = c.ref_start + (intersection.query_start - c.query_start);
                    intersection.ref_stop = c.ref_start + (intersection.query_stop - c.query_start);
                }
            }
                // Otherwise the ref/query advance independently
            else {
                intersection.ref_start = c.ref_start;
                intersection.ref_stop = c.ref_stop;
            }

            // In some cases, the window could be completely non overlapping
            if (l >= 0) {
                f_query(intersection, *query_iter);
            }

            bool window_exceeds_cigar;
            if (c.is_reverse) {
                window_exceeds_cigar = query_iter->first <= c.query_start;
            }
            else {
                window_exceeds_cigar = c.query_stop <= query_iter->second;
            }

            // Cases A and B indicate that more windows need to be consumed before advancing the cigar interval.
            // Otherwise, stop iterating
            if (window_exceeds_cigar) {

                // If the window is exactly caught up with the cigar interval, then advance to the next window
                if (c.is_reverse) {
                    if (query_iter->first == c.query_start) {
                        query_iter++;
                    }
                }
                else {
                    if (c.query_stop == query_iter->second) {
                        query_iter++;
                    }
                }

                // Otherwise, need to keep this window open for the next cigar
                break;
            }

            query_iter++;
        }

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

//        cerr << l_query << ',' << l_ref << ',' << a_query << ',' << a_query + l_ref << '\n';

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


}
