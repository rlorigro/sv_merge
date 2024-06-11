#pragma once

#include "Sequence.hpp"

#include <functional>
#include <utility>
#include <vector>
#include <string>
#include <array>

using std::function;
using std::string;
using std::vector;
using std::array;
using std::pair;

#include "misc.hpp"

namespace sv_merge {

/**
 * ****************************
 * *** CIGAR related macros ***
 * ****************************
 *
 * #define BAM_CMATCH      0
 * #define BAM_CINS        1
 * #define BAM_CDEL        2
 * #define BAM_CREF_SKIP   3
 * #define BAM_CSOFT_CLIP  4
 * #define BAM_CHARD_CLIP  5
 * #define BAM_CPAD        6
 * #define BAM_CEQUAL      7
 * #define BAM_CDIFF       8
 * #define BAM_CBACK       9
 *
**/


/**
 * Lookup table to fetch the char representation of a cigar code
 */
static const array <char,9> cigar_code_to_char = {
        'M', // 0
        'I', // 1
        'D', // 2
        'N', // 3
        'S', // 4
        'H', // 5
        'P', // 6
        '=', // 7
        'X', // 8
};


//                                                     0   1   2   3   4   5   6   7   8   9
static const array<uint8_t, 128> cigar_char_to_code = {10, 10, 10, 10, 10, 10, 10, 10, 10, 10,  // 0
                                                       10, 10, 10, 10, 10, 10, 10, 10, 10, 10,  // 10
                                                       10, 10, 10, 10, 10, 10, 10, 10, 10, 10,  // 20
                                                       10, 10, 10, 10, 10, 10, 10, 10, 10, 10,  // 30
                                                       10, 10, 10, 10, 10, 10, 10, 10, 10, 10,  // 40
                                                       10, 10, 10, 10, 10, 10, 10, 10, 10, 10,  // 50
                                                       10, 7,  10, 10, 10, 10, 10, 10, 2,  10,  // 60    =61, D68
                                                       10, 10, 5,  1,  10, 10, 10, 0,  3,  10,  // 70    H72, I73, M77, N78
                                                       6,  10, 10, 4,  10, 10, 10, 10, 8,  10,  // 80    P80, S83, X88
                                                       10, 10, 10, 10, 10, 10, 10, 10, 10, 10,  // 90
                                                       10, 10, 10, 10, 10, 10, 10, 10, 10, 10,  // 100
                                                       10, 10, 10, 10, 10, 10, 10, 10, 10, 10,  // 110
                                                       10, 10, 10, 10, 10, 10, 10, 10};         // 120


/**
 * Lookup table to fetch a formatted alignment character which indicates match/mismatch/gap
 */
static const array <char,9> cigar_code_to_format_char = {
        '|', // M
        ' ', // I
        ' ', // D
        ' ', // N
        ' ', // S
        ' ', // H
        ' ', // P
        '|', // =
        '*', // X
};


/**
 * Lookup table to check if a cigar code consumes query sequence
 */
static const array <bool,9> is_query_move = {
        1, // M
        1, // I
        0, // D
        0, // N
        1, // S
        0, // H
        0, // P
        1, // =
        1  // X
};


/**
 * Lookup table to check if a cigar code consumes reference sequence
 */
static const array <bool,9> is_ref_move = {
        1, // M
        0, // I
        1, // D
        1, // N
        0, // S
        0, // H
        0, // P
        1, // =
        1  // X
};


/**
 * Data class which represents a cigar tuple
 */
class CigarTuple{
public:
    int32_t length = -1;
    int8_t code = -1;

    CigarTuple(int64_t length, int8_t code);
    CigarTuple()=default;
};


/**
 * Data class which stores all the necessary cigar data and coordinates to map from a query interval to a ref interval
 */
class CigarInterval{
public:
    int32_t length = -1;
    int32_t ref_start = -1;
    int32_t ref_stop = -1;
    int32_t query_start = -1;
    int32_t query_stop = -1;
    int8_t code = -1;
    bool is_reverse;

    interval_t get_forward_ref_interval() const;
    interval_t get_forward_query_interval() const;
    int32_t get_ref_length() const;
    int32_t get_query_length() const;
    void set_ref_interval_forward();
    void set_query_interval_forward();
    void set_ref_interval_reverse();
    void set_query_interval_reverse();
    bool is_softclip() const;
    bool is_hardclip() const;
    bool is_clip() const;
};


/**
 * Abstract class (API) which defines the necessary functions for an Alignment.
 * To be overridden by BAM (htslib) and GAF implementations so that the cigar interval iterator can operate seamlessly
 * on both.
 */
class Alignment {
public:

    /**
     * This method allows the user to iterate a cigar string and fetch absolute intervals instead of a delta
     * The purpose is to factor out the boilerplate code associated with unrolling cigars.
     * @param f - lambda function to operate on each cigar interval
     */
    virtual void for_each_cigar_interval(bool unclip_coords, const function<void(const CigarInterval& cigar)>& f) = 0;


    /**
     * Very simple method to iterate the cigar tuples in a C++/OOP friendly way
     * @param f - lambda function to operate on each cigar tuple
     */
    virtual void for_each_cigar_tuple(const function<void(const CigarTuple& cigar)>& f) = 0;

    virtual int32_t get_query_length() const = 0;
    virtual void get_query_sequence(string& result, int32_t start, int32_t stop) = 0;
    virtual void get_query_sequence(string& result) = 0;
    virtual void get_qualities(vector<uint8_t>& result) = 0;
    virtual void get_query_name(string& result) const = 0;
    virtual void get_tag_as_string(const string& name, string& result, bool allow_missing=false) const = 0;
    virtual int32_t get_ref_start() const = 0;
    virtual int32_t get_ref_stop() const = 0;
    virtual int32_t get_query_start() const = 0;
    virtual bool is_unmapped() const = 0;
    virtual bool is_reverse() const = 0;
    virtual bool is_primary() const = 0;
    virtual bool is_supplementary() const = 0;
};


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
void for_cigar_interval_in_alignment(
        bool unclip_coords,
        Alignment& alignment,
        vector<interval_t>& ref_intervals,
        vector<interval_t>& query_intervals,
        const function<void(const CigarInterval& intersection, const interval_t& interval)>& f_ref,
        const function<void(const CigarInterval& intersection, const interval_t& interval)>& f_query
);


void get_formatted_sequence_of_cigar_interval(
        const CigarInterval& cigar,
        const string& ref_sequence,
        const string& query_sequence,
        string& s_ref,
        string& s_query,
        string& s_crossref
);

}
