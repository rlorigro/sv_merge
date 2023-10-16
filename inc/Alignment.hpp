#pragma once

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


class CigarTuple{
public:
    int64_t length = -1;
    int8_t code = -1;
};


class CigarInterval{
public:
    int64_t length = -1;
    int64_t ref_start = -1;
    int64_t ref_stop = -1;
    int64_t query_start = -1;
    int64_t query_stop = -1;
    int8_t code = -1;
    bool is_reverse;

    pair<int64_t,int64_t> get_forward_ref_interval() const;
    pair<int64_t,int64_t> get_forward_query_interval() const;
    int64_t get_ref_length() const;
    int64_t get_query_length() const;
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
    virtual void for_each_cigar_interval(const function<void(const CigarInterval& cigar)>& f) = 0;
    virtual void for_each_cigar_tuple(const function<void(const CigarTuple& cigar)>& f) = 0;
    virtual int64_t get_query_length() const = 0;
    virtual void get_query_sequence(string& result) = 0;
    virtual void get_query_name(string& result) const = 0;
    virtual int64_t get_ref_start() const = 0;
    virtual int64_t get_query_start() const = 0;
    virtual bool is_unmapped() const = 0;
    virtual bool is_reverse() const = 0;
};


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
        Alignment& alignment,
        vector<interval_t>& ref_intervals,
        vector<interval_t>& query_intervals,
        const function<void(const CigarInterval& intersection, const interval_t& interval)>& f_ref,
        const function<void(const CigarInterval& intersection, const interval_t& interval)>& f_query
);




}