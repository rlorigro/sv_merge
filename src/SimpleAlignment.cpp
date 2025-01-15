#include "SimpleAlignment.hpp"
#include "IntervalGraph.hpp"
#include "Alignment.hpp"
#include "Sequence.hpp"
#include "Region.hpp"

#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <ranges>
#include <deque>
#include <span>

using std::priority_queue;
using std::runtime_error;
using std::to_string;
using std::ofstream;
using std::deque;
using std::swap;
using std::cerr;
using std::cout;
using std::span;
using std::min;
using std::max;

using std::runtime_error;
using std::cerr;

namespace sv_merge{


SimpleAlignment::SimpleAlignment(const string& ref_sequence, const string& query_sequence, const string& cigar):
        ref_sequence(ref_sequence),
        query_sequence(query_sequence),
        cigar(cigar)
{
    // Check that the cigar string is not empty
    if (cigar.empty()){
        throw runtime_error("ERROR: CIGAR string is empty");
    }
}


void SimpleAlignment::for_each_cigar_tuple(const function<void(const CigarTuple& cigar)>& f){
    CigarTuple cigar_tuple;
    string length_token;

    for (auto c: cigar){
        if (isalpha(c) or c == '='){
            cigar_tuple.length = stol(length_token);
            cigar_tuple.code = cigar_char_to_code[c];

//            cerr << "----\n";
//            cerr << "cigar_tuple.length: " << cigar_tuple.length << '\n';
//            cerr << "cigar_tuple.code: " << c << '\n';

            f(cigar_tuple);

            length_token.clear();
        }
        else{
            length_token += c;
        }
    }
}


void SimpleAlignment::for_each_cigar_interval(bool unclip_coords, const function<void(const CigarInterval& cigar)>& f){
    CigarInterval c;

    // SimpleAlignment does not support reverse alignments
    c.is_reverse = false;

    // Initialize the cigar interval
    c.query_start = get_query_start();
    c.ref_start = get_ref_start();

    for_each_cigar_tuple([&](const CigarTuple& tuple){
        c.code = tuple.code;
        c.length = tuple.length;

        c.ref_stop = c.ref_start + is_ref_move[c.code]*c.length;
        c.query_stop = c.query_start + is_query_move[c.code]*c.length;

//        cerr << "----\n";
//        cerr << "c.ref_start: " << c.ref_start << '\n';
//        cerr << "c.ref_stop: " << c.ref_stop << '\n';
//        cerr << "c.query_start: " << c.query_start << '\n';
//        cerr << "c.query_stop: " << c.query_stop << '\n';

        f(c);

        c.ref_start = c.ref_stop;
        c.query_start = c.query_stop;
    });
}


void SimpleAlignment::get_query_sequence(string& result, int32_t start, int32_t stop) {
    result = query_sequence.substr(start, stop-start);
}

void SimpleAlignment::get_query_sequence(BinarySequence<uint64_t>& result) {
    throw runtime_error("ERROR: get_qualities not implemented for SimpleAlignment with BinarySequence");
}

void SimpleAlignment::get_query_sequence(string& result) {
    result = query_sequence;
}

void SimpleAlignment::get_qualities(vector<uint8_t>& result) {
    throw runtime_error("ERROR: get_qualities not implemented for SimpleAlignment");
}

void SimpleAlignment::get_query_name(string& result) const {
    throw runtime_error("ERROR: get_query_name not implemented for SimpleAlignment");
}

void SimpleAlignment::get_tag_as_string(const string& name, string& result, bool allow_missing) const {
    throw runtime_error("ERROR: get_tag_as_string not implemented for SimpleAlignment");
}

[[nodiscard]] int32_t SimpleAlignment::get_query_length() const {
    return query_sequence.size();
}

[[nodiscard]] int32_t SimpleAlignment::get_ref_start() const {
    return 0;
}

[[nodiscard]] int32_t SimpleAlignment::get_ref_stop() const {
    return ref_sequence.size();
}

[[nodiscard]] int32_t SimpleAlignment::get_query_start() const {
    return 0;
}

[[nodiscard]] bool SimpleAlignment::is_unmapped() const {
    return false;
}

[[nodiscard]] bool SimpleAlignment::is_reverse() const {
    return false;
}

[[nodiscard]] bool SimpleAlignment::is_primary() const {
    throw runtime_error("ERROR: is_primary not implemented for SimpleAlignment");
}

[[nodiscard]] bool SimpleAlignment::is_supplementary() const {
    throw runtime_error("ERROR: is_supplementary not implemented for SimpleAlignment");
}



}