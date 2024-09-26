#include <iostream>
#include <stdexcept>
#include <vector>

#include "HalfInterval.hpp"

using std::runtime_error;
using std::cerr;

using sv_merge::interval_t;
using sv_merge::deoverlap_intervals;

bool test_containment(){
    bool pass = true;

    // construct a vector of intervals in which one is contained by the other
    vector<interval_t> intervals = {
            {0, 10},
            {2, 8}
    };

    // construct a vector to store the deoverlapped intervals
    vector<interval_t> result_intervals;

    // construct a map to store the mapping of the deoverlapped intervals
    unordered_map<interval_t, size_t> result_mapping;

    // deoverlap the intervals
    deoverlap_intervals(intervals, result_intervals, result_mapping);

    // print the deoverlapped intervals
    for (const auto& interval : result_intervals){
        cerr << interval.first << ',' << interval.second << '\n';
    }

    // construct the expected result
    vector<interval_t> expected_result = {
            {0, 2},
            {2, 8},
            {8, 10}
    };

    // check if the result is as expected
    if (result_intervals != expected_result){
        pass = false;
    }

    return pass;
}


bool test_overlapping(){
    bool pass = true;

    // construct a vector of intervals in which one is contained by the other
    vector<interval_t> intervals = {
            {0, 10},
            {5, 15}
    };

    // construct a vector to store the deoverlapped intervals
    vector<interval_t> result_intervals;

    // construct a map to store the mapping of the deoverlapped intervals
    unordered_map<interval_t, size_t> result_mapping;

    // deoverlap the intervals
    deoverlap_intervals(intervals, result_intervals, result_mapping);

    // print the deoverlapped intervals
    for (const auto& interval : result_intervals){
        cerr << interval.first << ',' << interval.second << '\n';
    }

    // construct the expected result
    vector<interval_t> expected_result = {
            {0, 5},
            {5, 10},
            {10, 15}
    };

    // check if the result is as expected
    if (result_intervals != expected_result){
        pass = false;
    }

    return pass;
}


bool test_singleton(){
    bool pass = true;

    // construct a trivial vector of only one interval
    vector<interval_t> intervals = {
            {0, 10}
    };

    // construct a vector to store the deoverlapped intervals
    vector<interval_t> result_intervals;

    // construct a map to store the mapping of the deoverlapped intervals
    unordered_map<interval_t, size_t> result_mapping;

    // deoverlap the intervals
    deoverlap_intervals(intervals, result_intervals, result_mapping);

    // print the deoverlapped intervals
    for (const auto& interval : result_intervals){
        cerr << interval.first << ',' << interval.second << '\n';
    }

    // construct the expected result
    vector<interval_t> expected_result = {
            {0, 10}
    };

    // check if the result is as expected
    if (result_intervals != expected_result){
        pass = false;
    }

    return pass;
}


int main(){

    cerr << "Running tests for HalfInterval\n";

    cerr << "test_containment" << '\n';
    if (not test_containment()){
        throw runtime_error("ERROR: test_containment failed");
    }

    cerr << "test_overlapping" << '\n';
    if (not test_overlapping()){
        throw runtime_error("ERROR: test_overlapping failed");
    }

    cerr << "test_singleton" << '\n';
    if (not test_singleton()){
        throw runtime_error("ERROR: test_singleton failed");
    }

    return 0;
}
