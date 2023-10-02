#include <iostream>
#include <stdexcept>

using std::runtime_error;
using std::exception;
using std::cerr;


#include "IntervalGraph.hpp"
#include "pair_hash.hpp"

using sv_merge::IntervalNode;
using sv_merge::IntervalGraph;
using sv_merge::interval_t;


template <class T> void print_component(const IntervalGraph<T>& g, const unordered_set <interval_t>& component){
    for (const auto& interval: component){
        cerr << '\t' << interval.first << ',' << interval.second << ' ';

        g.for_value_in_interval(interval, [&](const T& value){
            cerr << value << ',';
        });

        cerr << '\n';
    }
}


template <class T> void print_components(const IntervalGraph<T>& g, const vector <unordered_set <interval_t> >& components){
    int c = 0;
    for (const auto& component: components){
        cerr << c++ << '\n';
        print_component(g, component);
        cerr << '\n';
    }
}


void test_standard(){
    vector <pair <interval_t,unordered_set<string> > > labeled_intervals = {
        {{-1,9}, {"a"}},    // 0
        {{0,10}, {"b"}},    // 0 collapse v
        {{0,10}, {"c"}},    // 0 collapse ^
        {{10,20}, {"d"}},   // 0
        {{11,20}, {"e"}},   // 0
        {{21,25}, {"f"}},   // 1
        {{22,26}, {"g"}}    // 1
    };

    IntervalGraph<string> g(labeled_intervals);

    vector <unordered_set <interval_t> > components;

    g.get_connected_components(components);

    vector <unordered_set <interval_t> > expected_result = {
        {
            {-1, 9},
            {0, 10},
            {10, 20},
            {11, 20}
        },
        {
            {21, 25},
            {22, 26},
        }
    };

    unordered_map <interval_t, unordered_set<string> > expected_values = {
            {{-1, 9},{"a"}},
            {{0, 10},{"b", "c"}},
            {{10, 20},{"d"}},
            {{11, 20},{"e"}},
            {{21, 25},{"f"}},
            {{22, 26},{"g"}}
    };

    int c = 0;
    for (const auto& component: components){
        cerr << "searching for matching component: " << c << '\n';
        print_component(g, component);

        bool pass = false;

        // Check that the components are structured as expected
        auto iter = expected_result.begin();
        while (iter != expected_result.end()){
            if (component == *iter){
                cerr << "found" << '\n';
                pass = true;
                break;
            }
            iter++;
        }

        expected_result.erase(iter);

        if (not pass){
            throw runtime_error("FAIL: components not expected");
        }

        // Check that all the values for each node in the component are as expected
        pass = true;
        for (const auto& interval: component){
            auto n = g.get_node(interval);

            auto result = expected_values.find(interval);

            if (result != expected_values.end()){
                if (result->second != n.values){
                    pass = false;
                }
            }
            else{
                pass = false;
            }
        }

        if (not pass){
            throw runtime_error("FAIL: values not expected");
        }

        c++;
    }

    if (not expected_result.empty()){
        throw runtime_error("FAIL: not all expected results are found in output");
    }

    cerr << "PASS: standard" << '\n';
}


void test_component_interval_iterator(){
    vector <pair <interval_t,unordered_set<string> > > labeled_intervals = {
        {{-1,9}, {"a"}},    // 0
        {{0,10}, {"b"}},    // 0 collapse v
        {{0,10}, {"c"}},    // 0 collapse ^
        {{10,20}, {"d"}},   // 0
        {{11,20}, {"e"}},   // 0
        {{21,25}, {"f"}},   // 1
        {{22,26}, {"g"}}    // 1
    };

    IntervalGraph<string> g(labeled_intervals);

    unordered_map <interval_t, unordered_set<string> > expected_results = {
        {{-1, 20},{"a","b","c","d","e"}},
        {{21, 26},{"f","g"}}
    };

    unordered_map <interval_t, unordered_set<string> > results;

    g.for_each_connected_component_interval([&](interval_t& interval, unordered_set<string>& values){
        cerr << interval.first << "," << interval.second;
        for (const auto& v: values){
            cerr << ' ' << v;
        }
        cerr << '\n';

        results[interval] = values;
    });

    for (const auto& [interval, values]: results){
        auto expected_result = expected_results.find(interval);

        if (expected_result == expected_results.end()){
            throw runtime_error("FAIL: not all results are found in expected results");
        }
        else{
            if (expected_result->second != values){
                throw runtime_error("FAIL: result for interval " + to_string(interval.first) + "," + to_string(interval.second) + " do not match expected results");
            }
        }

        expected_results.erase(expected_result);
    }

    if (not expected_results.empty()){
        throw runtime_error("FAIL: not all expected results are found in results");
    }

    cerr << "PASS: interval iterator" << '\n';
}


void test_bad_interval(){
    vector <pair <interval_t,unordered_set<string> > > labeled_intervals = {
        {{0,10}, {"a"}},
        {{0,10}, {"b"}},
        {{0,10}, {"c"}},
        {{10,20}, {"d"}},
        {{11,20}, {"e"}},
        {{11,0}, {"f"}}
    };

    try {
        IntervalGraph<string> g(labeled_intervals);
    }
    catch (exception& e){
        auto expected_error = "ERROR: interval start is greater than interval stop";
        string error = e.what();


        if (error.find(expected_error) != std::string::npos){
            cerr << "PASS: bad_interval" << '\n';
        }
        else{
            cerr << e.what() << '\n';
            throw runtime_error("FAIL: unexpected error");
        }
    }
}


int main(){
    test_standard();
    test_bad_interval();
    test_component_interval_iterator();

    return 0;
}
