#pragma once

#include "pair_hash.hpp"

#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <stdexcept>
#include <algorithm>
#include <ostream>
#include <utility>
#include <string>
#include <vector>
#include <memory>
#include <queue>
#include <set>

using std::unordered_map;
using std::unordered_set;
using std::runtime_error;
using std::to_string;
using std::set_union;
using std::function;
using std::ostream;
using std::vector;
using std::string;
using std::queue;
using std::sort;
using std::pair;
using std::set;


namespace sv_merge {

using interval_t = pair<int64_t,int64_t>;


template<class T> class IntervalNode{
public:
    unordered_set<T> values;
    unordered_set<interval_t> neighbors;

    IntervalNode()=default;
    explicit IntervalNode<T>(const unordered_set<T>& value);
    void update_values(const unordered_set<T>& value);
};


template <class T> void IntervalNode<T>::update_values(const unordered_set<T>& value){
    for(const auto& v: value){
        values.emplace(v);
    }
}


template <class T> IntervalNode<T>::IntervalNode(const unordered_set<T>& value):
    values(value),
    neighbors()
{}


template <class T> class IntervalGraph {
    unordered_map <interval_t, IntervalNode<T> > nodes;

public:
    explicit IntervalGraph<T>(vector <pair <interval_t,unordered_set<T> > >& intervals);
    void get_connected_components(vector <unordered_set<interval_t> >& components);
    void for_value_in_interval(const interval_t& interval, const function<void(const T& value)>& f) const;
    IntervalNode<T> get_node(const interval_t& interval);
};


template <class T> IntervalGraph<T>::IntervalGraph(vector <pair <interval_t,unordered_set<T> > >& labeled_intervals){
    // How to sort labeled intervals with structure ((a,b), label) by start (a)
    auto left_comparator = [](const pair <interval_t,unordered_set<T> >& a, const pair <interval_t,unordered_set<T> >& b){
        return a.first.first < b.first.first;
    };

    // How to sort intervals with structure (a,b) by end (b)
    auto right_comparator = [](const interval_t& a, const interval_t& b){
        return a.second < b.second;
    };

    sort(labeled_intervals.begin(), labeled_intervals.end(), left_comparator);

    // Only need to store/compare the interval (and not the data/label) for this DS
    set<interval_t, decltype(right_comparator)> active_intervals;

    for (const auto& [interval,value] : labeled_intervals){
        if (interval.first > interval.second){
            throw runtime_error("ERROR: interval start is greater than interval stop: " + to_string(interval.first) + ',' + to_string(interval.second));
        }

        // Add the interval nodes to the graph, mapped by their interval (pair)
        auto [iter,success] = nodes.try_emplace(interval,value);
        if (not success){
            iter->second.update_values(value);
        }

        vector<interval_t> to_be_removed;

        // Iterate active intervals, which are maintained sorted by interval stop
        for (const auto& other_interval: active_intervals){
            // Flag any expired intervals for deletion
            if (other_interval.second < interval.first){
                to_be_removed.emplace_back(other_interval);
            }
            // Add edges from all active intervals ot the current interval
            else{
                nodes.at(other_interval).neighbors.emplace(interval);
                nodes.at(interval).neighbors.emplace(other_interval);
            }
        }

        // Remove the intervals that have been passed already in the sweep
        for (const auto& item: to_be_removed){
            active_intervals.erase(item);
        }

        active_intervals.emplace(interval);
    }
}


template <class T> IntervalNode<T> IntervalGraph<T>::get_node(const interval_t& interval){
    return nodes.at(interval);
}


template <class T> void IntervalGraph<T>::for_value_in_interval(const interval_t& interval, const function<void(const T& value)>& f) const{
    for (const auto& value: nodes.at(interval).values){
        f(value);
    }
}


template <class T> void IntervalGraph<T>::get_connected_components(vector <unordered_set<interval_t> >& components) {
    unordered_set<interval_t> unvisited;
    components.clear();

    for (const auto& item: nodes){
        unvisited.emplace(item.first);
    }

    while (not unvisited.empty()){
        queue<interval_t> q;

        // Pop arbitrary first node from unvisited
        auto iter = unvisited.begin();
        q.emplace(*iter);
        unvisited.erase(*iter);

        // Initialize new connected component
        components.emplace_back();

        while (not q.empty()){
            auto n = q.front();
            q.pop();

            components.back().emplace(n);

            for (const auto& n_other: nodes.at(n).neighbors){
                if (unvisited.find(n_other) != unvisited.end()){
                    q.emplace(n_other);
                    unvisited.erase(n_other);
                }
            }
        }
    }
}







}

