#pragma once

#include "pair_hash.hpp"
#include "misc.hpp"

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

/**
 * Simple node object which represents the connections in an IntervalGraph.
 * @tparam T - whatever datatype the intervals are labeled with, usually 'string'
 */
template<class T> class IntervalNode{
public:
    unordered_set<T> values;
    unordered_set<interval_t> neighbors;

    IntervalNode()=default;
    explicit IntervalNode<T>(const unordered_set<T>& value);
    void update_values(const unordered_set<T>& value);
};


/**
 * Add labels to an interval (update its contents)
 * @tparam T - whatever datatype the intervals are labeled with, usually 'string'
 * @param value - set of labels that pertain to a given interval
 */
template <class T> void IntervalNode<T>::update_values(const unordered_set<T>& value){
    for(const auto& v: value){
        values.emplace(v);
    }
}


template <class T> IntervalNode<T>::IntervalNode(const unordered_set<T>& value):
    values(value),
    neighbors()
{}


/**
 * A conventional interval graph implemented with a simple sweep algorithm. Does not support dynamically updating the
 * graph, only iterates intervals once in the constructor.
 * @tparam T
 */
template <class T> class IntervalGraph {
    unordered_map <interval_t, IntervalNode<T> > nodes;

public:
    /// Constructor (builds the graph with sweep algorithm)
    explicit IntervalGraph<T>(vector <pair <interval_t,unordered_set<T> > >& intervals);

    /// Computing components
    void get_connected_components(vector <unordered_set<interval_t> >& components) const;
    void for_each_connected_component_interval(const function<void(interval_t& interval, unordered_set<T>& values)>& f) const;

    /// Iterating and accessing
    void for_value_in_interval(const interval_t& interval, const function<void(const T& value)>& f) const;
    IntervalNode<T> get_node(const interval_t& interval);
};


/**
 * Constructor for IntervalGraph, takes a fixed length vector of labeled intervals of the structure ( [a,b], label_set )
 * where label_set can be any size set which contains labels pertaining to an interval.
 * @tparam T - whatever datatype the intervals are labeled with, usually 'string'
 * @param labeled_intervals
 */
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
            // Add edges from all active intervals to the current interval
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


/**
 * Interval graph accessor/iterator function. Iterate over all the values in a set of values for a given node in the interval graph
 * @tparam T - whatever datatype the intervals are labeled with, usually 'string'
 * @param interval - interval_t
 * @param f - lambda function that acts on the values
 */
template <class T> void IntervalGraph<T>::for_value_in_interval(const interval_t& interval, const function<void(const T& value)>& f) const{
    for (const auto& value: nodes.at(interval).values){
        f(value);
    }
}


/**
 * For each connected component in the interval graph, build a set containing its intervals, and append to the
 * components object which is passed by reference. Uses BFS to find components.
 * @tparam T - whatever datatype the intervals are labeled with, usually 'string'
 * @param components - the object to be filled with results
 */
template <class T> void IntervalGraph<T>::get_connected_components(vector <unordered_set<interval_t> >& components) const{
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


/**
 * An iterator which acts on a simplified representation of a connected component of intervals, in which the upper
 * and lower bounds of the intervals participating in the component are returned as an interval_t and all their
 * values are combined into one set
 * @tparam T - whatever datatype the intervals are labeled with, usually 'string'
 * @param f - to be filled in with a corresponding lambda function that operates on the component intervals
 */
template <class T> void IntervalGraph<T>::for_each_connected_component_interval(const function<void(interval_t& interval, unordered_set<T>& values)> &f) const{
    unordered_set<interval_t> unvisited;
    vector <pair<interval_t, unordered_set<T> > > component_intervals;

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
        component_intervals.push_back({max_placeholder, {}});

        while (not q.empty()){
            auto n = q.front();
            q.pop();

            // Update the component interval by overwriting its bounds, if this interval extends outside them
            if (n.first < component_intervals.back().first.first){
                component_intervals.back().first.first = n.first;
            }
            if (n.second > component_intervals.back().first.second){
                component_intervals.back().first.second = n.second;
            }

            // Also store the labels of the component
            for (const auto& v: nodes.at(n).values) {
                component_intervals.back().second.emplace(v);
            }

            for (const auto& n_other: nodes.at(n).neighbors){
                if (unvisited.find(n_other) != unvisited.end()){
                    q.emplace(n_other);
                    unvisited.erase(n_other);
                }
            }
        }

        f(component_intervals.back().first, component_intervals.back().second);
    }

}


}
