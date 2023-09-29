#pragma once

#include "pair_hash.hpp"

#include <unordered_map>
#include <unordered_set>
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
using std::set_union;
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
    unordered_set<uint64_t> neighbors;

    IntervalNode<T>()=default;
    explicit IntervalNode<T>(const unordered_set<T>& values);
};


template <class T> IntervalNode<T>::IntervalNode(const unordered_set<T>& v):
    values(),
    neighbors()
{
    for(const auto& item: v){
        values.insert(item);
    }
}


template <class T> class IntervalGraph {
    unordered_map <interval_t, IntervalNode<T> > nodes;

public:
    IntervalGraph<T>(vector <pair <interval_t,unordered_set<T> > >& intervals);
    void get_connected_components(vector <unordered_set<interval_t> >& components);
    void print(ostream& out);
};


template <class T> void IntervalGraph<T>::print(ostream& out) {

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


template <class T> IntervalGraph<T>::IntervalGraph(vector <pair <interval_t,unordered_set<T> > >& intervals) :
        nodes()
{
    auto left_comparator = [](pair <interval_t,unordered_set<T> >& a, pair <interval_t,unordered_set<T> >& b){
        return a.first.first < b.first.first;
    };

    auto right_comparator = [](pair <interval_t,unordered_set<T> >& a, pair <interval_t,unordered_set<T> >& b){
        return a.first.second < b.first.second;
    };

    sort(intervals.begin(), intervals.end(), left_comparator);

    set<interval_t, decltype(right_comparator)> active_intervals;

    for (const auto& [interval,value] : intervals){
        // Add the interval nodes to the graph, mapped by their interval (pair)
        auto result = nodes.try_emplace(interval,IntervalNode(value));

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





}

