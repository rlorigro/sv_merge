#include "HalfInterval.hpp"


namespace sv_merge {


HalfInterval::HalfInterval(size_t id, int32_t position, bool is_start) :
        id(id),
        position(position),
        is_start(is_start)
{}


void deoverlap_intervals(const vector<interval_t>& intervals, vector<interval_t>& result_intervals, unordered_map<interval_t, size_t>& result_mapping){
    if (intervals.size() <= 1){
        return;
    }

    vector <HalfInterval> half_intervals;
    vector <int32_t> positions;
    vector <unordered_set<size_t> > ids;

    for (size_t i=0; i<intervals.size(); i++){
        const auto& a = intervals[i];
        half_intervals.emplace_back(i,a.first,true);
        half_intervals.emplace_back(i,a.second,false);
    }

    // How to sort labeled intervals with structure ((a,b), label) by start (a)
    auto left_comparator = [](const HalfInterval& a, const HalfInterval& b){
        if (a.position == b.position){
            return a.is_start > b.is_start;
        }
        else {
            return a.position < b.position;
        }
    };

    sort(half_intervals.begin(), half_intervals.end(), left_comparator);

    const auto& x = half_intervals[0];
    positions.emplace_back(x.position);
    ids.push_back({x.id});

    // With at least 2 intervals, should be guaranteed 4 half intervals in the vector
    for (size_t i=1; i<half_intervals.size(); i++){
        const auto& h = half_intervals[i];

//        cerr << "-- " << h.position << ',' << (h.is_start ? '[' : ')') << ',' << h.id << '\n';

        // Every time an element is added or removed, update the vector of intersected_intervals:
        //  1. Extend the previous interval
        //  2. Remove/add an element from the set (check if the start/stop is not identical to prev)
        if (h.is_start){
            // Only terminate the previous interval if this one does not have the same position
            if (h.position > positions.back()) {
                // Add new breakpoint
                positions.emplace_back(h.position);
                ids.emplace_back(ids.back());
            }
            // ADD item to last set
            ids.back().emplace(h.id);
        }
        else{
            // Only terminate the previous interval if this one does not have the same position
            if (h.position > positions.back()) {
                // Add new breakpoint
                positions.emplace_back(h.position);
                ids.emplace_back(ids.back());
            }
            // REMOVE item from last set
            ids.back().erase(h.id);
        }
    }

    for (size_t i=0; i<ids.size() - 1; i++) {
        if (ids[i].empty()) {
            continue;
        }

        auto a = positions[i];
        auto b = positions[i+1];

        for (auto& id: ids[i]){
            result_mapping.emplace(interval_t{a,b}, id);
        }

        result_intervals.emplace_back(a,b);
    }

}

}