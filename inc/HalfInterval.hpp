#pragma once

#include <cstddef>
#include <cstdint>

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>

#include "misc.hpp"

using std::vector;
using std::unordered_map;
using std::unordered_set;
using std::sort;


namespace sv_merge {

class HalfInterval {
public:
    size_t id;
    int32_t position;
    bool is_start;

    HalfInterval(size_t id, int32_t position, bool is_start);
};


void deoverlap_intervals(
        const vector<interval_t>& intervals,
        vector<interval_t>& result_intervals,
        unordered_map<interval_t, size_t>& result_mapping
    );

}
