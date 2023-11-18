#pragma once

#include <functional>
#include <cstdlib>
#include <string>
#include <array>

using std::function;
using std::string;
using std::array;

#include "Filesystem.hpp"
#include "pair_hash.hpp"

using ghc::filesystem::path;

namespace sv_merge {

class Region {
public:
    string name;
    int64_t start{};
    int64_t stop{};

    string to_string(char sep=':') const;

    explicit Region(string &region_string);
    Region(string& name, int64_t start, int64_t stop);
    Region(const string& name, int64_t start, int64_t stop);
    Region()=default;
    bool operator==(const Region& other) const;
};

}


// Simple recursive hash for Region
template <> struct std::hash<sv_merge::Region> {
    size_t operator()(const sv_merge::Region& r) const {
        size_t h = std::hash<string>()(r.name);
        hash_combine(h, r.start);
        hash_combine(h, r.stop);
        return h;
    }
};
