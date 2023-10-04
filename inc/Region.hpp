#pragma once

#include <functional>
#include <cstdlib>
#include <string>
#include <array>

using std::function;
using std::string;
using std::array;

#include "Filesystem.hpp"

using ghc::filesystem::path;

namespace sv_merge {

class Region {
public:
    string name;
    int64_t start{};
    int64_t stop{};

    Region(string &region_string);
    Region(string& name, int64_t start, int64_t stop);
    Region(const string& name, int64_t start, int64_t stop);
    Region()=default;
};

}
