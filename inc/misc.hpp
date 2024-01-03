#pragma once

#include "Filesystem.hpp"
#include <unordered_set>
#include <set>
#include <utility>
#include <chrono>
#include <string>
#include <array>
#include <limits>

using std::chrono::system_clock;
using ghc::filesystem::path;
using std::numeric_limits;
using std::unordered_set;
using std::set;
using std::string;
using std::array;
using std::pair;


namespace sv_merge{

using coord_t = pair<int32_t,int32_t>;

using named_coord_t = pair <string, coord_t>;

using interval_t = pair<int32_t,int32_t>;

using labeled_interval_t = pair <interval_t, unordered_set<string> >;

static const interval_t max_placeholder = {numeric_limits<int64_t>::max(),numeric_limits<int64_t>::min()};

inline std::string& ltrim(std::string& s, const char* t = " \t\n\r\f\v");

inline std::string& rtrim(std::string& s, const char* t = " \t\n\r\f\v");

inline std::string& trim(std::string& s, const char* t = " \t\n\r\f\v");

void run_command(string& command, string& result, bool trim_result=true);

void run_command(string& command, path output_path);

void run_command(string& command, bool redirect_stderr=true);

system_clock::time_point get_current_time();

bool files_equal(path p1, path p2);

}
