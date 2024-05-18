#pragma once
#include "Region.hpp"
#include <unordered_set>
#include <set>
#include <utility>
#include <chrono>
#include <string>
#include <array>
#include <limits>

using std::chrono::system_clock;
using std::numeric_limits;
using std::unordered_set;
using std::set;
using std::string;
using std::array;
using std::pair;

#include <filesystem>
using std::filesystem::path;


namespace sv_merge{

using coord_t = pair<int32_t,int32_t>;

using named_coord_t = pair <string, coord_t>;

using interval_t = pair<int32_t,int32_t>;

using labeled_interval_t = pair <interval_t, unordered_set<string> >;

static const interval_t max_placeholder = {numeric_limits<int64_t>::max(),numeric_limits<int64_t>::min()};

inline std::string& ltrim(std::string& s, const char* t = " \t\n\r\f\v");

inline std::string& rtrim(std::string& s, const char* t = " \t\n\r\f\v");

inline std::string& trim(std::string& s, const char* t = " \t\n\r\f\v");

void write_graphaligner_time_log(path output_dir, string vcf_name_prefix, string time_csv, bool success);

void run_command(string& command, string& result, bool trim_result=true);

void run_command(string& command, path output_path);

void run_command(string& command, bool redirect_stderr=true);

uint64_t get_peak_memory_usage();

/**
 * Run a command with a timeout
 * @param command
 * @param redirect_stderr whether or not to redirect the output of the command to the stderr of this executable
 * @param timeout duration in seconds to allow the command to run
 * @return a boolean indicating true if the command ran to completion, false if it was timed out. An exception is thrown
 * if the command itself failed, so there will be no return value in that case.
 */
bool run_command(string command, bool redirect_stderr, float timeout);

system_clock::time_point get_current_time();

bool files_equal(path p1, path p2);

bool equal_ignore_case(const string& str1, const string& str2);

void lowercase_string(string& str);

bool less_than(int32_t a, int32_t b, int32_t c, int32_t d, bool or_equal);

bool less_than(int32_t a, int32_t b, int32_t c, bool or_equal);

bool point_is_contained(int32_t p, const coord_t& i, bool or_equal);

bool point_is_contained(int32_t p, const Region& r, bool or_equal);


}
