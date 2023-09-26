#pragma once

#include <chrono>
#include <string>
#include <array>

using std::chrono::system_clock;
using std::string;
using std::array;

namespace sv_merge{

inline std::string& ltrim(std::string& s, const char* t = " \t\n\r\f\v");

inline std::string& rtrim(std::string& s, const char* t = " \t\n\r\f\v");

inline std::string& trim(std::string& s, const char* t = " \t\n\r\f\v");

void run_command(string& command, string& result, bool trim_result=true);

system_clock::time_point get_current_time();

}
