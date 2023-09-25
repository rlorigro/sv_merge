#pragma once

#include <chrono>
#include <string>
#include <array>

using std::chrono::system_clock;
using std::string;
using std::array;

namespace hapslap{

void run_command(string& command, string& result);
system_clock::time_point get_current_time();

}
