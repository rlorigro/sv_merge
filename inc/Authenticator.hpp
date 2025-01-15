#pragma once

#include <filesystem>
#include <functional>
#include <stdexcept>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <chrono>
#include <string>
#include <mutex>

using std::filesystem::path;
using std::runtime_error;
using std::chrono::system_clock;
using std::chrono::seconds;
using std::chrono::hours;
using std::to_string;
using std::string;
using std::cerr;
using std::cout;
using std::mutex;
using std::function;

#include "misc.hpp"

using sv_merge::run_command;
using sv_merge::get_current_time;


namespace sv_merge{


class GoogleAuthenticator{
    mutex m;
    string token = "NULL";

    // Initialize with guaranteed expired timepoint
    system_clock::time_point expiration = get_current_time() - hours(128);

public:
    /// Methods
    void update();
    void try_with_authentication(int64_t n_retries, const function<void()>& f);
};


}
