#pragma once

#include <string>
#include <iostream>
#include <cstdlib>
#include <stdexcept>
#include <chrono>
#include <iomanip>

using std::runtime_error;
using std::chrono::system_clock;
using std::chrono::seconds;
using std::chrono::hours;
using std::to_string;
using std::string;
using std::cerr;
using std::cout;

#include "Filesystem.hpp"
#include "misc.hpp"

using sv_merge::run_command;
using sv_merge::get_current_time;
using ghc::filesystem::path;


namespace sv_merge{


///
/// Consider using this: https://stackoverflow.com/a/57282480
/// To check for expiration
///
class GoogleAuthenticator{
    /// Attributes
    system_clock::time_point expiration = get_current_time() - hours(128);

    // 60min duration by default
    system_clock::duration token_lifetime = seconds(3600);

public:
    /// Methods
    void update();
};

}
