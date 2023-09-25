#include "misc.hpp"

#include <string>
#include <iostream>
#include <cstdio>
#include <array>
#include <stdexcept>
#include <iomanip>

using std::runtime_error;
using std::chrono::seconds;
using std::chrono::hours;
using std::to_string;
using std::string;
using std::cerr;
using std::cout;

#include "Filesystem.hpp"

using ghc::filesystem::path;
using std::runtime_error;


namespace hapslap{


void run_command(string& command, string& result, bool trim_result){
    array<char, 128> buffer;
    buffer.fill('\0');

    FILE* pipe = popen(command.c_str(), "r");

    if (!pipe){
        throw runtime_error("Pipe could not be opened: " + command);
    }

    while (fgets(buffer.data(), 128, pipe) != nullptr){
        result += buffer.data();
    }

    auto return_code = pclose(pipe);

    if (return_code != 0){
        throw runtime_error("Command failed: " + command);
    }

    if (trim_result){
        trim(result);
    }
}


system_clock::time_point get_current_time() {
    auto now = std::chrono::system_clock::now();
    return now;
}


/// TRIM FUNCTIONS taken from https://stackoverflow.com/a/25829233
// trim from left
inline std::string& ltrim(std::string& s, const char* t){
    s.erase(0, s.find_first_not_of(t));
    return s;
}

// trim from right
inline std::string& rtrim(std::string& s, const char* t){
    s.erase(s.find_last_not_of(t) + 1);
    return s;
}

// trim from left & right
inline std::string& trim(std::string& s, const char* t){
    return ltrim(rtrim(s, t), t);
}

}
