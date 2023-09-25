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


void run_command(string& command, string& result){
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
}


system_clock::time_point get_current_time() {
    auto now = std::chrono::system_clock::now();
    return now;
}


}
