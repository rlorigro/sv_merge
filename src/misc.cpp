#include "misc.hpp"

#include <fstream>
#include <string>
#include <iostream>
#include <cstdio>
#include <array>
#include <stdexcept>

using std::runtime_error;
using std::chrono::seconds;
using std::chrono::hours;
using std::to_string;
using std::ofstream;
using std::ifstream;
using std::string;
using std::cerr;
using std::cout;


using std::runtime_error;


namespace sv_merge{


void run_command(string& command, string& result, bool trim_result){
    result.clear();

    cerr << "RUNNING: " << command << '\n';

    array<char, 128> buffer;

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


void run_command(string& command, path output_path){
    cerr << "RUNNING: " << command << '\n';
    cerr << "REDIRECTING TO: " << output_path << '\n';

    ofstream file(output_path);

    array<char, 128> buffer;

    FILE* pipe = popen(command.c_str(), "r");

    if (!pipe){
        throw runtime_error("Pipe could not be opened: " + command);
    }

    while (fgets(buffer.data(), 128, pipe) != nullptr){
        file << buffer.data();
    }

    auto return_code = pclose(pipe);

    if (return_code != 0){
        throw runtime_error("Command failed: " + command);
    }
}


void run_command(string& command, bool redirect_stderr){
    if (redirect_stderr){
        command += " 2>&1";
    }

    cerr << "RUNNING: " << command << '\n';

    array<char, 128> buffer;

    FILE* pipe = popen(command.c_str(), "r");

    if (!pipe){
        throw runtime_error("Pipe could not be opened: " + command);
    }

    while (fgets(buffer.data(), 128, pipe) != nullptr){
        cerr << buffer.data();
    }

    auto return_code = pclose(pipe);

    if (return_code != 0){
        throw runtime_error("Command failed: " + command);
    }
}


bool run_command(string command, bool redirect_stderr, float timeout){
    command = "timeout " + to_string(timeout) + ' ' + command;

    cerr << "RUNNING: " << command << '\n';

    array<char, 128> buffer;

    FILE* pipe = popen(command.c_str(), "r");

    if (!pipe){
        throw runtime_error("Pipe could not be opened: " + command);
    }

    while (fgets(buffer.data(), 128, pipe) != nullptr){
        if (redirect_stderr) {
            cerr << buffer.data();
        }
    }

    int result = pclose(pipe);

    int return_code = WEXITSTATUS(result);

    if (return_code == 124){
        return false;
    }
    else if (return_code != 0){
        throw runtime_error("Command failed: " + command);
    }

    return true;
}


uint64_t get_peak_memory_usage() {
    uint64_t peakMemoryUsage = 0ULL;

    ifstream procStats("/proc/self/status");
    if (procStats) {
        string line;
        while (std::getline(procStats, line)) {
            if (string::npos == line.find("VmPeak")) {
                continue;
            }
            size_t pos = line.find(":");
            while (pos < line.size() && !isdigit(line[pos])) {
                pos++;
            }
            char* end;
            peakMemoryUsage = std::strtoull(line.c_str() + pos, &end, 10);
            // Convert from kB to bytes.
            peakMemoryUsage *= 1024;
            break;
        }
    }

    return peakMemoryUsage;
}


// https://stackoverflow.com/questions/6163611/compare-two-files
bool files_equal(path p1, path p2) {
    std::ifstream f1(p1.string(), std::ifstream::binary|std::ifstream::ate);
    std::ifstream f2(p2.string(), std::ifstream::binary|std::ifstream::ate);

    if (f1.fail() || f2.fail()) {
        return false; //file problem
    }

    if (f1.tellg() != f2.tellg()) {
        return false; //size mismatch
    }

    //seek back to beginning and use std::equal to compare contents
    f1.seekg(0, std::ifstream::beg);
    f2.seekg(0, std::ifstream::beg);
    return std::equal(std::istreambuf_iterator<char>(f1.rdbuf()),
                      std::istreambuf_iterator<char>(),
                      std::istreambuf_iterator<char>(f2.rdbuf()));
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


bool equal_ignore_case(const string& str1, const string& str2) {
    const size_t length1 = str1.length();
    const size_t length2 = str2.length();
    if (length1!=length2) return false;
    for (size_t i=0; i<length1; i++) {
        if (tolower(str1.at(i))!=tolower(str2.at(i))) return false;
    }
    return true;
}


void lowercase_string(string& str) {
    const size_t length = str.length();
    for (size_t i=0; i<length; i++) str.at(i)=(char)tolower(str.at(i));
}

}
