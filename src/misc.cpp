#include "misc.hpp"

#include <fstream>
#include <string>
#include <iostream>
#include <cstdio>
#include <array>
#include <stdexcept>
#include <random>
#include <vector>

using std::random_device;
using std::uniform_int_distribution;
using std::mt19937;

using std::runtime_error;
using std::chrono::seconds;
using std::chrono::hours;
using std::to_string;
using std::ofstream;
using std::ifstream;
using std::string;
using std::vector;
using std::cerr;
using std::cout;


using std::runtime_error;


namespace sv_merge{


/**
 * Append a log file and write the header if it hasn't been written yet
 * @param output_dir
 * @param name
 * @param time_csv result of calling Timer::to_csv() immediately after task exits
 * @param success whether or not the task timed out
 */
void write_time_log(const path& output_dir, const string& name, const string& time_csv, bool success, const string& notes){
    // Begin the logging process
    path log_path = output_dir / "log.csv";

    // Check if the log file needs to have a header written to it or not
    bool exists = std::filesystem::exists(log_path);

    ofstream file(log_path, std::ios_base::app);

    // Write the header
    if (not exists){
        file << "name,h,m,s,ms,success,notes" << '\n';
    }
    // Write the results for this region/tool
    file << name << ',' << time_csv << ',' << success << ',' << notes << '\n';
}


/**
 * Append a log file and write the header if it hasn't been written yet
 * @param output_dir
 * @param name
 * @param timer timer which has been running for the period of the task to be logged (will be polled immediately)
 * @param success whether or not the task timed out
 */
void write_time_log(const path& output_dir, const string& name, const Timer& timer, bool success, const string& notes){
    string time_csv = timer.to_csv();

    write_time_log(output_dir, name, time_csv, success, notes);
}


/**
 * Append a log file and write the header if it hasn't been written yet
 * @param output_dir
 * @param name
 * @param duration time of task in milliseconds
 * @param success whether or not the task timed out
 */
void write_time_log(const path& output_dir, const string& name, const milliseconds& duration, bool success, const string& notes){
    string time_csv;
    duration_to_csv(duration, time_csv);

    write_time_log(output_dir, name, time_csv, success, notes);
}


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
    uint64_t peak_memory_usage = 0ULL;

    ifstream procStats("/proc/self/status");
    if (procStats) {
        string line;
        while (std::getline(procStats, line)) {
            if (string::npos == line.find("VmHWM")) {
                continue;
            }
            size_t pos = line.find(':');
            while (pos < line.size() && !isdigit(line[pos])) {
                pos++;
            }
            char* end;
            peak_memory_usage = std::strtoull(line.c_str() + pos, &end, 10);
            break;
        }
    }

    return peak_memory_usage;
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


bool less_than(int32_t a, int32_t b, int32_t c, int32_t d, bool or_equal){
    if (or_equal){
        return (a <= b) and (b <= c) and (c <= d);
    }
    else{
        return (a < b) and (b < c) and (c < d);
    }
}


bool less_than(int32_t a, int32_t b, int32_t c, bool or_equal){
    if (or_equal){
        return (a <= b) and (b <= c);
    }
    else{
        return (a < b) and (b < c);
    }
}


bool point_is_contained(int32_t p, const coord_t& i, bool or_equal){
    return less_than(i.first, p, i.second, or_equal) or less_than(i.second, p, i.first, or_equal);
}


bool point_is_contained(int32_t p, const Region& r, bool or_equal){
    return less_than(r.start, p, r.stop, or_equal) or less_than(r.stop, p, r.start, or_equal);
}


// Taken from:
// https://stackoverflow.com/a/58467162
//
// WARNING: This function is not thread-safe!
//
//       ^
//      / \
//     / ! \
//    /_____\
//
//
string get_uuid() {
    static random_device dev;
    static mt19937 rng(dev());

    uniform_int_distribution<int> dist(0, 15);

    const char *v = "0123456789abcdef";
    const bool dash[] = {0,0,0,0,1,0,1,0,1,0,1,0,0,0,0,0};

    string res;
    res.reserve(16);
    for (int i = 0; i < 16; i++) {
        if (dash[i]) res += "-";
        res += v[dist(rng)];
        res += v[dist(rng)];
    }
    return res;
}


void for_each_row_in_csv(path csv_path, const function<void(const vector<string>& items)>& f){
    if (not (csv_path.extension() == ".csv")){
        throw runtime_error("ERROR: file does not have compatible csv extension: " + csv_path.string());
    }

    ifstream file(csv_path);

    if (not (file.is_open() and file.good())){
        throw runtime_error("ERROR: could not read file: " + csv_path.string());
    }

    char c;
    vector<string> items = {""};

    int64_t n_char_in_line = 0;
    char delimiter = ',';

    while (file.get(c)){
        if (c == delimiter){
            items.emplace_back();
            continue;
        }
        if (c == '\r'){
            throw runtime_error("ERROR: carriage return not supported: " + csv_path.string());
        }

        if (c == '\n'){
            if (n_char_in_line == 0){
                continue;
            }

            f(items);

            items.resize(1);
            items[0].clear();
            n_char_in_line = 0;
            continue;
        }

        items.back() += c;

        n_char_in_line++;
    }

    if (n_char_in_line > 0) {
        f(items);
    }
}


}
