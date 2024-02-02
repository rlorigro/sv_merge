#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <random>

using std::runtime_error;
using std::ofstream;
using std::cerr;
using std::string;
using std::random_device;
using std::uniform_int_distribution;
using std::mt19937;

#include "Filesystem.hpp"

using ghc::filesystem::path;
using ghc::filesystem::create_directories;

#include <vector>
#include <iostream>
#include <stdexcept>
#include <csignal>
#include <atomic>

using std::vector;
using std::runtime_error;
using std::cerr;
using std::signal;
using std::atomic;
using std::cerr;
using std::cout;


static const char* SHM_WARNING = "WARNING: this program uses shared memory and has terminated unexpectedly, remove the following directory manually to free allocated memory: ";
static const char* SHM_DIRECTORY = "/dev/shm/hapestry";
static const int SHM_WARNING_SIZE = 140;
static const int SHM_DIRECTORY_SIZE = 18;

extern "C" void signal_handler(int signal) {
    write(2,SHM_WARNING, SHM_WARNING_SIZE);
    write(2,SHM_DIRECTORY, SHM_DIRECTORY_SIZE);
    write(2,"\n", 1);
    _exit(signal);
}


// Taken from:
// https://stackoverflow.com/a/58467162
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


int main(){
    signal(SIGTERM, signal_handler);
    signal(SIGSEGV, signal_handler);
    signal(SIGINT, signal_handler);
    signal(SIGILL, signal_handler);
    signal(SIGABRT, signal_handler);
    signal(SIGFPE, signal_handler);

    path temp_dir = SHM_DIRECTORY;
    create_directories(temp_dir);

    path temp_file = temp_dir / (get_uuid() + ".txt");
    ofstream file(temp_file);

    file << "test\n";

    throw runtime_error("TEST");

    return 0;
}
