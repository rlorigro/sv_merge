
#include <iostream>
#include <stdexcept>

using std::runtime_error;
using std::cerr;

#include "misc.hpp"

using namespace sv_merge;


int main() {

    // This test can't work on APPLE because timeout doesn't exist
#ifndef __APPLE__
    string command = "sleep 2 && echo done";

    bool success;

    success = run_command(command, true, 1);
    if (not success) {
        cerr << "timed out" << '\n';
    }

    success = run_command(command, true, 3);
    if (not success) {
        cerr << "timed out" << '\n';
    }

    command = "asdjhkajhs";

    try {
        success = run_command(command, true, 3);
    }
    catch (std::exception& e){
        cerr << "operation failed successfully!"  <<'\n';
    }

    if (not success){
        cerr << "timed out" << '\n';
    }

#endif
    return 0;
}
