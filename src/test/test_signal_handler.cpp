#include <iostream>
using std::cerr;

#ifdef __linux__

#include <stdexcept>
#include <csignal>
#include <atomic>

using std::runtime_error;

using std::signal;
using std::atomic;
using std::cerr;
using std::cout;

extern "C" void signal_handler(int signal) {
    // UNSAFE USAGE OF STD NAMESPACE!!! See https://man7.org/linux/man-pages/man7/signal-safety.7.html
    switch(signal)
    {
        case SIGTERM:
            cerr << "signal SIGTERM received" << '\n';
            break;
        case SIGSEGV:
            cerr << "signal SIGSEGV received" << '\n';
            break;
        case SIGINT:
            cerr << "signal SIGINT received" << '\n';
            break;
        case SIGILL:
            cerr << "signal SIGILL received" << '\n';
            break;
        case SIGABRT:
            cerr << "signal SIGABRT received" << '\n';
            break;
        case SIGFPE:
            cerr << "signal SIGFPE received" << '\n';
            break;
        default:
            break;
    }

    _exit(signal);
}


int main(){
    signal(SIGTERM, signal_handler);
    signal(SIGSEGV, signal_handler);
    signal(SIGINT, signal_handler);
    signal(SIGILL, signal_handler);
    signal(SIGABRT, signal_handler);
    signal(SIGFPE, signal_handler);

    sleep(3);

    throw runtime_error("TEST ERROR");

    return 0;
}

#else

int main(){
    cerr << "This exe is only supported on Linux" << '\n';
    return 1;
}

#endif
