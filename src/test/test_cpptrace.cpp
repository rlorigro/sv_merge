#include <iostream>
#include <stdexcept>

using std::runtime_error;
using std::cerr;

#include "cpptrace/from_current.hpp"

void foo() {
    throw std::runtime_error("foo failed");
}
int main() {
    CPPTRACE_TRY {
        foo();
    } CPPTRACE_CATCH(const std::exception& e) {
        std::cerr<<"Exception: "<<e.what()<<std::endl;
        cpptrace::from_current_exception().print();
    }
}
