#include <iostream>
#include <stdexcept>

using std::runtime_error;
using std::cerr;


#include "IntervalGraph.hpp"

using sv_merge::IntervalNode;
using sv_merge::IntervalGraph;
using sv_merge::interval_t;


int main(){
    vector <pair <interval_t,unordered_set<string> > > intervals = {
            {{0,10}, {"a"}},
            {{0,10}, {"b"}},
            {{0,10}, {"c"}},
            {{10,20}, {"d"}},
            {{11,20}, {"e"}},
            {{11,0}, {"f"}},
    };

    IntervalGraph<string> g(intervals);



    return 0;
}
