#include <iostream>
#include <stdexcept>

using std::runtime_error;
using std::exception;
using std::cerr;

#include "bdsg/include/bdsg/internal/hash_map.hpp"
#include "TransitiveMap.hpp"
#include "optimize/path_optimizer.hpp"

using sv_merge::TransMap;
using sv_merge::HeteroNode;
using sv_merge::optimize_d;


void test_optimization(){
    TransMap transmap;

    transmap.add_sample("HG001");
    transmap.add_sample("HG002");

    transmap.add_read("read_01");
    transmap.add_read("read_02");
    transmap.add_read("read_03");
    transmap.add_read("read_04");
    transmap.add_read("read_05");
    transmap.add_read("read_06");

    transmap.add_path("a");
    transmap.add_path("b");
    transmap.add_path("c");

    transmap.add_edge("read_01", "HG001");
    transmap.add_edge("read_03", "HG001");
    transmap.add_edge("read_05", "HG001");

    transmap.add_edge("read_02", "HG002");
    transmap.add_edge("read_04", "HG002");
    transmap.add_edge("read_06", "HG002");

    // HG001 (odd)
    transmap.add_edge("read_01", "a", 1);
    transmap.add_edge("read_01", "b", 2);
    transmap.add_edge("read_01", "c", 3);

    transmap.add_edge("read_03", "a", 1);
    transmap.add_edge("read_03", "b", 2);
    transmap.add_edge("read_03", "c", 3);

    transmap.add_edge("read_05", "a", 1);
    transmap.add_edge("read_05", "b", 2);
    transmap.add_edge("read_05", "c", 3);

    // HG002 (even)
    transmap.add_edge("read_02", "a", 3);
    transmap.add_edge("read_02", "b", 2);
    transmap.add_edge("read_02", "c", 1);

    transmap.add_edge("read_04", "a", 2);
    transmap.add_edge("read_04", "b", 1);
    transmap.add_edge("read_04", "c", 3);

    transmap.add_edge("read_06", "a", 1);
    transmap.add_edge("read_06", "b", 2);
    transmap.add_edge("read_06", "c", 3);

    optimize_d(transmap);
    optimize_d_plus_n(transmap, 1, 10);
    optimize_d_plus_n(transmap, 1, 0);
    optimize_d_plus_n(transmap, 10, 1);
}


int main(){
    test_optimization();
}
