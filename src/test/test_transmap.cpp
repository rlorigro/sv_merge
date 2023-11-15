#include <iostream>
#include <stdexcept>

using std::runtime_error;
using std::exception;
using std::cerr;

#include "pair_hash.hpp"
#include "TransitiveMap.hpp"
#include "optimizer.hpp"

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
    TransMap transmap;

    transmap.add_read("read_01");
    transmap.add_sample("HG001");

    transmap.add_read("read_02");
    transmap.add_sample("HG002");

    transmap.add_path("a");
    transmap.add_path("b");
    transmap.add_path("c");

    transmap.add_edge("read_01", "HG001");
    transmap.add_edge("read_02", "HG002");

    transmap.add_edge("read_01", "a", 2);
    transmap.add_edge("read_01", "b", 1);

    transmap.add_edge("read_02", "b", 1);
    transmap.add_edge("read_02", "c", 2);

    cerr << '\n' << "TESTING: sample --> read" << '\n';
    transmap.for_each_read_of_sample("HG001", [&](const string& r, int64_t id){
        cerr << r << '\n';
        if (r != "read_01"){
            throw runtime_error("FAIL: unexpected neighbor");
        }
    });
    cerr << '\n';

    transmap.for_each_read_of_sample("HG002", [&](const string& r, int64_t id){
        cerr << r << '\n';
        if (r != "read_02"){
            throw runtime_error("FAIL: unexpected neighbor");
        }
    });
    cerr << '\n';

    cerr << '\n' << "TESTING: path --> sample" << '\n';
    transmap.for_each_sample_of_path("a", [&](const string& s, int64_t id){
        cerr << "a" << ' ' << s << '\n';
        if (s != "HG001"){
            throw runtime_error("FAIL: unexpected neighbor");
        }
    });
    cerr << '\n';

    transmap.for_each_sample_of_path("b", [&](const string& s, int64_t id){
        cerr << "b" << ' ' << s << '\n';
        if (not (s == "HG001" or s == "HG002")){
            throw runtime_error("FAIL: unexpected neighbor");
        }
    });
    cerr << '\n';

    transmap.for_each_sample_of_path("c", [&](const string& s, int64_t id){
        cerr << "c" << ' ' << s << '\n';
        if (s != "HG002"){
            throw runtime_error("FAIL: unexpected neighbor");
        }
    });
    cerr << '\n';

    cerr << '\n' << "TESTING: sample --> path" << '\n';
    transmap.for_each_path_of_sample("HG001", [&](const string& p, int64_t id){
        cerr << p << '\n';
        if (not (p == "a" or p == "b")){
            throw runtime_error("FAIL: unexpected neighbor");
        }
    });
    cerr << '\n';

    transmap.for_each_path_of_sample("HG002", [&](const string& p, int64_t id){
        cerr << p << '\n';
        if (not (p == "c" or p == "b")){
            throw runtime_error("FAIL: unexpected neighbor");
        }
    });
    cerr << '\n';

    cerr << '\n' << "TESTING: global iteration" << '\n';
    transmap.for_each_path([&](const string& p, int64_t id){
        cerr << p << '\n';
        if (not (p == "a" or p == "b" or p == "c")){
            throw runtime_error("FAIL: unexpected path");
        }
    });
    cerr << '\n';

    transmap.for_each_read([&](const string& r, int64_t id){
        cerr << r << '\n';
        if (not (r == "read_01" or r == "read_02")){
            throw runtime_error("FAIL: unexpected read");
        }
    });
    cerr << '\n';

    transmap.for_each_sample([&](const string& s, int64_t id){
        cerr << s << '\n';
        if (not (s == "HG001" or s == "HG002")){
            throw runtime_error("FAIL: unexpected sample");
        }
    });
    cerr << '\n';

    cerr << "TESTING try_get_ methods" << '\n';

    cerr << "id1" << '\n';
    string name1 = "read_01";
    auto [success1,id1] = transmap.try_get_id(name1);
    if (not (success1 == 1 and id1 == 3)){
        throw runtime_error("FAIL: non existent ID found in graph");
    }

    cerr << "id2" << '\n';
    string name2 = "read_02";
    auto [success2,id2] = transmap.try_get_id(name2);
    if (not (success2 == 1 and id2 == 5)){
        throw runtime_error("FAIL: non existent ID found in graph");
    }

    cerr << "e12" << '\n';
    auto [success_e12, w12] = transmap.try_get_edge_weight(id1, id2);
    if (not (success_e12 == 0 and w12 == 0)){
        throw runtime_error("FAIL: non existent edge found in graph");
    }

    transmap.add_edge(id1, id2, 0.666);

    cerr << "e12" << '\n';
    std::tie(success_e12, w12) = transmap.try_get_edge_weight(id1, id2);
    if (not (success_e12 == 1 and int(w12*1000) == 666)){
        throw runtime_error("FAIL: non existent ID found in graph");
    }

    cerr << success_e12 << ',' << w12 << '\n';

    test_optimization();

    cerr << "PASS" << '\n';
}
