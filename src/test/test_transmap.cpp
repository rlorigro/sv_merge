#include <iostream>
#include <stdexcept>

using std::runtime_error;
using std::exception;
using std::cerr;

#include "pair_hash.hpp"
#include "TransitiveMap.hpp"

using sv_merge::TransMap;
using sv_merge::HeteroNode;


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

    transmap.add_edge("read_01", "a");
    transmap.add_edge("read_01", "b");

    transmap.add_edge("read_02", "b");
    transmap.add_edge("read_02", "c");

    cerr << '\n' << "TESTING: sample --> read" << '\n';
    transmap.for_each_read_of_sample("HG001", [&](const HeteroNode& r){
        cerr << r.name << '\n';
        if (r.name != "read_01" and r.type != 'R'){
            throw runtime_error("FAIL: unexpected neighbor");
        }
    });
    cerr << '\n';

    transmap.for_each_read_of_sample("HG002", [&](const HeteroNode& r){
        cerr << r.name << '\n';
        if (r.name != "read_02" and r.type != 'R'){
            throw runtime_error("FAIL: unexpected neighbor");
        }
    });
    cerr << '\n';

    cerr << '\n' << "TESTING: sample --> path" << '\n';
    transmap.for_each_sample_of_path("a", [&](const HeteroNode& s){
        cerr << s.name << '\n';
        if (s.name != "read_01" and s.type != 'S'){
            throw runtime_error("FAIL: unexpected neighbor");
        }
    });
    cerr << '\n';

    transmap.for_each_sample_of_path("b", [&](const HeteroNode& s){
        cerr << s.name << '\n';
        if ((s.name != "read_01" or s.name != "read_02") and s.type != 'S'){
            throw runtime_error("FAIL: unexpected neighbor");
        }
    });
    cerr << '\n';

    transmap.for_each_sample_of_path("c", [&](const HeteroNode& s){
        cerr << s.name << '\n';
        if (s.name != "read_02" and s.type != 'S'){
            throw runtime_error("FAIL: unexpected neighbor");
        }
    });
    cerr << '\n';

    cerr << '\n' << "TESTING: path --> sample" << '\n';
    transmap.for_each_path_of_sample("HG001", [&](const HeteroNode& p){
        cerr << p.name << '\n';
        if ((p.name != "a" or p.name != "b") and p.type != 'P'){
            throw runtime_error("FAIL: unexpected neighbor");
        }
    });
    cerr << '\n';

    transmap.for_each_path_of_sample("HG002", [&](const HeteroNode& p){
        cerr << p.name << '\n';
        if ((p.name != "c" or p.name != "b") and p.type != 'P'){
            throw runtime_error("FAIL: unexpected neighbor");
        }
    });
    cerr << '\n';

    cerr << "PASS" << '\n';
}
