#include <iostream>
#include <stdexcept>

using std::runtime_error;
using std::exception;
using std::cerr;

#include "bdsg/include/bdsg/internal/hash_map.hpp"
#include "TransitiveMap.hpp"

using sv_merge::TransMap;
using sv_merge::HeteroNode;



int main(){
    TransMap transmap;

    transmap.add_read("read_01");
    transmap.add_sample("HG001");

    transmap.add_read("read_02");
    transmap.add_read("read_03");
    transmap.add_sample("HG002");

    transmap.add_path("a");
    transmap.add_path("b");
    transmap.add_path("c");

    transmap.add_edge("read_01", "HG001");
    transmap.add_edge("read_02", "HG002");
    transmap.add_edge("read_03", "HG002");

    transmap.add_edge("read_01", "a", 2);
    transmap.add_edge("read_01", "b", 1);
    transmap.add_edge("read_03", "b", 1);

    transmap.add_edge("read_02", "b", 1);
    transmap.add_edge("read_02", "c", 2);

    unordered_set<string> visited;

    cerr << '\n' << "TESTING: sample --> read" << '\n';
    transmap.for_each_read_of_sample("HG001", [&](const string& r, int64_t id){
        cerr << r << '\n';
        if (r != "read_01"){
            throw runtime_error("FAIL: unexpected neighbor");
        }

        // Verify no duplicates using visited set
        if (visited.find(r) != visited.end()){
            throw runtime_error("FAIL: duplicate neighbor");
        }
        visited.insert(r);
    });
    cerr << '\n';

    visited.clear();
    transmap.for_each_read_of_sample("HG002", [&](const string& r, int64_t id){
        cerr << r << '\n';
        if (r != "read_02" and r != "read_03"){
            throw runtime_error("FAIL: unexpected neighbor");
        }

        // Verify no duplicates using visited set
        if (visited.find(r) != visited.end()){
            throw runtime_error("FAIL: duplicate neighbor");
        }
        visited.insert(r);
    });
    cerr << '\n';

    cerr << '\n' << "TESTING: path --> sample" << '\n';

    visited.clear();
    transmap.for_each_sample_of_path("a", [&](const string& s, int64_t id){
        cerr << "a" << ' ' << s << '\n';
        if (s != "HG001"){
            throw runtime_error("FAIL: unexpected neighbor");
        }

        // Verify no duplicates using visited set
        if (visited.find(s) != visited.end()){
            throw runtime_error("FAIL: duplicate neighbor");
        }
        visited.insert(s);
    });
    cerr << '\n';

    visited.clear();
    transmap.for_each_sample_of_path("b", [&](const string& s, int64_t id){
        cerr << "b" << ' ' << s << '\n';
        if (not (s == "HG001" or s == "HG002")){
            throw runtime_error("FAIL: unexpected neighbor");
        }

        // Verify no duplicates using visited set
        if (visited.find(s) != visited.end()){
            throw runtime_error("FAIL: duplicate neighbor");
        }
        visited.insert(s);
    });
    cerr << '\n';

    visited.clear();
    transmap.for_each_sample_of_path("c", [&](const string& s, int64_t id){
        cerr << "c" << ' ' << s << '\n';
        if (s != "HG002"){
            throw runtime_error("FAIL: unexpected neighbor");
        }

        // Verify no duplicates using visited set
        if (visited.find(s) != visited.end()){
            throw runtime_error("FAIL: duplicate neighbor");
        }
        visited.insert(s);
    });
    cerr << '\n';

    cerr << '\n' << "TESTING: sample --> path" << '\n';

    visited.clear();
    transmap.for_each_path_of_sample("HG001", [&](const string& p, int64_t id){
        cerr << p << '\n';
        if (not (p == "a" or p == "b")){
            throw runtime_error("FAIL: unexpected neighbor");
        }

        // Verify no duplicates using visited set
        if (visited.find(p) != visited.end()){
            throw runtime_error("FAIL: duplicate neighbor");
        }
        visited.insert(p);

    });
    cerr << '\n';

    visited.clear();
    transmap.for_each_path_of_sample("HG002", [&](const string& p, int64_t id){
        cerr << p << '\n';
        if (not (p == "c" or p == "b")){
            throw runtime_error("FAIL: unexpected neighbor");
        }
        // Verify no duplicates using visited set
        if (visited.find(p) != visited.end()){
            throw runtime_error("FAIL: duplicate neighbor");
        }
        visited.insert(p);
    });
    cerr << '\n';

    cerr << '\n' << "TESTING: global iteration" << '\n';

    visited.clear();
    transmap.for_each_path([&](const string& p, int64_t id){
        cerr << p << '\n';
        if (not (p == "a" or p == "b" or p == "c")){
            throw runtime_error("FAIL: unexpected path");
        }

        // Verify no duplicates using visited set
        if (visited.find(p) != visited.end()){
            throw runtime_error("FAIL: duplicate neighbor");
        }
        visited.insert(p);
    });
    cerr << '\n';

    visited.clear();
    transmap.for_each_read([&](const string& r, int64_t id){
        cerr << r << '\n';
        if (not (r == "read_01" or r == "read_02" or r == "read_03")){
            throw runtime_error("FAIL: unexpected read");
        }

        // Verify no duplicates using visited set
        if (visited.find(r) != visited.end()){
            throw runtime_error("FAIL: duplicate neighbor");
        }
        visited.insert(r);
    });
    cerr << '\n';

    visited.clear();
    transmap.for_each_sample([&](const string& s, int64_t id){
        cerr << s << '\n';
        if (not (s == "HG001" or s == "HG002")){
            throw runtime_error("FAIL: unexpected sample");
        }

        // Verify no duplicates using visited set
        if (visited.find(s) != visited.end()){
            throw runtime_error("FAIL: duplicate neighbor");
        }
        visited.insert(s);
    });
    cerr << '\n';

    cerr << "TESTING try_get_ methods" << '\n';

    cerr << "id1" << '\n';
    string name1 = "read_01";
    auto [success1,id1] = transmap.try_get_id(name1);
    cerr << success1 << ',' <<  id1 << '\n';
    if (not (success1 == 1 and id1 == 4)){
        throw runtime_error("FAIL: incorrect ID found in graph");
    }

    cerr << "id2" << '\n';
    string name2 = "read_02";
    auto [success2,id2] = transmap.try_get_id(name2);
    cerr << success2 << ',' <<  id2 << '\n';
    if (not (success2 == 1 and id2 == 6)){
        throw runtime_error("FAIL: incorrect ID found in graph");
    }

    cerr << "e12" << '\n';
    auto [success_e12, w12] = transmap.try_get_edge_weight(id1, id2);
    if (not (success_e12 == 0 and w12 == 0)){
        throw runtime_error("FAIL: incorrect edge found in graph");
    }

    transmap.add_edge(id1, id2, 0.666);

    cerr << "e12" << '\n';
    std::tie(success_e12, w12) = transmap.try_get_edge_weight(id1, id2);
    if (not (success_e12 == 1 and int(w12*1000) == 666)){
        throw runtime_error("FAIL: incorrect ID found in graph");
    }

    cerr << success_e12 << ',' << w12 << '\n';

    cerr << "TESTING: exhaustive iteration" << '\n';
    transmap.for_each_sample([&](const string& sample_name, int64_t sample_id){
        cerr << "sample: " << sample_name << '\n';
        transmap.for_each_read_of_sample(sample_id, [&](int64_t read_id){
            cerr << "read: " << transmap.get_node(read_id).name << '\n';
            transmap.for_each_path_of_read(read_id, [&](int64_t path_id){
                cerr << "path: " << transmap.get_node(path_id).name << '\n';
                auto [success, weight] = transmap.try_get_edge_weight(read_id, path_id);

                if (not success){
                    throw runtime_error("ERROR: edge weight not found for read-path edge: " + to_string(read_id) + " " + to_string(path_id));
                }

                cerr << sample_name << ',' << transmap.get_node(read_id).name << ',' << transmap.get_node(path_id).name << ',' << weight << '\n';
            });
        });
    });

    cerr << "PASS" << '\n';

    // test removal of read nodes and verify that the remaining set of nodes/edges is correct

    auto id = transmap.get_id("read_01");
    transmap.remove_node(id);

    cerr << "TESTING: removal of read_01" << '\n';

    // first iterate nodes and verify that read_01 is not present
    bool read_found = false;
    transmap.for_each_read([&](const string& name, int64_t id){
        if (name == "read_01"){
            read_found = true;
        }

        cerr << name << '\n';
    });

    if (read_found){
        throw runtime_error("FAIL: read_01 was not removed");
    }

    read_found = false;

    // Then iterate all samples and verify that read_01 is not connected to any of them
    transmap.for_each_sample([&](const string& sample_name, int64_t sample_id){
        transmap.for_each_read_of_sample(sample_id, [&](int64_t read_id){
            if (transmap.get_node(read_id).name == "read_01"){
                read_found = true;
            }

            cerr << sample_name << ',' << transmap.get_node(read_id).name << '\n';
        });
    });

    if (read_found){
        throw runtime_error("FAIL: read_01 was not removed");
    }

}
