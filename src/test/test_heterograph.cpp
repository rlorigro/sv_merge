#include <iostream>
#include <stdexcept>

using std::runtime_error;
using std::exception;
using std::cerr;

#include "bdsg/include/bdsg/internal/hash_map.hpp"
#include "HeteroGraph.hpp"

using sv_merge::HeteroGraph;
using sv_merge::HeteroNode;


int main(){
    HeteroGraph<HeteroNode> g;

    string name_a = "a";
    auto node_a = g.add_node(name_a, '*');
    auto node_a2 = g.get_node(name_a);

    cerr << name_a << ' ' <<  node_a.name << ' ' <<  node_a2.name << '\n';

    cerr << '\n' << "TESTING: add node b" << '\n';

    string name_b = "b";
    auto node_b = g.add_node(name_b, '*');
    auto node_b2 = g.get_node(name_b);

    cerr << name_b << ' ' <<  node_b.name << ' ' <<  node_b2.name << '\n';

    if (not (name_b == node_b.name and name_b == node_b2.name and name_b == "b")){
        throw runtime_error("FAIL: b node does not have consistent 'b' name");
    }

    cerr << '\n' << "TESTING: add node a again" << '\n';

    try {
        node_a = g.add_node(name_a, '*');
    }
    catch (exception& e){
        auto expected_error = "ERROR: cannot add node with existing name: a";
        string error = e.what();

        if (error.find(expected_error) != std::string::npos){
            cerr << "PASS: duplicate node error" << '\n';
        }
        else{
            cerr << e.what() << '\n';
            throw runtime_error("FAIL: unexpected error");
        }
    }

    cerr << '\n' << "TESTING: value of node a" << '\n';

    node_a2 = g.get_node(name_a);

    cerr << name_a << ' ' <<  node_a.name << ' ' <<  node_a2.name << '\n';

    if (not (name_a == node_a.name and name_a == node_a2.name and name_a == "a")){
        throw runtime_error("FAIL: a node does not have consistent 'a' name");
    }

    cerr << '\n' << "TESTING: add edge" << '\n';

    g.add_edge(name_a, name_b, 1);

    unordered_set<pair<string,string> > expected_results = {
            {"a","b"},
            {"b","a"}
    };

    g.for_each_edge([&](const string& a, const string& b, float w){
        cerr << a << ',' << b << ':' << w << '\n';
        pair<string,string> e = {a,b};
        auto result = expected_results.find(e);

        if (result == expected_results.end()){
            throw runtime_error("FAIL: unexpected edge found in graph: " + a + "," + b);
        }

        expected_results.erase(result);
    });

    if (not expected_results.empty()){
        throw runtime_error("FAIL: not all expected edges found in graph");
    }

    cerr << '\n' << "TESTING: remove edge" << '\n';

    g.remove_edge(name_a, name_b);

    int counter = 0;
    g.for_each_edge([&](const string& a, const string& b, float w){
        cerr << a << ',' << b << ':' << w << '\n';
        counter++;
    });

    if (counter != 0){
        throw runtime_error("FAIL: not all edges removed from graph");
    }

    cerr << '\n' << "TESTING: add edge again" << '\n';

    g.add_edge(name_a, name_b, 2);

    expected_results = {
            {"a","b"},
            {"b","a"}
    };

    g.for_each_edge([&](const string& a, const string& b, float w){
        cerr << a << ',' << b << ':' << w << '\n';
        pair<string,string> e = {a,b};
        auto result = expected_results.find(e);

        if (result == expected_results.end()){
            throw runtime_error("FAIL: unexpected edge found in graph: " + a + "," + b);
        }

        expected_results.erase(result);
    });

    if (not expected_results.empty()){
        throw runtime_error("FAIL: not all expected edges found in graph");
    }

    cerr << "PASS" << '\n';

    cerr << '\n' << "TESTING: BFS" << '\n';

    g.add_node("c", '*');
    g.add_node("d", '*');
    g.add_node("e", '*');

    g.add_edge("a", "a",3);
    g.add_edge("a", "c",4);
    g.add_edge("c", "d",5);
    g.add_edge("d", "e",6);
    g.add_edge("e", "b",7);

    g.for_node_in_bfs("a", [&](const HeteroNode& node, int64_t id){
        cerr << node.name << '\n';
    });

    cerr << '\n' << "TESTING: type specific neighbor iteration" << '\n';

    g = {};

    g.add_node("a", 'A');
    g.add_node("b", 'B');
    g.add_node("c", 'C');
    g.add_edge("a", "b", 0);
    g.add_edge("c", "b", 0);
    g.add_edge("c", "a", 0);

    cerr << "all 'a' neighbors" << '\n';
    g.for_each_neighbor("a", [&](const HeteroNode& n, int64_t id){
        cerr << n.name << ' ' << n.type << '\n';
    });
    cerr << '\n';

    cerr << "all 'a' neighbors of type B" << '\n';
    g.for_each_neighbor_of_type("a", 'B', [&](const HeteroNode& n, int64_t id){
        cerr << n.name << ' ' << n.type << '\n';
        if (n.name != "b" and n.type != 'B'){
            throw runtime_error("FAIL: unexpected neighbor");
        }
    });
    cerr << '\n';

    cerr << "all 'b' neighbors of type A" << '\n';
    g.for_each_neighbor_of_type("b", 'A', [&](const HeteroNode& n, int64_t id){
        cerr << n.name << ' ' << n.type << '\n';
        if (n.name != "a" and n.type != 'A'){
            throw runtime_error("FAIL: unexpected neighbor");
        }
    });
    cerr << '\n';


    cerr << "PASS" <<'\n';
}
