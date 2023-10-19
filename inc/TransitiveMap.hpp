#pragma once

#include "HeteroGraph.hpp"

#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <string>
#include <vector>

using std::unordered_map;
using std::unordered_set;
using std::vector;
using std::string;
using std::pair;


namespace sv_merge {


class TransMap {
    /// Attributes
    HeteroGraph<HeteroNode> graph;
    const string sample_node_name;
    const string read_node_name;
    const string path_node_name;

public:
    TransMap();
    void add_sample(const string& name);
    void add_read(const string& name);
    void add_path(const string& name);
    void add_edge(const string& a, const string& b);
    void add_edge(const string& a, const string& b, float weight);

    void get_read_sample(const string& read_name, string& result);
    void for_each_read_of_sample(const string& sample_name, const function<void(const HeteroNode& node)>& f);
    void for_each_sample_of_read(const string& read_name, const function<void(const HeteroNode& node)>& f);
    void for_each_sample_of_path(const string& path_name, const function<void(const HeteroNode& node)>& f);
    void for_each_path_of_sample(const string& sample_name, const function<void(const HeteroNode& node)>& f);
};


TransMap::TransMap():
    sample_node_name("sample_node"),
    read_node_name("read_node"),
    path_node_name("path_node")
{
    auto& sample_node = graph.add_node(sample_node_name);
    sample_node.type = 'S';

    auto& read_node = graph.add_node(read_node_name);
    read_node.type = 'R';

    auto& path_node = graph.add_node(path_node_name);
    path_node.type = 'P';
}


void TransMap::add_read(const string& name){
    if (name == read_node_name or name == sample_node_name or name == path_node_name){
        throw runtime_error("ERROR: cannot add node with preset node name: " + name);
    }
    auto& node = graph.add_node(name);
    node.type = 'R';
}


void TransMap::add_sample(const string& name){
    if (name == read_node_name or name == sample_node_name or name == path_node_name){
        throw runtime_error("ERROR: cannot add node with preset node name: " + name);
    }
    auto& node = graph.add_node(name);
    node.type = 'S';
}


void TransMap::add_path(const string& name){
    if (name == read_node_name or name == sample_node_name or name == path_node_name){
        throw runtime_error("ERROR: cannot add node with preset node name: " + name);
    }
    auto& node = graph.add_node(name);
    node.type = 'P';
}


void TransMap::add_edge(const string& a, const string& b){
    graph.add_edge(a,b,0);
}


void TransMap::add_edge(const string& a, const string& b, float weight){
    graph.add_edge(a,b,weight);
}


void TransMap::get_read_sample(const string& read_name, string& result){
    result.clear();

    // Using the HeteroGraph back end means that we iterate, even for a 1:1 mapping.
    graph.for_each_neighbor_of_type(read_name, 'S', [&](const HeteroNode& neighbor){
        // Check for cases that should be impossible
        if (not result.empty()){
            throw runtime_error("ERROR: multiple samples found for read: " + read_name);
        }

        result = neighbor.name;
    });

    if (result.empty()){
        throw runtime_error("ERROR: no sample found for read: " + read_name);
    }
}


void TransMap::for_each_read_of_sample(const string& sample_name, const function<void(const HeteroNode& node)>& f){
    graph.for_each_neighbor_of_type(sample_name, 'R', [&](const HeteroNode& neighbor){
        f(neighbor);
    });
}


void TransMap::for_each_sample_of_read(const string& read_name, const function<void(const HeteroNode& node)>& f){
    graph.for_each_neighbor_of_type(read_name, 'S', [&](const HeteroNode& neighbor){
        f(neighbor);
    });
}


void TransMap::for_each_sample_of_path(const string& path_name, const function<void(const HeteroNode& node)>& f){
    graph.for_each_neighbor_of_type(path_name, 'R', [&](const HeteroNode& r){
        graph.for_each_neighbor_of_type(r.name, 'S', [&](const HeteroNode& s){
            f(s);
        });
    });
}


void TransMap::for_each_path_of_sample(const string& sample_name, const function<void(const HeteroNode& node)>& f){
    graph.for_each_neighbor_of_type(sample_name, 'R', [&](const HeteroNode& r){
        graph.for_each_neighbor_of_type(r.name, 'P', [&](const HeteroNode& p){
            f(p);
        });
    });
}



}
