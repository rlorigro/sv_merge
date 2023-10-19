#pragma once

#include "Graph.hpp"

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

class TransMapNode: public Node {
public:
    string name;
    char type;

    explicit TransMapNode(const string& name);
    explicit TransMapNode(string& name);
};


TransMapNode::TransMapNode(const string& name):
        Node(name),
        type('*')
{}


TransMapNode::TransMapNode(string& name):
        Node(name),
        type('*')
{}


class TransMap {
    /// Attributes
    Graph<TransMapNode> graph;
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

    const TransMapNode& get_read_sample(const string& read_id);
    void for_each_read_of_sample(const string& sample_name, const function<void(const TransMapNode& node)>& f);
    void for_each_sample_of_read(const string& read_name, const function<void(const TransMapNode& node)>& f);
    void for_each_sample_of_path(const string& path_name, const function<void(const TransMapNode& node)>& f);
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
    node.type;
}


void TransMap::add_sample(const string& name){
    if (name == read_node_name or name == sample_node_name or name == path_node_name){
        throw runtime_error("ERROR: cannot add node with preset node name: " + name);
    }
    auto& node = graph.add_node(name);
    node.type;
}


void TransMap::add_path(const string& name){
    if (name == read_node_name or name == sample_node_name or name == path_node_name){
        throw runtime_error("ERROR: cannot add node with preset node name: " + name);
    }
    auto& node = graph.add_node(name);
    node.type;
}


void TransMap::add_edge(const string& a, const string& b){
    graph.add_edge(a,b,0);
}


void TransMap::add_edge(const string& a, const string& b, float weight){
    graph.add_edge(a,b,weight);
}


}
