#pragma once

#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <stdexcept>
#include <iostream>
#include <utility>
#include <string>
#include <vector>
#include <queue>

using std::unordered_map;
using std::unordered_set;
using std::runtime_error;
using std::function;
using std::vector;
using std::string;
using std::queue;
using std::pair;
using std::cerr;


namespace sv_merge {


/**
 * The only requirement for the abstract Node class is that it has a string name, everything else is up to the user
 */
class Node {
public:
    string name;

    explicit Node(const string& name);
    explicit Node(string& name);
};


Node::Node(const string& name):
    name(name)
{}


Node::Node(string& name):
    name(name)
{}


template<class T> class Graph {
    unordered_map<int64_t, unordered_map <int64_t, float> > edges;
    unordered_map<int64_t, T > nodes;
    unordered_map<string, int64_t> id_map;
    int64_t id_counter = 0;

public:
    int64_t name_to_id(const string& name) const;

    T& add_node(const string& name);

    const T& get_node(const string& name) const;
    T& get_node(const string& name);
    T& get_node(int64_t id);

    void add_edge(const string& name_a, const string& name_b, float weight);
    void remove_edge(const string& name_a, const string& name_b);

    /// Global iterators
    void for_each_edge(const function<void(const string& a,const string& b, float weight)>& f) const;
    void for_node_in_bfs(const string& start_name, const function<void(const T& node)>& f) const;
    void for_node_in_bfs(
            const string& start_name,
            const function<bool(const T& node)>& criteria,
            const function<void(const T& node)>& f) const;

    /// Local iterators
    void for_each_neighbor(const string& name, const function<void(const T& neighbor)>& f) const;
};


template<class T> T& Graph<T>::add_node(const string& name) {
    if (id_map.find(name) != id_map.end()){
        throw runtime_error("ERROR: cannot add node with existing name: " + name);
    }

    // The id can be anything as long as it doesn't exist yet, because a map is used instead of a vector here.
    // The simplest way not to reuse an old ID is just incrementing from zero.
    auto id = id_counter;
    id_counter++;

    // Track the reverse mapping from name -> id;
    id_map.emplace(name,id);

    // Place a new node in the graph
    auto result = nodes.emplace(id, T(name));

    // Return a reference to the node
    return result.first->second;
}


template<class T> T& Graph<T>::get_node(const string& name) {
    auto id = name_to_id(name);
    return nodes.at(id);
}


template<class T> const T& Graph<T>::get_node(const string& name) const {
    auto id = name_to_id(name);
    return nodes.at(id);
}


template<class T> T& Graph<T>::get_node(int64_t id) {
    return nodes.at(id);
}


template<class T> int64_t Graph<T>::name_to_id(const string& name) const{
    auto result = id_map.find(name);

    if (result == id_map.end()) {
        throw runtime_error("ERROR: cannot find name: " + name);
    } else {
        return result->second;
    }
}


template<class T> void Graph<T>::add_edge(const string& name_a, const string& name_b, float weight) {
    auto id_a = name_to_id(name_a);
    auto id_b = name_to_id(name_b);

    // Undirected by default, add both directions. Also, don't care if it existed already
    edges[id_a][id_b] = weight;
    edges[id_b][id_a] = weight;
}


template<class T> void Graph<T>::remove_edge(const string& name_a, const string& name_b) {
    auto id_a = name_to_id(name_a);
    auto id_b = name_to_id(name_b);

    // Undirected by default
    auto& result_a = edges.at(id_a);
    result_a.erase(id_b);

    auto& result_b = edges.at(id_b);
    result_b.erase(id_a);

    // If there are no more IDs corresponding to a given ID, then don't leave an empty container dangling
    if (result_a.empty()){
        edges.erase(id_a);
    }
    if (result_b.empty()){
        edges.erase(id_b);
    }
}


template<class T> void Graph<T>::for_each_edge(const function<void(const string& a, const string& b, float weight)>& f) const{
    for (const auto& [id_a,item]: edges){
        for (const auto& [id_b, w]: item){
            auto a = nodes.at(id_a);
            auto b = nodes.at(id_b);
            f(a.name, b.name, w);
        }
    }
}


template<class T> void Graph<T>::for_each_neighbor(const string& name, const function<void(const T& neighbor)>& f) const{
    auto id = name_to_id(name);
    for (const auto& [id_b, w]: edges.at(id)){
        auto b = nodes.at(id_b);
        f(b);
    }
}


template<class T> void Graph<T>::for_node_in_bfs(const string& start_name, const function<void(const T& node)>& f) const{
    auto start_id = name_to_id(start_name);

    unordered_set<int64_t> visited;
    queue<int64_t> q;

    visited.emplace(start_id);
    q.emplace(start_id);

    while (not q.empty()){
        auto n = q.front();
        q.pop();

        f(nodes.at(n));

        for (auto& [n_other,weight]: edges.at(n)){
            if (visited.find(n_other) == visited.end()){
                q.emplace(n_other);
                visited.emplace(n_other);
            }
        }
    }
}


template<class T> void Graph<T>::for_node_in_bfs(
        const string& start_name,
        const function<bool(const T& node)>& criteria,
        const function<void(const T& node)>& f) const{
    auto start_id = name_to_id(start_name);

    unordered_set<int64_t> visited;
    queue<int64_t> q;

    visited.emplace(start_id);
    q.emplace(start_id);

    while (not q.empty()){
        auto n = q.front();
        q.pop();

        f(nodes.at(n));

        for (auto& [n_other,weight]: edges.at(n)){
            if (visited.find(n_other) == visited.end() and criteria(nodes.at(n_other)) == true){
                q.emplace(n_other);
                visited.emplace(n_other);
            }
        }
    }
}


}
