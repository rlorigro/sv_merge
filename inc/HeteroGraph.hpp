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
 * The only requirements for the abstract HeteroNode class are that it has a string name and char type,
 * everything else is up to the user
 */
class HeteroNode {
public:
    string name;
    char type;

    explicit HeteroNode(const string& name);
    explicit HeteroNode(string& name);
    explicit HeteroNode(const string& name, char type);
    explicit HeteroNode(string& name, char type);
};


HeteroNode::HeteroNode(const string& name):
        name(name),
        type('*')
{}


HeteroNode::HeteroNode(string& name):
        name(name),
        type('*')
{}


HeteroNode::HeteroNode(const string& name, char type):
        name(name),
        type(type)
{}


HeteroNode::HeteroNode(string& name, char type):
        name(name),
        type(type)
{}


template<class T> class HeteroGraph {
    // node_id --> type_of_neighbor --> neighbor_id --> weight
    unordered_map<int64_t, unordered_map<char, unordered_map <int64_t, float> > > edges;
    unordered_map<int64_t, T> nodes;
    unordered_map<string, int64_t> id_map;
    int64_t id_counter = 0;

public:
    int64_t name_to_id(const string& name) const;

    /// Building
    T& add_node(const string& name, char type);

    void add_edge(const string& name_a, const string& name_b, float weight);
    void remove_edge(const string& name_a, const string& name_b);

    /// Accessing
    const T& get_node(const string& name) const;
    T& get_node(const string& name);
    const T& get_node(int64_t id) const;
    T& get_node(int64_t id);

    /// Global iterators
    void for_each_edge(const function<void(const string& a,const string& b, float weight)>& f) const;
    void for_node_in_bfs(const string& start_name, const function<void(const T& node, int64_t id)>& f) const;

    // This iterator allows the user to filter nodes out within the inner loop of BFS before returning them
    void for_node_in_bfs(
            const string& start_name,
            const function<bool(const T& node)>& criteria,
            const function<void(const T& node, int64_t id)>& f) const;

    /// Local iterators
    void for_each_neighbor(const string& name, const function<void(const T& neighbor, int64_t id)>& f) const;
    void for_each_neighbor_of_type(const string& name, char type, const function<void(const T& neighbor, int64_t id)>& f) const;
    void for_each_neighbor_of_type(int64_t id, char type, const function<void(int64_t)>& f) const;
    void for_each_neighbor_of_type(int64_t id, char type, const function<void(int64_t, float w)>& f) const;
};

template<class T> T& HeteroGraph<T>::add_node(const string& name, char type) {
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
    auto result = nodes.emplace(id, T(name,type));

    // Return a reference to the node
    return result.first->second;
}


template<class T> T& HeteroGraph<T>::get_node(const string& name) {
    auto id = name_to_id(name);
    return nodes.at(id);
}


template<class T> const T& HeteroGraph<T>::get_node(const string& name) const {
    auto id = name_to_id(name);
    return nodes.at(id);
}


template<class T> T& HeteroGraph<T>::get_node(int64_t id) {
    return nodes.at(id);
}


template<class T> const T& HeteroGraph<T>::get_node(int64_t id) const{
    return nodes.at(id);
}


template<class T> int64_t HeteroGraph<T>::name_to_id(const string& name) const{
    auto result = id_map.find(name);

    if (result == id_map.end()) {
        throw runtime_error("ERROR: cannot find name: " + name);
    } else {
        return result->second;
    }
}


template<class T> void HeteroGraph<T>::add_edge(const string& name_a, const string& name_b, float weight) {
    auto id_a = name_to_id(name_a);
    auto id_b = name_to_id(name_b);

    auto type_a = nodes.at(id_a).type;
    auto type_b = nodes.at(id_b).type;

    // Undirected by default, add both directions. Also, overwrite if it existed already.
    // Because this is a heterogeneous graph we also want to sort the edges by type of neighbor
    edges[id_a][type_b][id_b] = weight;
    edges[id_b][type_a][id_a] = weight;
}


template<class T> void HeteroGraph<T>::remove_edge(const string& name_a, const string& name_b) {
    auto id_a = name_to_id(name_a);
    auto id_b = name_to_id(name_b);

    auto type_a = nodes.at(id_a).type;
    auto type_b = nodes.at(id_b).type;

    // Undirected by default
    auto& result_a_type = edges.at(id_a);
    auto& result_a = result_a_type.at(type_b);
    result_a.erase(id_b);

    auto& result_b_type = edges.at(id_b);
    auto& result_b = result_b_type.at(type_a);
    result_b.erase(id_a);

    // If there are no more IDs corresponding to a given ID, then don't leave an empty container dangling
    if (result_a.empty()){
        result_a_type.erase(id_a);
    }
    if (result_b.empty()){
        result_b_type.erase(id_b);
    }

    // If there are no more IDs corresponding to a given ID, then don't leave an empty container dangling
    if (result_a_type.empty()){
        edges.erase(type_a);
    }
    if (result_b_type.empty()){
        edges.erase(type_b);
    }
}


template<class T> void HeteroGraph<T>::for_each_edge(const function<void(const string& a, const string& b, float weight)>& f) const{
    for (const auto& [id_a,item]: edges){
        // Iterate all types indiscriminately
        for (const auto& [type_b,item2]: item) {
            for (const auto &[id_b,w]: item2) {
                auto a = nodes.at(id_a);
                auto b = nodes.at(id_b);
                f(a.name, b.name, w);
            }
        }
    }
}


template<class T> void HeteroGraph<T>::for_each_neighbor(const string& name, const function<void(const T& neighbor, int64_t id)>& f) const{
    auto id = name_to_id(name);

    // Iterate all types indiscriminately
    for (const auto& [type_b,item]: edges.at(id)) {
        for (const auto &[id_b,w]: item) {
            f(nodes.at(id_b), id_b);
        }
    }
}


template<class T> void HeteroGraph<T>::for_each_neighbor_of_type(const string& name, char type, const function<void(const T& neighbor, int64_t id)>& f) const{
    auto id = name_to_id(name);

    // Iterate only the type of neighbor specified
    for (const auto& [id_b,w]: edges.at(id).at(type)) {
        f(nodes.at(id_b), id_b);
    }
}


template<class T> void HeteroGraph<T>::for_each_neighbor_of_type(int64_t id, char type, const function<void(int64_t)>& f) const{
    // Iterate only the type of neighbor specified
    for (const auto& [id_b,w]: edges.at(id).at(type)) {
        f(id_b);
    }
}


template<class T> void HeteroGraph<T>::for_each_neighbor_of_type(int64_t id, char type, const function<void(int64_t, float w)>& f) const{
    // Iterate only the type of neighbor specified
    for (const auto& [id_b,w]: edges.at(id).at(type)) {
        f(id_b, w);
    }
}


template<class T> void HeteroGraph<T>::for_node_in_bfs(const string& start_name, const function<void(const T& node, int64_t id)>& f) const{
    auto start_id = name_to_id(start_name);

    unordered_set<int64_t> visited;
    queue<int64_t> q;

    visited.emplace(start_id);
    q.emplace(start_id);

    while (not q.empty()){
        auto n = q.front();
        q.pop();

        f(nodes.at(n), n);

        // Iterate all types indiscriminately
        for (const auto& [type_b,item]: edges.at(n)) {
            for (const auto &[n_other,w]: item) {
                if (visited.find(n_other) == visited.end()) {
                    q.emplace(n_other);
                    visited.emplace(n_other);
                }
            }
        }
    }
}


template<class T> void HeteroGraph<T>::for_node_in_bfs(
        const string& start_name,
        const function<bool(const T& node)>& criteria,
        const function<void(const T& node, int64_t id)>& f) const{
    auto start_id = name_to_id(start_name);

    unordered_set<int64_t> visited;
    queue<int64_t> q;

    visited.emplace(start_id);
    q.emplace(start_id);

    while (not q.empty()){
        auto n = q.front();
        q.pop();

        f(nodes.at(n), n);

        // Iterate all types indiscriminately
        for (const auto& [type_b,item]: edges.at(n)) {
            for (const auto &[n_other,w]: item) {
                if (visited.find(n_other) == visited.end() and criteria(nodes.at(n_other)) == true) {
                    q.emplace(n_other);
                    visited.emplace(n_other);
                }
            }
        }
    }
}


}
