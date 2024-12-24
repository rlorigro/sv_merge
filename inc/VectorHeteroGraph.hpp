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
using std::to_string;
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

/// The HeteroGraph class is a general-purpose graph that provides an interface that is convenient for accessing nodes
/// by their type and name, and iterating over edges and nodes. Nodes can have different types, and adjacencies can be
/// filtered by type. Node IDs and Names must be unique.
template<class T> class HeteroGraph {
    unordered_map<int64_t, vector<pair <int64_t,float> > > edges;
    unordered_map<int64_t, T> nodes;
    unordered_map<string, int64_t> id_map;

    /**
     * The first position of each type block in each sorted adjacency list.
     */
    unordered_map<char, unordered_map<int64_t,int64_t>> first_of_type;

    int64_t id_counter = 0;

    /**
     * Tells whether `first_of_type` reflects the current state of the adjacency lists.
     */
    bool is_first_of_type_updated = false;

public:

    /// Building
    void reserve_nodes(size_t n);
    void reserve_edges(size_t n);
    T& add_node(const string& name, char type);

    // Simple helper function to search a vector and update it (since edges are stored in vectors)
    void update_adjacency_list(int64_t id, float weight, vector<pair <int64_t,float> >& adjacencies);

    /**
     * Sorts the neighbors of every node by `(type,id)`. This order is the same for every set of neighbors.
     */
    void sort_adjacency_lists();

    /**
     * Updates `first_of_type`, assuming that `sort_adjacency_lists()` has already been called.
     */
    void update_first_of_type();

    /**
     * @return TRUE iff no duplicated ID exists in any vector of `edges`.
     */
    bool are_edges_distinct() const;

    void add_edge(const string& name_a, const string& name_b, float weight);
    void add_edge(int64_t id_a, int64_t id_b, float weight);
    void remove_edge(const string& name_a, const string& name_b);
    void remove_edge(int64_t id_a, int64_t id_b);
    void remove_node(int64_t id);

    /// Accessing
    int64_t name_to_id(const string& name) const;
    pair<bool,int64_t> try_name_to_id(const string& name) const;

    const T& get_node(const string& name) const;
    vector<pair<int64_t,float>> get_edges(int64_t id) const;
    T& get_node(const string& name);
    const T& get_node(int64_t id) const;
    T& get_node(int64_t id);
    int64_t get_node_count() const;
    int64_t get_edge_count() const;
    int64_t get_edge_count(int64_t id) const;
    pair<bool,float> try_get_edge_weight(int64_t id_a, int64_t id_b) const;

    /**
     * Assigns `new_weight` to both `id_a->id_b` and `id_b->id_a`, if they exist.
     */
    void update_edge_weight(int64_t id_a, int64_t id_b, float new_weight);

    bool has_edge(int64_t id_a, int64_t id_b) const;
    bool has_node(const string& name) const;

    /// Global iterators
    void for_each_edge(const function<void(const string& a,const string& b, float weight)>& f) const;
    void for_each_node(const function<void(int64_t id, const T& node)>& f) const;
    void for_node_in_bfs(const string& start_name, const function<void(const T& node, int64_t id)>& f) const;

    // This iterator allows the user to filter nodes out within the inner loop of BFS before returning them
    void for_node_in_bfs(
            const string& start_name,
            float min_edge_weight,
            const function<bool(const T& node)>& criteria,
            const function<void(const T& node, int64_t id)>& f) const;

    /// Local iterators
    void for_each_neighbor(const string& name, const function<void(const T& neighbor, int64_t id)>& f) const;
    void for_each_neighbor(int64_t n, const function<void(const T& neighbor, int64_t id)>& f) const;
    void for_each_neighbor_of_type(const string& name, char type, const function<void(const T& neighbor, int64_t id)>& f) const;
    void for_each_neighbor_of_type(int64_t id, char type, const function<void(int64_t)>& f) const;
    void for_each_neighbor_of_type(int64_t id, char type, const function<void(int64_t, float w)>& f) const;

    /// Helpers
    pair<int64_t,int64_t> canonicalize_edge(int64_t a, int64_t b) const;
};


template<class T> void HeteroGraph<T>::reserve_nodes(size_t n){
    id_map.reserve(n);
    nodes.reserve(n);
}


template<class T> void HeteroGraph<T>::reserve_edges(size_t n){
    edges.reserve(n);
}


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


template<class T> int64_t HeteroGraph<T>::get_node_count() const{
    return nodes.size();
}


template<class T> int64_t HeteroGraph<T>::get_edge_count() const {
    int64_t out = 0;
    for (const auto& [id_a,item]: edges) { out+=(int64_t)item.size(); }
    return out/2;
}


template<class T> const T& HeteroGraph<T>::get_node(int64_t id) const{
    return nodes.at(id);
}


template<class T> vector<pair<int64_t,float>> HeteroGraph<T>::get_edges(int64_t node_id) const {
    return edges.at(node_id);
}


template<class T> int64_t HeteroGraph<T>::name_to_id(const string& name) const{
    auto result = id_map.find(name);

    if (result == id_map.end()) {
        throw runtime_error("ERROR: cannot find name: " + name);
    } else {
        return result->second;
    }
}


template<class T> pair<bool,int64_t> HeteroGraph<T>::try_name_to_id(const string& name) const{
    auto result = id_map.find(name);

    if (result == id_map.end()) {
        return {false,-1};
    } else {
        return {true,result->second};
    }
}


template<class T> pair<int64_t,int64_t> HeteroGraph<T>::canonicalize_edge(int64_t a, int64_t b) const{
    if (a > b){
        return {b,a};
    }
    else{
        return {a,b};
    }
}


template<class T> void HeteroGraph<T>::update_adjacency_list(
        int64_t id,
        float weight,
        vector<pair <int64_t,float> >& adjacencies) {

    // Search for the destination node id in the adjacency list
    auto ab = std::find_if(adjacencies.begin(), adjacencies.end(), [&](pair<int64_t, float>& p){
        return p.first == id;
    });

    // Only append the vector if id is not found (edge doesn't already exist)
    if (ab == adjacencies.end()){
        adjacencies.emplace_back(id, weight);
        is_first_of_type_updated=false;
    }
    // Otherwise just overwrite the weight
    else{
        ab->second = weight;
    }
}


template<class T> void HeteroGraph<T>::sort_adjacency_lists() {
    for (auto& element: edges) {
        std::sort(element.second.begin(), element.second.end(), [&](const pair <int64_t,float>& a, const pair <int64_t,float>& b) {
            return (nodes.at(a.first).type<nodes.at(b.first).type || (nodes.at(a.first).type==nodes.at(b.first).type && a.first<b.first));
        });
    }
}


template<class T> void HeteroGraph<T>::update_first_of_type() {
    char type, current_type;
    size_t i;
    size_t length;

    if (is_first_of_type_updated) return;
    first_of_type.clear();
    for (auto& element: edges) {
        length=element.second.size();
        type='_';
        for (i=0; i<length; i++) {
            current_type=nodes.at(element.second.at(i).first).type;
            if (current_type!=type) {
                if (first_of_type.contains(current_type)) first_of_type.at(current_type).emplace(element.first,(int64_t)i);
                else {
                    unordered_map<int64_t,int64_t> new_map;
                    new_map.emplace(element.first,(int64_t)i);
                    first_of_type.emplace(current_type,new_map);
                }
                type=current_type;
            }
        }
    }
    is_first_of_type_updated=true;
}


template<class T> void HeteroGraph<T>::add_edge(const string& name_a, const string& name_b, float weight) {
    auto id_a = name_to_id(name_a);
    auto id_b = name_to_id(name_b);

    auto& result_a = edges[id_a];
    auto& result_b = edges[id_b];

    update_adjacency_list(id_b, weight, result_a);
    update_adjacency_list(id_a, weight, result_b);
}


// TODO: Untested
template<class T> void HeteroGraph<T>::add_edge(int64_t id_a, int64_t id_b, float weight) {
    auto result_a = nodes.find(id_a);
    auto result_b = nodes.find(id_b);

    if (result_a == nodes.end()){
        throw runtime_error("ERROR: cannot find id: " + to_string(id_a));
    }

    if (result_b == nodes.end()){
        throw runtime_error("ERROR: cannot find id: " + to_string(id_b));
    }

    if ((result_a->second.type == 'V' and result_b->second.type != 'P') or (result_a->second.type != 'P' and result_b->second.type == 'V')){
        throw runtime_error("ERROR: cannot add edge from variant to non-path: " + to_string(id_a) + " " + to_string(id_b));
    }

    auto& adjacencies_a = edges[id_a];
    auto& adjacencies_b = edges[id_b];

    update_adjacency_list(id_b, weight, adjacencies_a);
    update_adjacency_list(id_a, weight, adjacencies_b);
}


template<class T> bool HeteroGraph<T>::has_edge(int64_t id_a, int64_t id_b) const {
    return try_get_edge_weight(id_a, id_b).first;
}


template<class T> bool HeteroGraph<T>::has_node(const string& name) const {
    return id_map.find(name) != id_map.end();
}


// TODO: Untested
template<class T> pair<bool,float> HeteroGraph<T>::try_get_edge_weight(int64_t id_a, int64_t id_b) const{
    pair<bool,float> value = {false, 0};

    auto result = edges.find(id_a);

    if (result != edges.end()){
        auto& adjacencies = result->second;

        auto ab = std::find_if(adjacencies.begin(), adjacencies.end(), [&](const pair<int64_t, float>& p){
            return p.first == id_b;
        });

        if (ab != adjacencies.end()) {
            value = {true, ab->second};
        }
    }

    return value;
}


template<class T> void HeteroGraph<T>::update_edge_weight(int64_t id_a, int64_t id_b, float new_weight) {
    // Updating `a->b`, if it exists.
    auto result = edges.find(id_a);
    if (result!=edges.end()) {
        auto& adjacencies = result->second;
        auto ab = std::find_if(adjacencies.begin(),adjacencies.end(),[&](const pair<int64_t, float>& p){ return p.first==id_b; });
        if (ab!=adjacencies.end()) (*ab).second=new_weight;
    }

    // Updating `b->a`, if it exists.
    result=edges.find(id_b);
    if (result!=edges.end()) {
        auto& adjacencies = result->second;
        auto ab = std::find_if(adjacencies.begin(),adjacencies.end(),[&](const pair<int64_t, float>& p){ return p.first==id_a; });
        if (ab!=adjacencies.end()) (*ab).second=new_weight;
    }
}


template<class T> void HeteroGraph<T>::remove_edge(const string& name_a, const string& name_b) {
    auto id_a = name_to_id(name_a);
    auto id_b = name_to_id(name_b);

    remove_edge(id_a, id_b);
}


template<class T> void HeteroGraph<T>::remove_node(int64_t id){
    vector <pair <int64_t, int64_t> > edges_to_remove;

    for_each_neighbor(id, [&](const T& neighbor, int64_t id_b){
        edges_to_remove.emplace_back(id, id_b);
    });

    for (const auto& edge: edges_to_remove){
        remove_edge(edge.first, edge.second);
    }

    id_map.erase(nodes.at(id).name);
    nodes.erase(id);
}


template<class T> void HeteroGraph<T>::remove_edge(int64_t id_a, int64_t id_b) {
    auto result_a = edges.find(id_a);

    if (result_a != edges.end()){
        auto& adjacencies = result_a->second;

        auto ab = std::find_if(adjacencies.begin(), adjacencies.end(), [&](const pair<int64_t, float>& p){
            return p.first == id_b;
        });

        if (ab != adjacencies.end()){
            adjacencies.erase(ab);
            is_first_of_type_updated=false;
        }
    }

    auto result_b = edges.find(id_b);

    if (result_b != edges.end()){
        auto& adjacencies = result_b->second;

        auto ba = std::find_if(adjacencies.begin(), adjacencies.end(), [&](const pair<int64_t, float>& p){
            return p.first == id_a;
        });

        if (ba != adjacencies.end()) {
            adjacencies.erase(ba);
            is_first_of_type_updated=false;
        }
    }
}


template<class T> void HeteroGraph<T>::for_each_edge(const function<void(const string& a, const string& b, float weight)>& f) const{
    for (const auto& [id_a,item]: edges){
        // Iterate all types indiscriminately
        for (const auto& [id_b,w]: item) {
            auto a = nodes.at(id_a);
            auto b = nodes.at(id_b);
            f(a.name, b.name, w);
        }
    }
}


template<class T> void HeteroGraph<T>::for_each_node(const function<void(int64_t id, const T& node)>& f) const{
    for (const auto& [id,node]: nodes){
        f(id, node);
    }
}


template<class T> void HeteroGraph<T>::for_each_neighbor(const string& name, const function<void(const T& neighbor, int64_t id)>& f) const{
    auto id = name_to_id(name);

    auto result = edges.find(id);

    // No neighbors
    if (result == edges.end()){
        return;
    }

    // Iterate all types indiscriminately
    for (const auto& [id_b,w]: result->second) {
        f(nodes.at(id_b), id_b);
    }
}


template<class T> void HeteroGraph<T>::for_each_neighbor(int64_t n, const function<void(const T& neighbor, int64_t id)>& f) const{
    auto result = edges.find(n);

    // Iterate all types indiscriminately
    for (const auto& [id_b,w]: result->second) {
        f(nodes.at(id_b), id_b);
    }
}


template<class T> int64_t HeteroGraph<T>::get_edge_count(int64_t id) const{
    auto result = edges.find(id);

    // No neighbors
    if (result == edges.end()){
        return 0;
    }

    return int64_t(result->second.size());
}


template<class T> void HeteroGraph<T>::for_each_neighbor_of_type(const string& name, char type, const function<void(const T& neighbor, int64_t id)>& f) const{
    auto id = name_to_id(name);

    // Check if there are any edges from this node
    const auto result = edges.find(id);
    if (result==edges.end()) return;

    // Iterate all edges, but only operate on the specified type
    if (is_first_of_type_updated) {
        int64_t first = first_of_type.at(type).at(id);
        const size_t length = result->second.size();
        for (size_t i=first; i<length; i++) {
            const int64_t id_b = result->second.at(i).first;
            const auto &node = nodes.at(id_b);
            if (node.type!=type) break;
            f(node,id_b);
        }
    }
    else {
        for (const auto &[id_b, w]: result->second) {
            const auto &node = nodes.at(id_b);
            if (node.type==type) f(node,id_b);
        }
    }
}


template<class T> void HeteroGraph<T>::for_each_neighbor_of_type(int64_t id, char type, const function<void(int64_t)>& f) const{
    // Check if there are any edges from this node
    const auto result = edges.find(id);
    if (result==edges.end()) return;

    // Iterate all edges, but only operate on the specified type
    if (is_first_of_type_updated) {
        int64_t first = first_of_type.at(type).at(id);
        const size_t length = result->second.size();
        for (size_t i=first; i<length; i++) {
            const int64_t id_b = result->second.at(i).first;
            if (nodes.at(id_b).type!=type) break;
            f(id_b);
        }
    }
    else {
        for (const auto& [id_b,w]: result->second) {
            if (nodes.at(id_b).type==type) f(id_b);
        }
    }
}


template<class T> void HeteroGraph<T>::for_each_neighbor_of_type(int64_t id, char type, const function<void(int64_t, float w)>& f) const{
    // Check if there are any edges from this node
    const auto result = edges.find(id);
    if (result==edges.end()) return;

    // Iterate all edges, but only operate on the specified type
    if (is_first_of_type_updated) {
        int64_t first = first_of_type.at(type).at(id);
        const size_t length = result->second.size();
        for (size_t i=first; i<length; i++) {
            const int64_t id_b = result->second.at(i).first;
            if (nodes.at(id_b).type!=type) break;
            f(id_b,result->second.at(i).second);
        }
    }
    else {
        for (const auto &[id_b, w]: result->second) {
            if (nodes.at(id_b).type==type) f(id_b,w);
        }
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
        for (const auto& [n_other,w]: edges.at(n)) {
            if (visited.find(n_other) == visited.end()) {
                q.emplace(n_other);
                visited.emplace(n_other);
            }
        }
    }
}


template<class T> void HeteroGraph<T>::for_node_in_bfs(
        const string& start_name,
        float min_edge_weight,
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
        for (const auto& [n_other,w]: edges.at(n)) {
            if (w < min_edge_weight){
                continue;
            }

            if (visited.find(n_other) == visited.end() and criteria(nodes.at(n_other)) == true) {
                q.emplace(n_other);
                visited.emplace(n_other);
            }
        }
    }
}


template<class T> bool HeteroGraph<T>::are_edges_distinct() const {
    size_t i, j;
    for (const auto& [id_a,item]: edges){
        for (i=0; i<item.size(); i++) {
            for (j=i+1; j<item.size(); j++) {
                if (item.at(i).first==item.at(j).first) {
                    cerr << "are_edges_distinct> DUPLICATED NEIGHBOR OF id_a=" << to_string(id_a) << ": i=" << to_string(i) << " (" << to_string(item.at(i).first) << "," << to_string(item.at(i).second) << ") " << ", j=" << to_string(j) << " (" << to_string(item.at(j).first) << "," << to_string(item.at(j).second) << ")\n";
                    return false;
                }
            }
        }
    }
    return true;
}


}
