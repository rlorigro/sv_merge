#include "TransitiveMap.hpp"
#include <fstream>

using std::ofstream;


namespace sv_merge{

TransMap::TransMap():
        sample_node_name("sample_node"),
        read_node_name("read_node"),
        path_node_name("path_node")
{
    // Source nodes use lower case types to avoid being confused with the types they point to
    graph.add_node(sample_node_name, 's');
    graph.add_node(read_node_name, 'r');
    graph.add_node(path_node_name, 'p');
}


void TransMap::reserve_nodes(size_t n){
    graph.reserve_nodes(n);
}


void TransMap::reserve_edges(size_t n){
    graph.reserve_edges(n);
}


void TransMap::reserve_sequences(size_t n){
    sequences.reserve(n);
}


int64_t TransMap::get_id(const string& name) const{
    return graph.name_to_id(name);
}


pair<bool,int64_t> TransMap::try_get_id(const string& name) const{
    return graph.try_name_to_id(name);
}


int64_t TransMap::get_node_count() const{
    return graph.get_node_count();
}


int64_t TransMap::get_edge_count(int64_t id, char type) const{
    return graph.get_edge_count(id, type);
}


const HeteroNode& TransMap::get_node(int64_t id) const{
    return graph.get_node(id);
}


const HeteroNode& TransMap::get_node(const string& name) const{
    auto id = graph.name_to_id(name);
    return graph.get_node(id);
}


const string& TransMap::get_sequence(const string& name) const{
    auto id = graph.name_to_id(name);
    return sequences.at(id);
}


const string& TransMap::get_sequence(int64_t id) const{
    return sequences.at(id);
}


void TransMap::add_flank_coord(const string& name, int32_t start, int32_t stop){
    auto id = graph.name_to_id(name);
    sequence_flanks[id] = {start,stop};
}


coord_t TransMap::get_flank_coord(const string& name) const{
    auto id = graph.name_to_id(name);
    return sequence_flanks.at(id);
}


coord_t TransMap::get_flank_coord(int64_t id) const{
    return sequence_flanks.at(id);
}


const unordered_map<int64_t,interval_t>& TransMap::get_flank_map() const {
    return sequence_flanks;
}


void TransMap::construct_named_flank_map(unordered_map<string,interval_t>& flank_map) const {
    flank_map.clear();
    flank_map.reserve(sequence_flanks.size());
    for (const auto& [id,item]: sequence_flanks){
        auto& name = get_node(id).name;
        flank_map.emplace(name, item);
    }
}


void TransMap::add_read(const string& name){
    if (name == read_node_name or name == sample_node_name or name == path_node_name){
        throw runtime_error("ERROR: cannot add node with preset node name: " + name);
    }
    graph.add_node(name, 'R');
    graph.add_edge(read_node_name, name, 0);
}


void TransMap::add_read(const string& name, const string& sequence){
    if (name == read_node_name or name == sample_node_name or name == path_node_name){
        throw runtime_error("ERROR: cannot add node with preset node name: " + name);
    }
    graph.add_node(name, 'R');
    graph.add_edge(read_node_name, name, 0);

    sequences.emplace(graph.name_to_id(name), sequence);
}


void TransMap::add_read_with_move(string& name, string& sequence){
    if (name == read_node_name or name == sample_node_name or name == path_node_name){
        throw runtime_error("ERROR: cannot add node with preset node name: " + name);
    }

    graph.add_node(name, 'R');
    graph.add_edge(read_node_name, name, 0);

    sequences.emplace(graph.name_to_id(name), std::move(sequence));
}


void TransMap::add_read_with_move(string& name, string& sequence, bool is_reverse){
    if (name == read_node_name or name == sample_node_name or name == path_node_name){
        throw runtime_error("ERROR: cannot add node with preset node name: " + name);
    }

    graph.add_node(name, 'R');
    graph.add_edge(read_node_name, name, 0);

    auto id = graph.name_to_id(name);
    sequences.emplace(id, std::move(sequence));
    sequence_reversals.emplace(id, is_reverse);
}


void TransMap::add_sample(const string& name){
    if (name == read_node_name or name == sample_node_name or name == path_node_name){
        throw runtime_error("ERROR: cannot add node with preset node name: " + name);
    }
    graph.add_node(name, 'S');
    graph.add_edge(sample_node_name, name, 0);
}


void TransMap::add_path(const string& name){
    if (name == read_node_name or name == sample_node_name or name == path_node_name){
        throw runtime_error("ERROR: cannot add node with preset node name: " + name);
    }
    graph.add_node(name, 'P');
    graph.add_edge(path_node_name, name, 0);
}


void TransMap::add_edge(const string& a, const string& b){
    graph.add_edge(a,b,0);
}


void TransMap::add_edge(int64_t a, int64_t b, float weight){
    graph.add_edge(a,b,weight);
}


pair<bool,float> TransMap::try_get_edge_weight(int64_t id_a, int64_t id_b) const{
    return graph.try_get_edge_weight(id_a, id_b);
}


bool TransMap::has_edge(int64_t a, int64_t b) const{
    return graph.has_edge(a,b);
}


void TransMap::remove_edge(int64_t a, int64_t b){
    return graph.remove_edge(a,b);
}


void TransMap::add_edge(const string& a, const string& b, float weight){
    graph.add_edge(a,b,weight);
}


void TransMap::get_read_sample(const string& read_name, string& result) const{
    result.clear();

    // Using the HeteroGraph back end means that we iterate, even for a 1:1 mapping.
    graph.for_each_neighbor_of_type(read_name, 'S', [&](const HeteroNode& neighbor, int64_t id){
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


void TransMap::get_read_sample(int64_t read_id, string& result) const{
    result.clear();

    if (graph.get_node(read_id).type != 'R'){
        throw runtime_error("ERROR: non-read ID provided for get_read_sample: " + to_string(read_id) + " " + graph.get_node(read_id).name);
    }

    // Using the HeteroGraph back end means that we iterate, even for a 1:1 mapping.
    graph.for_each_neighbor_of_type(read_id, 'S', [&](int64_t id){
        // Check for cases that should be impossible
        if (not result.empty()){
            throw runtime_error("ERROR: multiple samples found for id: " + to_string(read_id));
        }

        result = graph.get_node(id).name;
    });

    // For this project, every read must have exactly one sample.
    if (result.empty()){
        throw runtime_error("ERROR: no sample found for id: " + to_string(read_id));
    }
}


void TransMap::for_each_read(const function<void(const string& name, const string& sequence)>& f) const{
    graph.for_each_neighbor_of_type(read_node_name, 'R', [&](const HeteroNode& neighbor, int64_t id){
        f(neighbor.name, sequences.at(id));
    });
}


void TransMap::for_each_read(const function<void(const string& name, int64_t id)>& f) const{
    graph.for_each_neighbor_of_type(read_node_name, 'R', [&](const HeteroNode& neighbor, int64_t id){
        f(neighbor.name, id);
    });
}


void TransMap::for_each_read_id(const function<void(int64_t id)>& f) const{
    graph.for_each_neighbor_of_type(read_node_name, 'R', [&](const HeteroNode& neighbor, int64_t id){
        f(id);
    });
}


void TransMap::for_each_sample(const function<void(const string& name, int64_t id)>& f) const{
    graph.for_each_neighbor_of_type(sample_node_name, 'S', [&](const HeteroNode& neighbor, int64_t id){
        f(neighbor.name, id);
    });
}


void TransMap::for_each_neighbor_of_type(int64_t id, char type, const function<void(int64_t id)>& f) const{
    graph.for_each_neighbor_of_type(id, type, [&](int64_t id){
        f(id);
    });
}


void TransMap::for_each_path(const function<void(const string& name, int64_t id)>& f) const{
    graph.for_each_neighbor_of_type(path_node_name, 'P', [&](const HeteroNode& neighbor, int64_t id){
        f(neighbor.name, id);
    });
}


void TransMap::for_each_read_of_sample(const string& sample_name, const function<void(const string& name, int64_t id)>& f) const{
    graph.for_each_neighbor_of_type(sample_name, 'R', [&](const HeteroNode& neighbor, int64_t id){
        f(neighbor.name, id);
    });
}


void TransMap::for_each_read_of_sample(int64_t sample_id, const function<void(int64_t read_id)>& f) const{
    graph.for_each_neighbor_of_type(sample_id, 'R', [&](int64_t id){
        f(id);
    });
}


void TransMap::for_each_path_of_read(int64_t read_id, const function<void(int64_t path_id)>& f) const{
    graph.for_each_neighbor_of_type(read_id, 'P', [&](int64_t id){
        f(id);
    });
}


void TransMap::for_each_read_of_path(int64_t path_id, const function<void(int64_t id)>& f) const{
    graph.for_each_neighbor_of_type(path_id, 'R', [&](int64_t r){
        f(r);
    });
}


void TransMap::for_each_sample_of_read(const string& read_name, const function<void(const string& name, int64_t id)>& f) const{
    graph.for_each_neighbor_of_type(read_name, 'S', [&](const HeteroNode& neighbor, int64_t id){
        f(neighbor.name, id);
    });
}


void TransMap::for_each_sample_of_path(const string& path_name, const function<void(const string& name, int64_t id)>& f) const{
    auto id = graph.name_to_id(path_name);
    graph.for_each_neighbor_of_type(id, 'R', [&](int64_t r){
        graph.for_each_neighbor_of_type(r, 'S', [&](int64_t s){
            f(graph.get_node(s).name, s);
        });
    });
}


void TransMap::for_each_path_of_sample(const string& sample_name, const function<void(const string& name, int64_t id)>& f) const{
    auto id = graph.name_to_id(sample_name);
    graph.for_each_neighbor_of_type(id, 'R', [&](int64_t r){
        graph.for_each_neighbor_of_type(r, 'P', [&](int64_t p){
            f(graph.get_node(p).name, p);
        });
    });
}


void TransMap::for_each_read_to_path_edge(const function<void(int64_t read_id, int64_t path_id, float weight)>& f) const{
    auto id = graph.name_to_id(read_node_name);

    // Starting from the source node which connects to all reads, find all neighbors (read nodes)
    graph.for_each_neighbor_of_type(id, 'R', [&](int64_t read_id){
        // Find all path-type neighbors
        graph.for_each_neighbor_of_type(read_id, 'P', [&](int64_t path_id, float w){
            f(read_id, path_id, w);
        });
    });
}


void TransMap::for_node_in_bfs(
        const string& start_name,
        float min_edge_weight,
        const function<bool(const HeteroNode& node)>& criteria,
        const function<void(const HeteroNode& node, int64_t id)>& f) const{
    graph.for_node_in_bfs(start_name,min_edge_weight,criteria,f);
}


void TransMap::write_edge_info_to_csv(path output_path) const{
    ofstream out(output_path);

    out << "sample,read,path,weight" << '\n';

    for_each_sample([&](const string& sample_name, int64_t sample_id){
        for_each_read_of_sample(sample_id, [&](int64_t read_id){
            for_each_path_of_read(read_id, [&](int64_t path_id){
                auto [success, weight] = try_get_edge_weight(read_id, path_id);

                if (not success){
                    throw runtime_error("ERROR: edge weight not found for read-path edge: " + to_string(read_id) + " " + to_string(path_id));
                }

                out << sample_name << ',' << get_node(read_id).name << ',' << get_node(path_id).name << ',' << weight << '\n';
            });
        });
    });
}


}
