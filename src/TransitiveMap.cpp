#include "TransitiveMap.hpp"
#include "gaf.hpp"
#include <fstream>

using std::ofstream;
using std::tuple;


namespace sv_merge{

TransMap::TransMap(){
    // Source nodes use lower case types to avoid being confused with the types they point to
    graph.add_node(sample_node_name, 's');
    graph.add_node(read_node_name, 'r');
    graph.add_node(path_node_name, 'p');
    graph.add_node(variant_node_name, 'v');
}

const string TransMap::sample_node_name = "sample_node";
const string TransMap::read_node_name = "read_node";
const string TransMap::path_node_name = "path_node";
const string TransMap::variant_node_name = "variant_node";

void TransMap::reserve_nodes(size_t n){
    graph.reserve_nodes(n);
}


void TransMap::reserve_edges(size_t n){
    graph.reserve_edges(n);
}


void TransMap::reserve_sequences(size_t n){
    sequences.reserve(n);
}


bool TransMap::empty() const{
    // In the context of a TransMap, empty means that the graph has only the source/sink nodes
    return graph.get_node_count() == 4 and
        graph.get_edge_count(0) == 0 and
        graph.get_edge_count(1) == 0 and
        graph.get_edge_count(2) == 0 and
        graph.get_edge_count(3) == 0;
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


int64_t TransMap::get_edge_count(int64_t id) const{
    return graph.get_edge_count(id);
}


const HeteroNode& TransMap::get_node(int64_t id) const{
    return graph.get_node(id);
}


const HeteroNode& TransMap::get_node(const string& name) const{
    auto id = graph.name_to_id(name);
    return graph.get_node(id);
}


void TransMap::get_sequence(const string& name, string& result) const{
    auto id = graph.name_to_id(name);
    sequences.at(id).to_string(result);
}


void TransMap::get_sequence(int64_t id, string& result) const{
    sequences.at(id).to_string(result);
}


size_t TransMap::get_sequence_size(int64_t id) const{
    return sequences.at(id).size();
}


size_t TransMap::get_sequence_size(const string& name) const{
    auto id = graph.name_to_id(name);
    return sequences.at(id).size();
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
    if (name == read_node_name or name == sample_node_name or name == path_node_name or name == variant_node_name){
        throw runtime_error("ERROR: cannot add node with preset node name: " + name);
    }
    graph.add_node(name, 'R');
    graph.add_edge(read_node_name, name, 0);
}


void TransMap::add_read(const string& name, const string& sequence){
    if (name == read_node_name or name == sample_node_name or name == path_node_name or name == variant_node_name){
        throw runtime_error("ERROR: cannot add node with preset node name: " + name);
    }
    graph.add_node(name, 'R');
    graph.add_edge(read_node_name, name, 0);

    sequences.emplace(graph.name_to_id(name), sequence);
}


void TransMap::add_read_with_move(string& name, BinarySequence<uint64_t>& sequence){
    if (name == read_node_name or name == sample_node_name or name == path_node_name or name == variant_node_name){
        throw runtime_error("ERROR: cannot add node with preset node name: " + name);
    }

    graph.add_node(name, 'R');
    graph.add_edge(read_node_name, name, 0);

    sequences.emplace(graph.name_to_id(name), std::move(sequence));
}


void TransMap::add_read_with_move(string& name, BinarySequence<uint64_t>& sequence, bool is_reverse){
    if (name == read_node_name or name == sample_node_name or name == path_node_name or name == variant_node_name){
        throw runtime_error("ERROR: cannot add node with preset node name: " + name);
    }

    graph.add_node(name, 'R');
    graph.add_edge(read_node_name, name, 0);

    auto id = graph.name_to_id(name);
    sequences.emplace(id, std::move(sequence));
    sequence_reversals.emplace(id, is_reverse);
}


void TransMap::add_sample(const string& name){
    if (name == read_node_name or name == sample_node_name or name == path_node_name or name == variant_node_name){
        throw runtime_error("ERROR: cannot add node with preset node name: " + name);
    }
    graph.add_node(name, 'S');
    graph.add_edge(sample_node_name, name, 0);
}


void TransMap::add_path(const string& name){
    if (name == read_node_name or name == sample_node_name or name == path_node_name or name == variant_node_name){
        throw runtime_error("ERROR: cannot add node with preset node name: " + name);
    }
    graph.add_node(name, 'P');
    graph.add_edge(path_node_name, name, 0);
}


void TransMap::add_variant(const string& name){
    if (name == read_node_name or name == sample_node_name or name == path_node_name or name == variant_node_name){
        throw runtime_error("ERROR: cannot add node with preset node name: " + name);
    }
    graph.add_node(name, 'V');
    graph.add_edge(variant_node_name, name, 0);
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


bool TransMap::has_node(const string& name) const{
    return graph.has_node(name);
}


bool TransMap::has_node(int64_t id) const{
    return graph.has_node(id);
}


void TransMap::remove_edge(int64_t a, int64_t b){
    graph.remove_edge(a,b);
}


void TransMap::remove_node(int64_t id){
    // Additionally erase the sequence if it was a Read type node
    auto result = sequences.find(id);

    if (result != sequences.end()) {
        sequences.erase(result);
    }

    graph.remove_node(id);
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




int64_t TransMap::get_read_sample(int64_t read_id) const{
    if (graph.get_node(read_id).type != 'R'){
        throw runtime_error("ERROR: non-read ID provided for get_read_sample: " + to_string(read_id) + " " + graph.get_node(read_id).name);
    }

    int64_t result = -1;

    graph.for_each_neighbor_of_type(read_id, 'S', [&](int64_t id){
        // Check for cases that should be impossible
        if (result != -1){
            throw runtime_error("ERROR: multiple samples found for id: " + to_string(read_id));
        }

        result = id;
    });

    // For this project, every read must have exactly one sample.
    if (result == -1){
        throw runtime_error("ERROR: no sample found for id: " + to_string(read_id));
    }

    return result;
}


void TransMap::for_each_read(const function<void(const string& name, const BinarySequence<uint64_t>& sequence)>& f) const{
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
    graph.for_each_neighbor_of_type(id, type, [&](int64_t id2){
        f(id2);
    });
}


void TransMap::for_each_neighbor_of_type(int64_t id, char type, const function<void(int64_t id, float w)>& f) const{
    graph.for_each_neighbor_of_type(id, type, [&](int64_t id2, float w){
        f(id2, w);
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


void TransMap::for_each_path_of_read(int64_t read_id, const function<void(const string& path_name, int64_t path_id)>& f) const{
    graph.for_each_neighbor_of_type(read_id, 'P', [&](int64_t id){
        f(graph.get_node(id).name, id);
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


void TransMap::for_each_variant_of_path(int64_t path_id, const function<void(int64_t id)>& f) const{
    graph.for_each_neighbor_of_type(path_id, 'V', [&](int64_t id){
        f(id);
    });
}


string TransMap::get_sample_of_read(const string& read_name) const{
    string result;
    size_t i=0;
    graph.for_each_neighbor_of_type(read_name, 'S', [&](const HeteroNode& neighbor, int64_t id){
        if (i > 0){
            throw runtime_error("ERROR: multiple samples found for read: " + read_name);
        }
        result = neighbor.name;
        i++;
    });

    return result;
}


void TransMap::for_each_sample_of_path(const string& path_name, const function<void(const string& name, int64_t id)>& f) const{
    unordered_set<int64_t> visited;

    auto id = graph.name_to_id(path_name);
    graph.for_each_neighbor_of_type(id, 'R', [&](int64_t r){
        graph.for_each_neighbor_of_type(r, 'S', [&](int64_t s){
            if (not visited.count(s)){
                f(graph.get_node(s).name, s);
                visited.emplace(s);
            }
        });
    });
}


void TransMap::for_each_sample_of_path(int64_t id, const function<void(const string& name, int64_t id)>& f) const{
    unordered_set<int64_t> visited;

    graph.for_each_neighbor_of_type(id, 'R', [&](int64_t r){
        graph.for_each_neighbor_of_type(r, 'S', [&](int64_t s){
            if (not visited.count(s)){
                f(graph.get_node(s).name, s);
                visited.emplace(s);
            }
        });
    });
}


void TransMap::for_each_path_of_sample(const string& sample_name, const function<void(const string& name, int64_t id)>& f) const{
    unordered_set<int64_t> visited;

    auto id = graph.name_to_id(sample_name);
    graph.for_each_neighbor_of_type(id, 'R', [&](int64_t r){
        graph.for_each_neighbor_of_type(r, 'P', [&](int64_t p){
            if (not visited.count(p)){
                f(graph.get_node(p).name, p);
                visited.emplace(p);
            }
        });
    });
}


void TransMap::for_each_path_of_sample(int64_t sample_id, const function<void(const string& name, int64_t id)>& f) const{
    unordered_set<int64_t> visited;

    graph.for_each_neighbor_of_type(sample_id, 'R', [&](int64_t r){
        graph.for_each_neighbor_of_type(r, 'P', [&](int64_t p){
            if (not visited.count(p)){
                f(graph.get_node(p).name, p);
                visited.emplace(p);
            }
        });
    });
}


void TransMap::for_each_phased_variant_of_sample(const string& sample_name, const function<void(const string& name, int64_t id, bool phase)>& f) const {
    auto id = graph.name_to_id(sample_name);
    for_each_phased_variant_of_sample(id, f);
}


void TransMap::for_each_phased_variant_of_sample(int64_t sample_id, const function<void(const string& name, int64_t id, bool phase)>& f) const{
    set<int64_t> path_ids;

    // First find an arbitrary assignment of phases, which is sorted by path id
    for_each_path_of_sample(sample_id, [&](const string& path_name, int64_t p){
        path_ids.emplace(p);
    });

    if (path_ids.empty()){
//        cerr << "WARNING: no paths found for sample " << get_node(sample_id).name << '\n';
        return;
    }

    if (path_ids.size() > 2){
        string s;

        for (auto p: path_ids){
            s += graph.get_node(p).name + " ";
        }

        throw runtime_error("ERROR: more than two paths found for sample " + get_node(sample_id).name + ": " + s);
    }

    // Use begin and rbegin to return first and last element, which may be the same element if homozygous
    auto p = path_ids.begin();
    graph.for_each_neighbor_of_type(*p, 'V', [&](int64_t v) {
        f(graph.get_node(v).name, v, 0);
    });

    auto p2 = path_ids.rbegin();
    graph.for_each_neighbor_of_type(*p2, 'V', [&](int64_t v) {
        f(graph.get_node(v).name, v, 1);
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


void TransMap::for_edge_in_bfs(
        const string& start_name,
        float min_edge_weight,
        const function<bool(const HeteroNode& node)>& criteria,
        const function<void(const HeteroNode& a, int64_t a_id, const HeteroNode& b, int64_t b_id)>& f) const{
    graph.for_edge_in_bfs(start_name,min_edge_weight,criteria,f);
}


void TransMap::write_edge_info_to_csv(path output_path, const VariantGraph& variant_graph) const{
    ofstream out(output_path);

    // Write header
    out << "sample,read,read_length,path,path_length,weight\n";

    for_each_sample([&](const string& sample_name, int64_t sample_id){
        for_each_read_of_sample(sample_id, [&](int64_t read_id){
            for_each_path_of_read(read_id, [&](int64_t path_id){
                auto [success, weight] = try_get_edge_weight(read_id, path_id);

                if (not success){
                    throw runtime_error("ERROR: edge weight not found for read-path edge: " + to_string(read_id) + " " + to_string(path_id));
                }

                auto read_length = sequences.at(read_id).size();

                // Get the length of the path by summing the lengths of the nodes in the variant graph
                auto path_length = 0;

                // Convert path name to vector of ID/orientation pairs
                vector <pair <string,bool> > path;
                GafAlignment::parse_string_as_path(graph.get_node(path_id).name, path);

                // Sum the lengths of the nodes in the path
                for (const auto& [node_name, is_reverse]: path){
                    // Get hashgraph handle from name/orientation pair
                    auto h = variant_graph.graph.get_handle(stoll(node_name), is_reverse);
                    auto node_length = variant_graph.graph.get_length(h);
                    path_length += int32_t(node_length);
                }

                out << sample_name << ',' << get_node(read_id).name << ',' << read_length << ',' << get_node(path_id).name << ',' << path_length << ',' << weight << '\n';
            });
        });
    });
}


void TransMap::extract_sample_as_transmap(const string& sample_name, TransMap& result) {
    result.add_sample(sample_name);

    for_each_read_of_sample(sample_name, [&](const string& read_name, int64_t read_id) {
        result.add_read(read_name);

        // When generating the new transmap, must use string names, NOT IDs, because the IDs will differ betw transmaps
        result.add_edge(sample_name, read_name);

        for_each_path_of_read(read_id, [&](const string& path_name, int64_t path_id) {
            if (not result.has_node(path_name)) {
                result.add_path(path_name);
            }

            result.add_edge(read_name, path_name);
        });
    });
}


void TransMap::detangle_sample_paths(unordered_map<string,string>& hapmap){
    vector <pair <int64_t, int64_t> > edges_to_be_removed;
    vector <tuple <string, string, float> > edges_to_add;
    string result;

    for_each_path([&](const string& path_name, int64_t path_id) {
        unordered_map <int64_t,vector<int64_t> > sample_reads_of_path;

        // Find the reads that connect to this path and group them by sample. Each sample will get their own path.
        for_each_read_of_path(path_id, [&](int64_t read_id) {
            auto sample_id = get_read_sample(read_id);
            sample_reads_of_path[sample_id].push_back(read_id);
        });

        // Only care about paths that are used by multiple samples
        if (sample_reads_of_path.size() > 1) {
            for (const auto& [sample_id,reads]: sample_reads_of_path) {
                // Construct a child node to replace the parent node, because it needs to be untangled w.r.t. samples
                auto new_name = path_name + "_" + to_string(sample_id);
                hapmap.emplace(new_name, path_name);

                // Connect the reads to the new child path, and stage the old edges to be deleted (leaving the parent
                // path stranded as an island, will be reconnected later)
                for (auto read_id: reads) {
                    auto [_,w] = try_get_edge_weight(read_id, path_id);
                    edges_to_add.emplace_back(get_node(read_id).name, new_name, w);
                    edges_to_be_removed.emplace_back(read_id, path_id);
                }
            }
        }
    });

    for (auto [a,b]: edges_to_be_removed) {
        remove_edge(a,b);
    }

    for (const auto& [child_name,parent_name]: hapmap) {
        add_path(child_name);
    }

    for (const auto& [a,b,w]: edges_to_add) {
        add_edge(a,b,w);
    }
}


void TransMap::retangle_sample_paths(const unordered_map<string,string>& hapmap){
    // Reconnect any reads from children paths to their original parent path
    for (const auto& [child_path_name,parent_path_name]: hapmap) {
        auto [a_success,child_path_id] = try_get_id(child_path_name);

        // It's ok if some of the child paths have been deleted
        if (not a_success) {
            continue;
        }

        auto [b_success,parent_path_id] = try_get_id(parent_path_name);

        // It's NOT OK (!!) if the parent path corresponding to the child path has been deleted !!
        if (not b_success) {
            throw runtime_error("ERROR: cannot retangle TransMap with missing parent path: " + parent_path_name);
        }

        // Copy the incoming read-->child_path edges to the parent_path
        for_each_neighbor_of_type(child_path_id, 'R', [&](int64_t read_id, float w) {
            add_edge(read_id, parent_path_id, w);
        });

        // Delete the naughty children
        remove_node(child_path_id);
    }
}


bool operator==(const TransMap& a, const TransMap& b) {
    // All data must be compared on a name level to avoid inconsistencies in name<->id mapping, which are arbitrary
    unordered_set <pair <string,string> > a_edges;
    unordered_set <pair <string,string> > b_edges;
    unordered_set<string> a_nodes;
    unordered_set<string> b_nodes;

    a.for_edge_in_bfs(
        TransMap::sample_node_name,
        0,
        [&](const HeteroNode& node){return true;},
        [&](const HeteroNode& a_node, int64_t a_id, const HeteroNode& b_node, int64_t b_id) {
            a_edges.emplace(a_node.name, b_node.name);
            a_nodes.emplace(a_node.name);
            a_nodes.emplace(b_node.name);
    });

    b.for_edge_in_bfs(
        TransMap::sample_node_name,
        0,
        [&](const HeteroNode& node){return true;},
        [&](const HeteroNode& a_node, int64_t a_id, const HeteroNode& b_node, int64_t b_id) {
            b_edges.emplace(a_node.name, b_node.name);
            b_nodes.emplace(a_node.name);
            b_nodes.emplace(b_node.name);
    });

    return (a_nodes == b_nodes) and (a_edges == b_edges);
}


}
