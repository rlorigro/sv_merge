#include "TransitiveMap.hpp"
#include "gaf.hpp"
#include <fstream>
#include <algorithm>

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


bool TransMap::get_is_compressed() const{
    return is_compressed;
}


size_t TransMap::get_read_count() const{
    auto id = graph.name_to_id(read_node_name);
    return graph.get_edge_count(id);
}


size_t TransMap::get_path_count() const{
    auto id = graph.name_to_id(path_node_name);
    return graph.get_edge_count(id);
}


size_t TransMap::get_sample_count() const{
    auto id = graph.name_to_id(sample_node_name);
    return graph.get_edge_count(id);
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


size_t TransMap::get_edge_count() const{
    return graph.get_edge_count();
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


void TransMap::get_sample_of_read(const string& read_name, string& result) const{
    result.clear();

    // Using the HeteroGraph back end means that we iterate, even for a 1:1 mapping.
    graph.for_each_neighbor_of_type(read_name, 'S', [&](const HeteroNode& neighbor, int64_t id){
        // Check for cases that should be impossible
        if (not result.empty()){
            if (is_compressed){
                cerr << "WARNING: cannot expect 1:1 mapping from read->sample after using compression, even after decompression. Exiting..." << '\n';
            }

            throw runtime_error("ERROR: multiple samples found for read: " + read_name);
        }

        result = neighbor.name;
    });

    if (result.empty()){
        throw runtime_error("ERROR: no sample found for read: " + read_name);
    }
}


string TransMap::get_sample_of_read(const string& read_name) const{
    string result;

    // Using the HeteroGraph back end means that we iterate, even for a 1:1 mapping.
    graph.for_each_neighbor_of_type(read_name, 'S', [&](const HeteroNode& neighbor, int64_t id){
        // Check for cases that should be impossible
        if (not result.empty()){
            if (is_compressed){
                cerr << "WARNING: cannot expect 1:1 mapping from read->sample after using compression, even after decompression. Exiting..." << '\n';
            }

            throw runtime_error("ERROR: multiple samples found for read: " + read_name);
        }

        result = neighbor.name;
    });

    if (result.empty()){
        throw runtime_error("ERROR: no sample found for read: " + read_name);
    }

    return result;
}


void TransMap::get_sample_of_read(int64_t read_id, string& result) const{
    result.clear();

    if (graph.get_node(read_id).type != 'R'){
        throw runtime_error("ERROR: non-read ID provided for get_read_sample: " + to_string(read_id) + " " + graph.get_node(read_id).name);
    }

    // Using the HeteroGraph back end means that we iterate, even for a 1:1 mapping.
    graph.for_each_neighbor_of_type(read_id, 'S', [&](int64_t id){
        // Check for cases that should be impossible
        if (not result.empty()){
            if (is_compressed){
                cerr << "WARNING: cannot expect 1:1 mapping from read->sample after using compression, even after decompression. Exiting..." << '\n';
            }

            throw runtime_error("ERROR: multiple samples found for id: " + to_string(read_id));
        }

        result = graph.get_node(id).name;
    });

    // For this project, every read must have exactly one sample.
    if (result.empty()){
        throw runtime_error("ERROR: no sample found for id: " + to_string(read_id));
    }
}


int64_t TransMap::get_sample_of_read(int64_t read_id) const{
    if (graph.get_node(read_id).type != 'R'){
        throw runtime_error("ERROR: non-read ID provided for get_read_sample: " + to_string(read_id) + " " + graph.get_node(read_id).name);
    }

    int64_t result = -1;

    graph.for_each_neighbor_of_type(read_id, 'S', [&](int64_t id){
        // Check for cases that should be impossible
        if (result != -1){
            if (is_compressed){
                cerr << "WARNING: cannot expect 1:1 mapping from read->sample after using compression, even after decompression. Exiting..." << '\n';
            }

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


bool TransMap::has_sample(const string& sample_name) const {
    const auto result = graph.get_edges(graph.name_to_id(sample_node_name));
    for (const auto& [id_b, w]: result) {
        const auto& node = graph.get_node(id_b);
        if (node.type=='S' && node.name==sample_name) return true;
    }
    return false;
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


void TransMap::for_each_sample_of_read(const string& read_name, const function<void(const string& name, int64_t sample_id)>& f) const{
    graph.for_each_neighbor_of_type(read_name, 'S', [&](const HeteroNode& neighbor, int64_t sample_id){
        f(neighbor.name, sample_id);
    });
}


void TransMap::for_each_sample_of_read(const int64_t& read_id, const function<void(int64_t sample_id)>& f) const{
    graph.for_each_neighbor_of_type(read_id, 'S', [&](int64_t sample_id){ f(sample_id); });
}


void TransMap::for_each_variant_of_path(int64_t path_id, const function<void(int64_t id)>& f) const{
    graph.for_each_neighbor_of_type(path_id, 'V', [&](int64_t id){
        f(id);
    });
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


void TransMap::write_edge_info_to_csv(path output_path, const VariantGraph& variant_graph, bool use_sample_id) const{
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

                string s = sample_name;

                if (use_sample_id) {
                    s = to_string(sample_id);
                }

                out << s << ',' << get_node(read_id).name << ',' << read_length << ',' << get_node(path_id).name << ',' << path_length << ',' << weight << '\n';
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
            auto sample_id = get_sample_of_read(read_id);
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



/// WARNING: DOES NOT COMPARE EDGE WEIGHTS, ONLY COMPARES EDGE PRESENCE/ABSENCE
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


pair<int64_t,int64_t> TransMap::get_n_paths_of_read(int64_t read_id) const {
    int64_t n_paths = 0;
    int64_t last_path = -1;
    for_each_path_of_read(read_id, [&](int64_t path_id) { n_paths++; last_path=path_id; });
    return std::make_pair(n_paths,last_path);
}


bool TransMap::are_edges_distinct() const {
    return graph.are_edges_distinct();
}


void TransMap::sort_adjacency_lists() {
    graph.sort_adjacency_lists();
}


void TransMap::update_first_of_type() {
    graph.update_first_of_type();
}


void TransMap::partition(vector<TransMap>& maps) {
    size_t i;
    size_t set_size;
    int64_t id, type, component_id;
    string sample_name, s_name;
    vector<int64_t> stack, component_size;
    unordered_set<int64_t> set;
    unordered_map<int64_t, int64_t> read_component, path_component;

    graph.update_first_of_type();

    // Computing connected components and building the corresponding transmaps
    maps.clear();
    component_id=-1;
    for_each_read([&](const string& read_name, int64_t read_id) {
        if (read_component.contains((read_id))) return;
        read_component.emplace(read_id,++component_id);
        maps.emplace_back();
        TransMap& new_map = maps.back();
        new_map.add_read(read_name);
        sample_name=get_sample_of_read(read_name);
        new_map.add_sample(sample_name);
        new_map.add_edge(read_name,sample_name);
        stack.clear(); stack.emplace_back(read_id); stack.emplace_back(0);
        while (!stack.empty()) {
            type=stack.at(stack.size()-1); stack.pop_back();
            id=stack.at(stack.size()-1); stack.pop_back();
            if (type==0) {
                const string& r_name = graph.get_node(id).name;
                for_each_path_of_read(id, [&](int64_t p_id) {
                    if (path_component.contains(p_id)) return;
                    path_component.emplace(p_id,component_id);
                    stack.emplace_back(p_id); stack.emplace_back(1);
                    const string& p_name = graph.get_node(p_id).name;
                    new_map.add_path(p_name);
                    new_map.add_edge(r_name,p_name,try_get_edge_weight(id,p_id).second);
                });
            }
            else {
                const string& p_name = graph.get_node(id).name;
                for_each_read_of_path(id, [&](int64_t r_id) {
                    if (read_component.contains(r_id)) return;
                    read_component.emplace(r_id,component_id);
                    stack.emplace_back(r_id); stack.emplace_back(0);
                    const string& r_name = graph.get_node(r_id).name;
                    new_map.add_read(r_name);
                    new_map.add_edge(r_name,p_name,try_get_edge_weight(id,r_id).second);
                    s_name=get_sample_of_read(r_name);
                    if (!new_map.has_sample(s_name)) new_map.add_sample(s_name);
                    new_map.add_edge(s_name,r_name);
                });
            }
        }
        component_size.emplace_back(new_map.get_read_count());
        component_size.emplace_back(new_map.get_path_count());
        component_size.emplace_back(new_map.get_sample_count());
        component_size.emplace_back(new_map.get_edge_count());
    });
    // cerr << "Number of connected components: " << to_string(component_id+1) << '\n';
    // cerr << "Component \t n_reads \t n_paths \t n_samples \t n_edges\n";
    // for (i=0; i<component_size.size(); i+=4) cerr << to_string(i/4) << '\t' << to_string(component_size.at(i)) << '\t' << to_string(component_size.at(i+1)) << '\t' << to_string(component_size.at(i+2)) << '\t' << to_string(component_size.at(i+3)) << '\n';

    // Collecting samples assigned to 2 components
    partitioned_samples.clear();
    for_each_sample([&](const string& sample_name, int64_t sample_id) {
        set.clear();
        for_each_read_of_sample(sample_name, [&](const string& read_name, int64_t read_id) {
            set.emplace(read_component.at(read_id));
        });
        set_size=set.size();
        if (set_size>2) throw runtime_error("ERROR: the reads of sample "+sample_name+" were partitioned into "+to_string(set.size())+" connected components");
        else if (set_size==2) partitioned_samples.emplace_back(sample_name);
    });
    // cerr << "Number of samples assigned to 2 components: " << to_string(partitioned_samples.size()) << '\n';
}


TransMap TransMap::partition_get_test_transmap() {
    TransMap out;
    out.add_sample("sample1");
    out.add_read("s1r1"); out.add_read("s1r2"); out.add_read("s1r3");
    out.add_edge("sample1","s1r1"); out.add_edge("sample1","s1r2"); out.add_edge("sample1","s1r3");
    out.add_sample("sample2");
    out.add_read("s2r1"); out.add_read("s2r2"); out.add_read("s2r3");
    out.add_edge("sample2","s2r1"); out.add_edge("sample2","s2r2"); out.add_edge("sample2","s2r3");
    out.add_sample("sample3");
    out.add_read("s3r1"); out.add_read("s3r2"); out.add_read("s3r3");
    out.add_edge("sample3","s3r1"); out.add_edge("sample3","s3r2"); out.add_edge("sample3","s3r3");

    out.add_path("path1"); out.add_path("path2"); out.add_path("path3"); out.add_path("path4");
    out.add_edge("s1r1","path1");
    out.add_edge("s1r2","path2");
    out.add_edge("s1r3","path3");

    out.add_edge("s2r1","path2");
    out.add_edge("s2r2","path3");
    out.add_edge("s2r3","path3");

    out.add_edge("s3r1","path1");
    out.add_edge("s3r2","path4");
    out.add_edge("s3r3","path4");
    return out;
}


void TransMap::compress_reads(float weight_quantum, bool verbose) {
    is_compressed = true;
    size_t length;
    int64_t i, j, k;
    int64_t read_id, n_clusters, n_reads;
    vector<bool> is_redundant, is_representative;
    vector<int64_t> node_ids, cluster_representative, cluster_size, sample_ids;
    vector<vector<int64_t>> neighbors;
    vector<vector<float>> weights, compared_weights;

    graph.update_first_of_type();

    // Collecting all read-haplotype edges, with reads grouped by sample.
    n_reads=get_read_count();
    neighbors.reserve(n_reads);
    for (i=0; i<n_reads; i++) neighbors.emplace_back();
    weights.reserve(n_reads);
    for (i=0; i<n_reads; i++) weights.emplace_back();
    compared_weights.reserve(n_reads);
    for (i=0; i<n_reads; i++) compared_weights.emplace_back();
    node_ids.reserve(n_reads); sample_ids.reserve(n_reads);
    i=-1;
    for_each_sample([&](const string &sample_name, int64_t sample_id) {
        if (solved_samples.contains(sample_id)) return;
        for_each_read_of_sample(sample_id, [&](int64_t read_id) {
            node_ids.emplace_back(read_id);
            sample_ids.emplace_back(sample_id);
            i++;
            graph.for_each_neighbor_of_type(read_id, 'P', [&](int64_t path_id) {
                // For each read, its neighbors lie in the same global order.
                auto [success, weight] = try_get_edge_weight(read_id, path_id);
                neighbors.at(i).emplace_back(path_id);
                weights.at(i).emplace_back(weight);
                compared_weights.at(i).emplace_back(get_edge_weight(weight,weight_quantum));
            });
        });
    });
    n_reads=i+1;
    if (verbose) {
        cerr << "Read-hap weights before compression: \n";
        for (i=0; i<n_reads; i++) {
            read_id=node_ids.at(i);
            cerr << "read=" << to_string(read_id) << " weights=";
            for (j=0; j<weights.at(i).size(); j++) cerr << "(" << to_string(neighbors.at(i).at(j)) << "," << to_string(weights.at(i).at(j)) << "), ";
            cerr << '\n';
        }
    }

    // Clustering reads; removing redundant reads; computing new read-hap weights.
    is_redundant.reserve(n_reads);
    for (i=0; i<n_reads; i++) is_redundant.emplace_back(false);
    is_representative.reserve(n_reads);
    for (i=0; i<n_reads; i++) is_representative.emplace_back(false);
    cluster_size.reserve(n_reads);
    for (i=0; i<n_reads; i++) cluster_size.emplace_back(0);
    n_clusters=0;
    for (i=0; i<n_reads; i++) {
        if (is_redundant.at(i)) continue;
        n_clusters++;
        cluster_size.at(i)=1;
        read_id=node_ids.at(i);
        length=neighbors.at(i).size();
        for (j=i+1; j<n_reads; j++) {
            if (sample_ids.at(j)!=sample_ids.at(i)) break;
            if (is_redundant.at(j) || neighbors.at(j)!=neighbors.at(i) || compared_weights.at(j)!=compared_weights.at(i)) continue;
            is_redundant.at(j)=true; is_representative.at(i)=true;
            cluster_size.at(i)++;
            for (k=0; k<length; k++) weights.at(i).at(k)+=weights.at(j).at(k);
        }
    }
    compared_weights.clear();
    // cerr << "n_reads=" << to_string(n_reads) << " -> n_read_clusters=" << to_string(n_clusters) <<  " (after compressing reads)\n";
    if (verbose) {
        cerr << "Read-hap weights after compression: \n";
        for (i=0; i<n_reads; i++) {
            read_id=node_ids.at(i);
            cerr << "read=" << to_string(read_id) << " weights=";
            for (j=0; j<weights.at(i).size(); j++) cerr << "(" << to_string(neighbors.at(i).at(j)) << "," << to_string(weights.at(i).at(j)) << "), ";
            cerr << '\n';
        }
    }

    // Updating the transmap.
    // Remark: edge removal could be implemented much faster.
    for (i=0; i<n_reads; i++) {
        if (is_redundant.at(i) || !is_representative.at(i)) continue;
        read_id=node_ids.at(i); length=neighbors.at(i).size();
        for (j=0; j<length; j++) graph.update_edge_weight(read_id, neighbors.at(i).at(j), weights.at(i).at(j));
    }
    for (i=0; i<n_reads; i++) {
        if (is_redundant.at(i)) remove_node(node_ids.at(i));
    }
}


/**
 * Implemented as a quadratic scan over all the reads in the cohort. Might be too slow for large cohorts.
 */
void TransMap::compress_samples(float weight_quantum) {
    is_compressed = true;
    bool contained, trivial;
    size_t length;
    int64_t i, j, k;
    int64_t read_id, sample_id, s_id, next_i, next_j, n_clusters, n_reads;
    int64_t contained_first, contained_last, container_id, container_first, container_last;
    vector<bool> used;
    vector<int64_t> sample_ids, to_remove, cluster_size, cluster_size_prime;
    vector<vector<int64_t>> neighbors;
    vector<vector<float>> weights, compared_weights;
    unordered_map<int64_t,tuple<int64_t,int64_t,int64_t,int64_t,int64_t>> sample_to_identical_sample, sample_to_container_sample;

    graph.update_first_of_type();

    // Collecting all read-haplotype edges, with reads grouped by unsolved sample.
    n_reads=get_read_count();
    neighbors.reserve(n_reads);
    for (i=0; i<n_reads; i++) neighbors.emplace_back();
    weights.reserve(n_reads);
    for (i=0; i<n_reads; i++) weights.emplace_back();
    compared_weights.reserve(n_reads);
    for (i=0; i<n_reads; i++) compared_weights.emplace_back();
    read_ids.clear(); read_ids.reserve(n_reads); sample_ids.reserve(n_reads);
    i=-1;
    for_each_sample([&](const string& sample_name, int64_t sample_id) {
        if (solved_samples.contains(sample_id)) return;
        for_each_read_of_sample(sample_id, [&](int64_t read_id) {
            read_ids.emplace_back(read_id); sample_ids.emplace_back(get_id(sample_name));
            i++;
            graph.for_each_neighbor_of_type(read_id, 'P', [&](int64_t path_id) {
                // For each read, its neighbors lie in the same global order.
                auto [success, weight] = try_get_edge_weight(read_id,path_id);
                neighbors.at(i).emplace_back(path_id);
                weights.at(i).emplace_back(weight);
                compared_weights.at(i).emplace_back(get_edge_weight(weight,weight_quantum));
            });
        });
    });
    n_reads=i+1;

    // Finding equivalent reads across all unsolved samples
    cluster_ids.clear(); cluster_ids.reserve(n_reads);
    for (i=0; i<n_reads; i++) cluster_ids.emplace_back(0);
    n_clusters=0;
    for (i=0; i<n_reads; i++) {
        if (cluster_ids.at(i)!=0) continue;
        n_clusters++;
        cluster_ids.at(i)=n_clusters;
        read_id=read_ids.at(i);
        for (j=i+1; j<n_reads; j++) {
            if (cluster_ids.at(j)!=0 || neighbors.at(j)!=neighbors.at(i) || compared_weights.at(j)!=compared_weights.at(i)) continue;
            cluster_ids.at(j)=n_clusters;
        }
    }
    compared_weights.clear();
    // cerr << "n_reads=" << to_string(n_reads) << " -> n_read_clusters=" << to_string(n_clusters) << " (across all unsolved samples)\n";

    cluster_size.reserve(n_clusters);
    for (i=0; i<n_clusters; i++) cluster_size.emplace_back(0);
    cluster_size_prime.reserve(n_clusters);
    for (i=0; i<n_clusters; i++) cluster_size_prime.emplace_back(0);

    // Removing identical samples
    i=0;
    while (i<n_reads) {
        sample_id=sample_ids.at(i);
        next_i=i+1;
        while (next_i<n_reads) {
            if (sample_ids.at(next_i)!=sample_id) break;
            next_i++;
        }
        if (sample_to_identical_sample.contains(sample_id)) { i=next_i; continue; }
        for (j=0; j<n_clusters; j++) cluster_size.at(j)=0;
        for (j=i; j<next_i; j++) cluster_size.at(cluster_ids.at(j)-1)++;
        j=next_i;
        while (j<n_reads) {
            s_id=sample_ids.at(j);
            next_j=j+1;
            while (next_j<n_reads) {
                if (sample_ids.at(next_j)!=s_id) break;
                next_j++;
            }
            if (sample_to_identical_sample.contains(s_id)) { j=next_j; continue; }
            for (k=0; k<n_clusters; k++) cluster_size_prime.at(k)=0;
            for (k=j; k<next_j; k++) cluster_size_prime.at(cluster_ids.at(k)-1)++;
            if (cluster_size_prime==cluster_size) {
                to_remove.emplace_back(s_id);
                for (k=j; k<next_j; k++) to_remove.emplace_back(read_ids.at(k));
                sample_to_identical_sample.emplace(s_id,std::make_tuple(sample_id,j,next_j-1,i,next_i-1));
            }
            j=next_j;
        }
        i=next_i;
    }
    if (!sample_to_identical_sample.empty()) cerr << "Removed " << to_string(sample_to_identical_sample.size()) << " identical samples\n";

    // Removing contained samples where every read can be assigned to only one hap
    i=0;
    while (i<n_reads) {
        sample_id=sample_ids.at(i);
        next_i=i+1;
        while (next_i<n_reads) {
            if (sample_ids.at(next_i)!=sample_id) break;
            next_i++;
        }
        if (sample_to_identical_sample.contains(sample_id)) { i=next_i; continue; }
        trivial=true;
        for (j=i; j<next_i; j++) {
            if (neighbors.at(j).size()>1) { trivial=false; break; }
        }
        if (!trivial) { i=next_i; continue; }
        for (j=0; j<n_clusters; j++) cluster_size.at(j)=0;
        for (j=i; j<next_i; j++) cluster_size.at(cluster_ids.at(j)-1)++;
        j=next_i;
        while (j<n_reads) {
            s_id=sample_ids.at(j);
            next_j=j+1;
            while (next_j<n_reads) {
                if (sample_ids.at(next_j)!=s_id) break;
                next_j++;
            }
            if (sample_to_identical_sample.contains(s_id)) { j=next_j; continue; }
            for (k=0; k<n_clusters; k++) cluster_size_prime.at(k)=0;
            for (k=j; k<next_j; k++) cluster_size_prime.at(cluster_ids.at(k)-1)++;
            contained=true;
            for (k=0; k<n_clusters; k++) {
                if (cluster_size.at(k)>cluster_size_prime.at(k)) { contained=false; break; }
            }
            if (contained) {
                to_remove.emplace_back(sample_id);
                for (k=i; k<next_i; k++) to_remove.emplace_back(read_ids.at(k));
                sample_to_container_sample.emplace(sample_id,std::make_tuple(s_id,i,next_i-1,j,next_j-1));
                break;
            }
            j=next_j;
        }
        i=next_i;
    }
    // if (!sample_to_container_sample.empty()) cerr << "Removed " << to_string(sample_to_container_sample.size()) << " contained samples with trivial haplotypes\n";

    // Preparing data structures for decompression
    sample_to_sample.clear();
    for (auto& pair: sample_to_identical_sample) {
        auto& t = sample_to_identical_sample.at(pair.first);
        container_id=std::get<0>(t);
        contained_first=std::get<1>(t); contained_last=std::get<2>(t);
        container_first=std::get<3>(t); container_last=std::get<4>(t);
        while (sample_to_container_sample.contains(container_id)) {
            t=sample_to_container_sample.at(container_id);
            container_id=std::get<0>(sample_to_container_sample.at(container_id));
            container_first=std::get<3>(t); container_last=std::get<4>(t);
        }
        sample_to_sample.emplace(get_node(pair.first).name,std::make_tuple(contained_first,contained_last,container_first,container_last));
    }
    for (auto& pair: sample_to_container_sample) {
        auto& t = sample_to_container_sample.at(pair.first);
        container_id=std::get<0>(t);
        contained_first=std::get<1>(t); contained_last=std::get<2>(t);
        container_first=std::get<3>(t); container_last=std::get<4>(t);
        while (sample_to_container_sample.contains(container_id)) {
            t=sample_to_container_sample.at(container_id);
            container_id=std::get<0>(sample_to_container_sample.at(container_id));
            container_first=std::get<3>(t); container_last=std::get<4>(t);
        }
        sample_to_sample.emplace(get_node(pair.first).name,std::make_tuple(contained_first,contained_last,container_first,container_last));
    }

    // Updating the transmap.
    // Remark: edge removal could be implemented much faster.
    for (auto& pair: sample_to_sample) compress_samples_update_weights(std::get<0>(pair.second),std::get<1>(pair.second),std::get<2>(pair.second),std::get<3>(pair.second),weights,used);
    for (i=0; i<n_reads; i++) {
        sample_id=sample_ids.at(i);
        if (sample_to_identical_sample.contains(sample_id) || sample_to_container_sample.contains(sample_id)) continue;
        read_id=read_ids.at(i);
        length=neighbors.at(i).size();
        for (j=0; j<length; j++) graph.update_edge_weight(read_id,neighbors.at(i).at(j),weights.at(i).at(j));
    }
    for (auto node_id: to_remove) remove_node(node_id);
}


void TransMap::compress_samples_update_weights(int64_t from_first, int64_t from_last, int64_t to_first, int64_t to_last, vector<vector<float>>& weights, vector<bool>& used) {
    is_compressed = true;
    int64_t i, j, k;
    int64_t cluster_id;

    used.clear();
    for (i=to_first; i<=to_last; i++) used.emplace_back(false);
    for (i=from_first; i<=from_last; i++) {
        cluster_id=cluster_ids.at(i);
        for (j=to_first; j<=to_last; j++) {
            if (used.at(j-to_first) || cluster_ids.at(j)!=cluster_id) continue;
            used.at(j-to_first)=true;
            for (k=0; k<weights.at(i).size(); k++) weights.at(j).at(k)+=weights.at(i).at(k);
            break;
        }
    }
}


/**
 * Essentially the same as `compress_samples_update_weights()`.
 */
void TransMap::decompress_samples() {
    int64_t i, j;
    int64_t cluster_id, sample_id, from_first, from_last, to_first, to_last;
    string sample_name;
    vector<bool> used;

    for (auto& pair: sample_to_sample) {
        sample_name=pair.first;
        add_sample(sample_name);
        sample_id=get_id(sample_name);
        auto& t = pair.second;
        from_first=std::get<0>(t); from_last=std::get<1>(t); to_first=std::get<2>(t); to_last=std::get<3>(t);
        used.clear();
        for (i=to_first; i<=to_last; i++) used.emplace_back(false);
        for (i=from_first; i<=from_last; i++) {
            cluster_id=cluster_ids.at(i);
            for (j=to_first; j<=to_last; j++) {
                if (used.at(j-to_first) || cluster_ids.at(j)!=cluster_id) continue;
                used.at(j-to_first)=true;
                add_edge(sample_id,read_ids.at(j),0);  // Updates both sides of the edge
                break;
            }
        }
    }
    read_ids.clear(); cluster_ids.clear(); sample_to_sample.clear();
}


int64_t TransMap::get_mandatory_haplotypes() {
    graph.update_first_of_type();

    int64_t out = 0;
    for_each_read([&](const string& read_name, int64_t read_id) {
        auto p = get_n_paths_of_read(read_id);
        if (p.first==1) {
            present_haps.emplace(p.second);
            present_edges.emplace(read_id,p.second);
            out++;
        }
    });
    return out;
}


void TransMap::compress_haplotypes_global(float weight_quantum) {
    is_compressed = true;
    const int64_t n_samples = get_sample_count();
    int64_t n_haps = get_path_count();

    int64_t i, j;
    int64_t n_equivalent, n_contained, hap_id;
    vector<bool> removed;
    vector<int64_t> hap_ids;
    unordered_set<int64_t> to_remove;
    vector<vector<int64_t>> neighbors;
    vector<vector<float>> weights;
    unordered_map<int64_t,unordered_set<int64_t>> edges_to_remove;

    graph.update_first_of_type();

    // Collecting all haplotype-read edges
    neighbors.reserve(n_haps);
    for (i=0; i<n_haps; i++) neighbors.emplace_back();
    weights.reserve(n_haps);
    for (i=0; i<n_haps; i++) weights.emplace_back();
    hap_ids.clear(); hap_ids.reserve(n_haps);
    i=-1;
    for_each_path([&](const string& path_name, int64_t path_id) {
        hap_ids.emplace_back(path_id);
        i++;
        for_each_read_of_path(path_id, [&](int64_t read_id) {
            // For each hap, its neighbors lie in the same global order.
            auto [success, weight] = try_get_edge_weight(read_id,path_id);
            neighbors.at(i).emplace_back(read_id);
            weights.at(i).emplace_back(get_edge_weight(weight,weight_quantum));
        });
    });

    // Removing globally-equivalent haps
    n_equivalent=0;
    removed.clear(); removed.reserve(n_haps);
    for (i=0; i<n_haps; i++) removed.emplace_back(false);
    for (i=0; i<n_haps; i++) {
        if (removed.at(i)) continue;
        for (j=i+1; j<n_haps; j++) {
            if (removed.at(j) || neighbors.at(j)!=neighbors.at(i) || weights.at(j)!=weights.at(i)) continue;
            n_equivalent++;
            removed.at(j)=true;
            to_remove.emplace(hap_ids.at(j));
        }
    }
    // if (n_equivalent!=0) cerr << "Removed " << to_string(n_equivalent) << " globally-equivalent haplotypes (out of " << to_string(n_haps) << " total)\n";

    // Removing globally-contained haps
    n_contained=0;
    for (i=0; i<n_haps; i++) {
        if (removed.at(i)) continue;
        for (j=i+1; j<n_haps; j++) {
            if (removed.at(j)) continue;
            if (is_haplotype_contained(i,j,neighbors,weights)) {
                n_contained++;
                removed.at(i)=true;
                to_remove.emplace(hap_ids.at(i));
                break;
            }
        }
    }
    neighbors.clear(); weights.clear(); hap_ids.clear();
    // if (n_contained!=0) cerr << "Removed " << to_string(n_contained) << " globally-contained haplotypes (out of " << to_string(n_haps) << " total)\n";

    // Updating the transmap.
    // Remark: edge removal could be implemented much faster.
    for (auto node_id: to_remove) remove_node(node_id);
}


bool TransMap::is_haplotype_contained(int64_t from, int64_t to, const vector<vector<int64_t>>& neighbors, const vector<vector<float>>& weights) {
    const size_t length_from = neighbors.at(from).size();
    const size_t length_to = neighbors.at(to).size();
    if (length_from>length_to) return false;

    bool found;
    size_t i, j;
    int64_t neighbor_from, neighbor_to;

    j=0;
    for (i=0; i<length_from; i++) {
        neighbor_from=neighbors.at(from).at(i);
        found=false;
        while (j<length_to) {
            neighbor_to=neighbors.at(to).at(j);
            if (neighbor_to>neighbor_from) break;
            else if (neighbor_to<neighbor_from) j++;
            else {
                found=true;
                if (weights.at(to).at(j)>weights.at(from).at(i)) return false;
                break;
            }
        }
        if (!found) return false;
    }
    return true;
}


bool TransMap::has_large_weight(float n_weight, float d_weight) {
    const float DELTA = n_weight/d_weight;
    return graph.has_large_weight(get_id(read_node_name),'P',DELTA);
}


void TransMap::compress_haplotypes_local(float n_weight, float d_weight, float weight_quantum) {
    is_compressed = true;
    const float DELTA = n_weight/d_weight;
    const int64_t N_SAMPLES = get_sample_count();
    // cerr << "compress_haplotypes_local> DELTA=" << to_string(DELTA) << '\n';

    bool large_weight;
    int64_t i, j;
    int64_t n_dominated, hap_id, n_removed_edges, n_reads, read_id;
    float weight;
    size_t length, n_haps;
    vector<bool> removed;
    vector<int64_t> hap_ids, read_ids, tmp_vector;
    vector<pair<int64_t,int64_t>> edges_to_remove_prime;
    unordered_set<int64_t> to_remove;
    vector<vector<int64_t>> neighbors, neighbors_prime;
    vector<vector<float>> weights, weights_prime;
    unordered_map<int64_t,unordered_set<int64_t>> edges_to_remove;

    graph.update_first_of_type();

    // Removing locally-dominated haps
    n_dominated=0;
    for_each_sample([&](const string& sample_name, int64_t sample_id) {
        // Building the read-hap matrix
        read_ids.clear(); neighbors_prime.clear(); weights_prime.clear();
        i=-1; large_weight=false;
        for_each_read_of_sample(sample_id, [&](int64_t read_id) {
            // Reads are assumed to be enumerated in sorted order by ID.
            i++; read_ids.emplace_back(read_id);
            neighbors_prime.emplace_back(); weights_prime.emplace_back();
            for_each_path_of_read(read_id, [&](int64_t path_id) {
                // Haps are assumed to be enumerated in sorted order by ID.
                auto [success, w] = try_get_edge_weight(read_id,path_id);
                neighbors_prime.at(i).emplace_back(path_id);
                w=get_edge_weight(w,weight_quantum);
                weights_prime.at(i).emplace_back(w);
                if (w>=DELTA) large_weight=true;
            });
        });
        n_reads=i+1;
        if (n_reads==0 || !large_weight) return;
        // Building the hap-read matrix
        tmp_vector.clear();
        for (i=0; i<n_reads; i++) tmp_vector.insert(tmp_vector.end(),neighbors_prime.at(i).begin(),neighbors_prime.at(i).end());
        length=tmp_vector.size();
        if (length==0) return;
        if (length>1) std::sort(tmp_vector.begin(),tmp_vector.end());
        hap_ids.clear(); hap_ids.emplace_back(tmp_vector.at(0));
        for (i=1; i<tmp_vector.size(); i++) {
            if (tmp_vector.at(i)!=tmp_vector.at(i-1)) hap_ids.emplace_back(tmp_vector.at(i));
        }
        n_haps=hap_ids.size();
        neighbors.clear(); weights.clear();
        for (i=0; i<n_haps; i++) neighbors.emplace_back();
        for (i=0; i<n_haps; i++) weights.emplace_back();
        for (i=0; i<n_reads; i++) {
            read_id=read_ids.at(i);
            length=neighbors_prime.at(i).size();
            for (j=0; j<length; j++) {
                hap_id=neighbors_prime.at(i).at(j);
                weight=weights_prime.at(i).at(j);
                auto p = std::lower_bound(hap_ids.begin(),hap_ids.end(),hap_id);  // Could have been implemented as a mergesort since both arrays are sorted
                neighbors.at(p-hap_ids.begin()).emplace_back(read_id);
                weights.at(p-hap_ids.begin()).emplace_back(weight);
            }
        }
        neighbors_prime.clear(); weights_prime.clear(); read_ids.clear();
        // Removing dominated haps
        removed.clear();
        for (i=0; i<n_haps; i++) removed.emplace_back(false);
        for (i=0; i<n_haps; i++) {
            if (removed.at(i)) continue;
            for (j=i+1; j<n_haps; j++) {
                if (removed.at(j)) continue;
                if (is_haplotype_dominated(DELTA,i,j,neighbors,weights)) {
                    n_dominated++;
                    removed.at(i)=true;
                    hap_id=hap_ids.at(i);
                    if (edges_to_remove.contains(hap_id)) edges_to_remove.at(hap_id).emplace(sample_id);
                    else {
                        unordered_set<int64_t> v;
                        v.emplace(sample_id);
                        edges_to_remove.emplace(hap_id,v);
                    }
                    break;
                }
            }
        }
        neighbors.clear(); weights.clear(); hap_ids.clear();
    });
    // if (n_dominated!=0) cerr << "Found " << to_string(n_dominated) << " locally-dominated haplotypes (" << to_string(((float)n_dominated)/N_SAMPLES) << " avg per sample)\n";

    // Updating the transmap.
    // Remark: edge removal could be implemented much faster.
    edges_to_remove_prime.clear();
    n_removed_edges=0;
    for (auto& pair: edges_to_remove) {
        hap_id=pair.first;
        for_each_read_of_path(hap_id, [&](int64_t read_id) {
            if (pair.second.contains(get_sample_of_read(read_id))) {
                edges_to_remove_prime.emplace_back(read_id,hap_id);
                n_removed_edges++;
            }
        });
    }
    for (auto& pair: edges_to_remove_prime) remove_edge(pair.first,pair.second);
    // if (n_removed_edges!=0) cerr << "Removed " << to_string(n_removed_edges) << " edges to locally-dominated haplotypes\n";
}


bool TransMap::is_haplotype_dominated(float delta, int64_t from, int64_t to, const vector<vector<int64_t>>& neighbors, const vector<vector<float>>& weights) {
    const size_t length_from = neighbors.at(from).size();
    const size_t length_to = neighbors.at(to).size();
    if (length_from>length_to) return false;

    bool found;
    size_t i, j;
    int64_t neighbor_from, neighbor_to;

    j=0;
    for (i=0; i<length_from; i++) {
        neighbor_from=neighbors.at(from).at(i);
        found=false;
        while (j<length_to) {
            neighbor_to=neighbors.at(to).at(j);
            if (neighbor_to>neighbor_from) break;
            else if (neighbor_to<neighbor_from) j++;
            else {
                found=true;
                if (weights.at(to).at(j)>weights.at(from).at(i)-delta) return false;
                break;
            }
        }
        if (!found) return false;
    }
    return true;
}


void TransMap::solve_easy_samples(float n_weight, float d_weight, float weight_quantum) {
    const float DELTA = n_weight/d_weight;
    const int64_t N_HAPS = get_path_count();

    bool large_weight, not_worse, favored;
    int64_t i, j, k;
    int64_t read_id, hap_id, min_hap_id, n_reads, hap1, hap2, n_one_hap_samples, n_two_hap_samples;
    size_t n_haps;
    float weight, min_weight;
    vector<int64_t> read_ids;
    unordered_set<int64_t> tmp_set;
    vector<vector<int64_t>> neighbors;
    vector<vector<float>> weights;
    vector<pair<int64_t,int64_t>> edges_to_remove;
    unordered_set<int64_t> favored_haps;
    unordered_map<int64_t,unordered_set<int64_t>> hap_to_reads;

    graph.update_first_of_type();

    // Solving easy samples
    solved_samples.clear();
    n_one_hap_samples=0; n_two_hap_samples=0;
    edges_to_remove.clear();
    hap_to_reads.reserve(N_HAPS); favored_haps.reserve(N_HAPS); tmp_set.reserve(N_HAPS);
    for_each_sample([&](const string& sample_name, int64_t sample_id) {
        read_ids.clear(); neighbors.clear(); weights.clear();
        i=-1; large_weight=false;
        for_each_read_of_sample(sample_id, [&](int64_t read_id) {
            i++; read_ids.emplace_back(read_id);
            neighbors.emplace_back(); weights.emplace_back();
            for_each_path_of_read(read_id, [&](int64_t path_id) {
                // For each read, its neighbors lie in the same global order.
                auto [success, w] = try_get_edge_weight(read_id,path_id);
                neighbors.at(i).emplace_back(path_id);
                w=get_edge_weight(w,weight_quantum);
                weights.at(i).emplace_back(w);
                if (w>=DELTA) large_weight=true;
            });
        });
        n_reads=read_ids.size();
        if (n_reads==0 || !large_weight) return;
        hap_to_reads.clear(); favored_haps.clear();
        for (i=0; i<n_reads; i++) {
            read_id=read_ids.at(i);
            n_haps=neighbors.at(i).size();
            for (j=0; j<n_haps; j++) {
                hap_id=neighbors.at(i).at(j);
                weight=weights.at(i).at(j);
                not_worse=true;
                for (k=0; k<n_haps; k++) {
                    if (k!=j && weight>weights.at(i).at(k)) { not_worse=false; break; }
                }
                if (not_worse) {
                    if (hap_to_reads.contains(hap_id)) hap_to_reads.at(hap_id).insert(read_id);
                    else {
                        unordered_set<int64_t> new_set;
                        new_set.insert(read_id);
                        hap_to_reads[hap_id]=new_set;
                    }
                    favored=true;
                    for (k=0; k<n_haps; k++) {
                        if (k!=j && weight>weights.at(i).at(k)-DELTA) { favored=false; break; }
                    }
                    if (favored) favored_haps.insert(hap_id);
                }
            }
        }
        // One-hap sample
        hap1=-1;
        for (auto& hap_id: favored_haps) {
            if (hap_to_reads.at(hap_id).size()==n_reads) { hap1=hap_id; break; }
        }
        if (hap1!=-1) {
            for (i=0; i<n_reads; i++) {
                read_id=read_ids.at(i);
                n_haps=neighbors.at(i).size();
                for (j=0; j<n_haps; j++) {
                    hap_id=neighbors.at(i).at(j);
                    if (hap_id!=hap1) edges_to_remove.emplace_back(read_id,hap_id);
                    else present_edges.emplace(read_id,hap_id);
                }
            }
            present_haps.emplace(hap1);
            solved_samples.insert(sample_id);
            n_one_hap_samples++;
            return;
        }
        // Two-hap sample
        hap1=-1; hap2=-1;
        for (auto& hap_id: favored_haps) {
            for (auto& hap_id_prime: favored_haps) {
                if (hap_id_prime==hap_id) continue;
                tmp_set.clear();
                for (auto& element: hap_to_reads.at(hap_id)) tmp_set.emplace(element);
                for (auto& element: hap_to_reads.at(hap_id_prime)) tmp_set.emplace(element);
                if (tmp_set.size()!=n_reads) continue;
                hap1=hap_id; hap2=hap_id_prime;
                break;
            }
            if (hap1!=-1 && hap2!=-1) break;
        }
        if (hap1!=-1 && hap2!=-1) {
            for (i=0; i<n_reads; i++) {
                read_id=read_ids.at(i);
                min_hap_id=-1; min_weight=INT32_MAX;
                n_haps=neighbors.at(i).size();
                for (j=0; j<n_haps; j++) {
                    hap_id=neighbors.at(i).at(j);
                    weight=weights.at(i).at(j);
                    if (hap_id!=hap1 && hap_id!=hap2) edges_to_remove.emplace_back(read_id,hap_id);
                    else if (weight<min_weight) { min_weight=weight; min_hap_id=hap_id; }
                }
                for (j=0; j<n_haps; j++) {
                    hap_id=neighbors.at(i).at(j);
                    if (hap_id==hap1) {
                        if (min_hap_id!=hap1) edges_to_remove.emplace_back(read_id,hap_id);
                        else present_edges.emplace(read_id,hap_id);
                    }
                    else if (hap_id==hap2) {
                        if (min_hap_id!=hap2) edges_to_remove.emplace_back(read_id,hap_id);
                        else present_edges.emplace(read_id,hap_id);
                    }
                }
            }
            present_haps.emplace(hap1); present_haps.emplace(hap2);
            solved_samples.insert(sample_id);
            n_two_hap_samples++;
            return;
        }
    });
    // if (n_one_hap_samples+n_two_hap_samples!=0) cerr << "Solved " << to_string(n_one_hap_samples) << " one-hap samples and " << to_string(n_two_hap_samples) << " two-hap samples. Removed " << to_string(edges_to_remove.size()) << " edges. Set to one " << to_string(present_haps.size()) << " haplotypes and " << to_string(present_edges.size()) << " edges.\n";

    // Updating the transmap.
    // Remark: edge removal could be implemented much faster.
    for (auto& pair: edges_to_remove) remove_edge(pair.first,pair.second);
}


TransMap TransMap::solve_easy_samples_get_test_transmap(float n_weight, float d_weight) {
    const float DELTA = n_weight/d_weight;
    const float SMALL_WEIGHT = DELTA+1;
    const float LARGE_WEIGHT = SMALL_WEIGHT+DELTA;

    TransMap out;
    out.add_sample("sample1");
    out.add_read("s1r1"); out.add_read("s1r2"); out.add_read("s1r3");
    out.add_edge("sample1","s1r1"); out.add_edge("sample1","s1r2"); out.add_edge("sample1","s1r3");
    out.add_sample("sample2");
    out.add_read("s2r1"); out.add_read("s2r2"); out.add_read("s2r3");
    out.add_edge("sample2","s2r1"); out.add_edge("sample2","s2r2"); out.add_edge("sample2","s2r3");
    out.add_sample("sample3");
    out.add_read("s3r1"); out.add_read("s3r2"); out.add_read("s3r3");
    out.add_edge("sample3","s3r1"); out.add_edge("sample3","s3r2"); out.add_edge("sample3","s3r3");

    out.add_path("path1"); out.add_path("path2"); out.add_path("path3");
    out.add_edge("s1r1","path1",SMALL_WEIGHT); out.add_edge("s1r1","path2",LARGE_WEIGHT); out.add_edge("s1r1","path3",LARGE_WEIGHT);
    out.add_edge("s1r2","path1",LARGE_WEIGHT); out.add_edge("s1r2","path2",LARGE_WEIGHT); out.add_edge("s1r2","path3",LARGE_WEIGHT);
    out.add_edge("s1r3","path1",LARGE_WEIGHT); out.add_edge("s1r3","path2",LARGE_WEIGHT); out.add_edge("s1r3","path3",LARGE_WEIGHT);

    out.add_edge("s2r1","path2",SMALL_WEIGHT); out.add_edge("s2r1","path1",LARGE_WEIGHT); out.add_edge("s2r1","path3",LARGE_WEIGHT);
    out.add_edge("s2r2","path3",SMALL_WEIGHT); out.add_edge("s2r2","path1",LARGE_WEIGHT); out.add_edge("s2r2","path2",LARGE_WEIGHT);
    out.add_edge("s2r3","path1",LARGE_WEIGHT); out.add_edge("s2r3","path2",LARGE_WEIGHT); out.add_edge("s2r3","path3",LARGE_WEIGHT);

    out.add_edge("s3r1","path1",SMALL_WEIGHT); out.add_edge("s3r1","path2",LARGE_WEIGHT); out.add_edge("s3r1","path3",LARGE_WEIGHT);
    out.add_edge("s3r2","path1",LARGE_WEIGHT); out.add_edge("s3r2","path2",LARGE_WEIGHT); out.add_edge("s3r2","path3",LARGE_WEIGHT);
    out.add_edge("s3r3","path2",LARGE_WEIGHT); out.add_edge("s3r3","path3",LARGE_WEIGHT);
    return out;
}


float TransMap::get_edge_weight(float weight, float weight_quantum) {
    return weight_quantum!=0?(float)(round(weight/weight_quantum)*weight_quantum):weight;
}


void TransMap::clear_present_haps_edges() {
    present_haps.clear(); present_edges.clear();
}


}
