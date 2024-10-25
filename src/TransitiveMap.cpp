#include "TransitiveMap.hpp"
#include "gaf.hpp"
#include <fstream>

using std::ofstream;


namespace sv_merge{

TransMap::TransMap():
        sample_node_name("sample_node"),
        read_node_name("read_node"),
        path_node_name("path_node"),
        variant_node_name("variant_node")
{
    // Source nodes use lower case types to avoid being confused with the types they point to
    graph.add_node(sample_node_name, 's');
    graph.add_node(read_node_name, 'r');
    graph.add_node(path_node_name, 'p');
    graph.add_node(variant_node_name, 'v');
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


void TransMap::add_read_with_move(string& name, string& sequence){
    if (name == read_node_name or name == sample_node_name or name == path_node_name or name == variant_node_name){
        throw runtime_error("ERROR: cannot add node with preset node name: " + name);
    }

    graph.add_node(name, 'R');
    graph.add_edge(read_node_name, name, 0);

    sequences.emplace(graph.name_to_id(name), std::move(sequence));
}


void TransMap::add_read_with_move(string& name, string& sequence, bool is_reverse){
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


int64_t TransMap::get_n_reads() const {
    return graph.get_edge_count(graph.name_to_id(read_node_name));
}


int64_t TransMap::get_n_samples() const {
    return graph.get_edge_count(graph.name_to_id(sample_node_name));
}


int64_t TransMap::get_n_paths() const {
    return graph.get_edge_count(graph.name_to_id(path_node_name));
}


bool TransMap::are_edges_distinct() const {
    return graph.are_edges_distinct();
}


/**
 * Currently implemented as a quadratic scan, probably too slow for large cohorts.
 */
void TransMap::compress(float weight_quantum, uint64_t mode) {
    const int64_t n_reads = get_n_reads();
    const int64_t n_samples = get_n_samples();

    size_t length;
    int64_t i, j, k;
    int64_t read_id, n_clusters;
    vector<bool> is_redundant;
    vector<int64_t> node_ids, cluster_representative, cluster_size;
    vector<string> sample_names;
    vector<vector<int64_t>> neighbors, compared_weights;
    vector<vector<float>> weights;

    // Making sure that the neighbors of all nodes lie in the same global order
    graph.sort_adjacency_lists();

    // Collecting all read-haplotype edges
    neighbors.reserve(n_reads);
    for (i=0; i<n_reads; i++) neighbors.emplace_back();
    weights.reserve(n_reads);
    for (i=0; i<n_reads; i++) weights.emplace_back();
    if (weight_quantum!=0) {
        compared_weights.reserve(n_reads);
        for (i=0; i<n_reads; i++) compared_weights.emplace_back();
    }
    node_ids.reserve(n_reads);
    i=-1;
    for_each_read([&](const string& name, int64_t read_id) {
        // The order in which reads are enumerated is not important
        node_ids.emplace_back(read_id);
        i++;
        graph.for_each_neighbor_of_type(read_id,'P',[&](int64_t path_id) {
            // For every read, its neighbors lie in the same global order.
            auto [success, weight] = try_get_edge_weight(read_id, path_id);
            neighbors.at(i).emplace_back(path_id);
            weights.at(i).emplace_back(weight);
            if (weight_quantum!=0) compared_weights.at(i).emplace_back((int64_t)floor(weight/weight_quantum));
        });
    });


/*
    cerr << "Read-hap weights before compression: \n";
    for (i=0; i<n_reads; i++) {
        read_id=node_ids.at(i);
        cerr << "read=" << to_string(read_id) << " weights=";
        for (j=0; j<weights.at(i).size(); j++) cerr << "(" << to_string(neighbors.at(i).at(j)) << "," << to_string(weights.at(i).at(j)) << "), ";
        cerr << "\n";
    }
*/

    cerr << "HG01175 before compression: \n";
    for_each_sample([&](const string& sample_name, int64_t sample_id) {
        if (sample_name!="HG01175") return;
        for_each_read_of_sample(sample_name, [&](const string& read_name, int64_t read_id) {
            for_each_path_of_read(read_id, [&](int64_t path_id) {
                auto [success, weight] = try_get_edge_weight(read_id, path_id);
                string path_name = get_node(path_id).name;
                cerr << sample_name << "," << read_name << "," << path_name << "," << to_string(weight) << '\n';
            });
        });
    });






    // Clustering reads; adding sample-read edges; removing redundant reads; computing new read-hap weights.


vector<int64_t> representative;
representative.reserve(n_reads);
for (i=0; i<n_reads; i++) representative.emplace_back();


    is_redundant.reserve(n_reads);
    for (i=0; i<n_reads; i++) is_redundant.emplace_back(false);
    if (mode==3) {
        cluster_size.reserve(n_reads);
        for (i=0; i<n_reads; i++) cluster_size.emplace_back(0);
    }
    n_clusters=0;
    for (i=0; i<n_reads; i++) {
        if (is_redundant.at(i)) continue;
        n_clusters++;
        if (mode==3) cluster_size.at(i)=1;
        read_id=node_ids.at(i); length=neighbors.at(i).size();


representative.at(i)=read_id;


        for (j=i+1; j<n_reads; j++) {
            if (is_redundant.at(j) || neighbors.at(j)!=neighbors.at(i) || (weight_quantum==0 && weights.at(j)!=weights.at(i)) || (weight_quantum!=0 && compared_weights.at(j)!=compared_weights.at(i))) continue;
            is_redundant.at(j)=true;

representative.at(j)=read_id;

            if (mode==3) cluster_size.at(i)++;
            for (k=0; k<length; k++) {
                switch (mode) {
                    case 0: weights.at(i).at(k)=std::max(weights.at(i).at(k),weights.at(j).at(k)); break;
                    case 1: weights.at(i).at(k)=std::min(weights.at(i).at(k),weights.at(j).at(k)); break;
                    case 2: weights.at(i).at(k)=weights.at(i).at(k)+weights.at(j).at(k); break;
                    case 3: weights.at(i).at(k)=weights.at(i).at(k)+weights.at(j).at(k); break;
                }
            }
            for_each_sample_of_read(node_ids.at(j),[&](int64_t sample_id) {
                graph.add_edge(sample_id,read_id,0);  // Overwrites existing edge if present. Updates both sides of the edge.
            });
            remove_node(node_ids.at(j));
        }
    }
    cerr << "Read clusters: " << to_string(n_clusters) << " N. reads: " << to_string(n_reads) << "\n";
    cerr << "Read-hap weights after compression: \n";
    for (i=0; i<n_reads; i++) {
        read_id=node_ids.at(i);
        cerr << "read=" << to_string(read_id) << " representative=" << to_string(representative.at(i)) << " weights=";
        for (j=0; j<weights.at(i).size(); j++) cerr << "(" << to_string(neighbors.at(i).at(j)) << "," << to_string(weights.at(i).at(j)) << "), ";
        cerr << "\n";
    }

    // Updating read-hap weights
    for (i=0; i<n_reads; i++) {
        if (is_redundant.at(i)) continue;
        if (mode==3) {
            length=weights.at(i).size();
            for (j=0; j<length; j++) weights.at(i).at(j)/=cluster_size.at(i);
        }
        read_id=node_ids.at(i); length=neighbors.at(i).size();
        for (j=0; j<length; j++) graph.update_edge_weight(read_id,neighbors.at(i).at(j),weights.at(i).at(j));
    }
    neighbors.clear(); weights.clear(); compared_weights.clear(); is_redundant.clear(); cluster_size.clear(); node_ids.clear();

    // Making sure that the neighbors of all nodes lie in the same global order after read compression
    graph.sort_adjacency_lists();

    // Collecting all sample-compressedRead edges
    sample_names.reserve(n_samples); node_ids.reserve(n_samples); neighbors.reserve(n_samples);
    i=-1;
    for_each_sample([&](string sample_name, int64_t sample_id) {
        i++; sample_names.emplace_back(sample_name); node_ids.emplace_back(sample_id); neighbors.emplace_back();
        for_each_read_of_sample(sample_id, [&](int64_t read_id) {
            // For every sample, its compressed read neighbors lie in the same global order.
            neighbors.at(i).emplace_back(read_id);
        });
    });

    // Removing redundant samples
    sample_to_compressed_sample.clear();
    is_redundant.clear(); is_redundant.reserve(n_samples);
    for (i=0; i<n_samples; i++) is_redundant.emplace_back(false);
    n_clusters=0;
    for (i=0; i<n_samples; i++) {
        if (is_redundant.at(i)) continue;
        n_clusters++;
        for (j=i+1; j<n_samples; j++) {
            if (is_redundant.at(j) || neighbors.at(j)!=neighbors.at(i)) continue;
            is_redundant.at(j)=true;
            sample_to_compressed_sample.emplace(sample_names.at(j),sample_names.at(i));
            remove_node(node_ids.at(j));
        }
    }
    neighbors.clear(); is_redundant.clear(); node_ids.clear();
    cerr << "Sample clusters: " << to_string(n_clusters) << " N. samples: " << to_string(n_samples) << "\n";




    cerr << "Elements in transmap after compression: \n";
    int64_t x = 0;
    for_each_sample([&](string sample_name, int64_t sample_id) { x++; });
    cerr << "Transmap samples: " << to_string(x) << "\n";
    x=0;
    for_each_read([&](string read_name, int64_t read_id) { x++; });
    cerr << "Transmap reads: " << to_string(x) << "\n";
    x=0;
    for_each_path([&](string path_name, int64_t path_id) { x++; });
    cerr << "Transmap paths: " << to_string(x) << "\n";




    cerr << "HG01175 after compression: \n";
    for_each_sample([&](const string& sample_name, int64_t sample_id) {
        if (sample_name!="HG01175") return;
        for_each_read_of_sample(sample_name, [&](const string& read_name, int64_t read_id) {
            for_each_path_of_read(read_id, [&](int64_t path_id) {
                auto [success, weight] = try_get_edge_weight(read_id, path_id);
                string path_name = get_node(path_id).name;
                cerr << sample_name << "," << read_name << "," << path_name << "," << to_string(weight) << '\n';
            });
        });
    });


}


void TransMap::decompress_samples() {
    for (auto& pair: sample_to_compressed_sample) {
        add_sample(pair.first);
        for_each_read_of_sample(pair.second, [&](const string& read_name, int64_t read_id) {
            add_edge(pair.first,read_name,0);  // Updates both sides of the edge
        });
    }
    sample_to_compressed_sample.clear();
}


}
