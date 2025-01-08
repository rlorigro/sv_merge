#include <queue>

using std::priority_queue;

#include "SteinerTree.hpp"


namespace sv_merge {

SteinerTree::SteinerTree(TransMap& transmap, float hap_cost, float edge_cost_multiplier):
    solutions(),
    solution_samples(),
    solution_sample_weights(),
    solution_sample_reads()
{
    this->hap_cost=hap_cost;
    this->edge_cost_multiplier=edge_cost_multiplier;
    const int64_t N_HAPS = transmap.get_n_paths();

    size_t n_haps, n_reads, n_solutions1, n_solutions2, n_solutions;
    int64_t i, j, k;
    int64_t read_id, hap_id, hap1, hap2;
    double weight1, weight2;
    vector<int64_t> read_ids;
    unordered_set<int64_t> tmp_set;
    vector<vector<int64_t>> neighbors;
    vector<vector<float>> weights;
    unordered_map<int64_t, unordered_set<int64_t>> hap_to_reads;  // hap -> reads (for one sample)
    unordered_map<int64_t, double> solutions1_prime;  // hap1 -> weight (for one sample)
    unordered_map<tuple<int64_t,int64_t>, double> solutions2_prime;  // (hap1,hap2) -> weight (for one sample)
    unordered_map<int64_t, vector<pair<int64_t,double>>> solutions1;  // hap1 -> (sample,weight)
    unordered_map<tuple<int64_t,int64_t>, vector<pair<int64_t,double>>> solutions2;  // (hap1,hap2) -> (sample,weight)

    n_samples=transmap.get_n_samples();
    transmap.update_first_of_type();

    // Enumerating solution-sample edges
    hap_to_reads.reserve(N_HAPS);
    solutions1.reserve(N_HAPS); solutions2.reserve(N_HAPS);
    solutions1_prime.reserve(N_HAPS); solutions2_prime.reserve(N_HAPS);
    solution_sample_reads.reserve(N_HAPS);
    transmap.for_each_sample([&](const string& sample_name, int64_t sample_id) {
        // Building read-hap and hap-read maps for this sample
        read_ids.clear(); neighbors.clear(); weights.clear();
        i=-1;
        transmap.for_each_read_of_sample(sample_id, [&](int64_t read_id) {  // Reads are assumed to be enumerated in sorted order
            i++; read_ids.emplace_back(read_id);
            neighbors.emplace_back(); weights.emplace_back();
            transmap.for_each_path_of_read(read_id, [&](int64_t path_id) {  // Haps are assumed to be enumerated in sorted order
                auto [success, weight] = transmap.try_get_edge_weight(read_id,path_id);
                neighbors.at(i).emplace_back(path_id);
                weights.at(i).emplace_back(weight);
            });
        });
        n_reads=read_ids.size();
        if (n_reads==0) return;
        tmp_set.reserve(n_reads);
        hap_to_reads.clear();
        for (i=0; i<n_reads; i++) {
            read_id=read_ids.at(i);
            n_haps=neighbors.at(i).size();
            for (j=0; j<n_haps; j++) {
                hap_id=neighbors.at(i).at(j);
                if (hap_to_reads.contains(hap_id)) hap_to_reads.at(hap_id).insert(read_id);
                else {
                    unordered_set<int64_t> new_set;
                    new_set.insert(read_id);
                    hap_to_reads[hap_id]=new_set;
                }
            }
        }
        // Computing one-hap solutions and their weight for this sample
        solutions1_prime.clear();
        for (auto& element: hap_to_reads) {
            if (element.second.size()==n_reads) solutions1_prime[element.first]=0.0;
        }
        for (i=0; i<n_reads; i++) {
            read_id=read_ids.at(i);
            n_haps=neighbors.at(i).size();
            for (j=0; j<n_haps; j++) {
                hap_id=neighbors.at(i).at(j);
                if (!solutions1_prime.contains(hap_id)) continue;
                solutions1_prime.at(hap_id)+=weights.at(i).at(j);
                auto key = std::make_tuple(hap_id,hap_id,sample_id);
                if (solution_sample_reads.contains(key)) solution_sample_reads.at(key).insert(std::make_pair(read_id,hap_id));
                else {
                    set<pair<int64_t,int64_t>> new_set;
                    new_set.insert(std::make_pair(read_id,hap_id));
                    solution_sample_reads[key]=new_set;
                }
            }
        }
        for (auto& element: solutions1_prime) {
            if (solutions1.contains(element.first)) solutions1.at(element.first).emplace_back(sample_id,element.second);
            else {
                vector<pair<int64_t,double>> new_vector;
                new_vector.emplace_back(sample_id,element.second);
                solutions1[element.first]=new_vector;
                n_solutions1++;
            }
        }
        solutions1_prime.clear();
        // Computing two-hap solutions and their weight for this sample
        solutions2_prime.clear();
        for (auto& element: hap_to_reads) {
            hap1=element.first;
            for (auto& element_prime: hap_to_reads) {
                hap2=element_prime.first;
                if (hap2<=hap1) continue;
                tmp_set.clear();
                for (auto& r: element.second) tmp_set.insert(r);
                for (auto& r: element_prime.second) tmp_set.insert(r);
                if (tmp_set.size()==n_reads) solutions2_prime[std::make_pair(hap1,hap2)]=0.0;
            }
        }
        hap_to_reads.clear();
        for (auto& element: solutions2_prime) {
            hap1=std::get<0>(element.first); hap2=std::get<1>(element.first);
            for (i=0; i<n_reads; i++) {
                n_haps=neighbors.at(i).size();
                weight1=-1; weight2=-1;
                for (j=0; j<n_haps; j++) {
                    hap_id=neighbors.at(i).at(j);
                    if (hap_id==hap1) weight1=weights.at(i).at(j);
                    else if (hap_id==hap2) weight2=weights.at(i).at(j);
                }
                if (weight1==-1 && weight2==-1) throw runtime_error("ERROR: a pair of haplotypes, that was considered a solution for a sample, does not cover all the reads in the sample.");
                auto key = std::make_tuple(hap1,hap2,sample_id);
                if (weight1<weight2) {
                    solutions2_prime.at(element.first)+=weight1;
                    if (solution_sample_reads.contains(key)) solution_sample_reads.at(key).insert(std::make_pair(read_id,hap1));
                    else {
                        set<pair<int64_t,int64_t>> new_set;
                        new_set.insert(std::make_pair(read_id,hap1));
                        solution_sample_reads[key]=new_set;
                    }
                }
                else {
                    solutions2_prime.at(element.first)+=weight2;
                    if (solution_sample_reads.contains(key)) solution_sample_reads.at(key).insert(std::make_pair(read_id,hap2));
                    else {
                        set<pair<int64_t,int64_t>> new_set;
                        new_set.insert(std::make_pair(read_id,hap2));
                        solution_sample_reads[key]=new_set;
                    }
                }
            }
        }
        neighbors.clear(); weights.clear();
        for (auto& element: solutions2_prime) {
            hap1=std::get<0>(element.first); hap2=std::get<1>(element.first);
            auto key = std::make_tuple(hap1,hap2);
            if (solutions2.contains(key)) solutions2[key].emplace_back(sample_id,element.second);
            else {
                vector<pair<int64_t,double>> new_vector;
                new_vector.emplace_back(sample_id,element.second);
                solutions2[key]=new_vector;
            }
        }
        solutions2_prime.clear();
    });
    n_solutions1=solutions1.size(); n_solutions2=solutions2.size(); n_solutions=n_solutions1+n_solutions2;
    cerr << "Found " << to_string(n_solutions1) << " one-hap candidate solutions (" << to_string((100.0*n_solutions1)/N_HAPS) << "% of haps) and " << to_string(n_solutions2) << " two-hap candidate solutions (" << to_string((100.0*n_solutions2)/(N_HAPS*N_HAPS)) << "% of hap pairs)\n";

    // Building output matrices
    solutions.clear(); solutions.reserve(n_solutions);
    solution_sample_weights.clear(); solution_sample_weights.reserve(n_solutions);
    i=-1;
    for (auto& element: solutions1) {  // Order of solutions not important
        i++;
        solutions.emplace_back(element.first,element.first);
        solution_samples.emplace_back(); solution_sample_weights.emplace_back();
        for (auto& element_prime: element.second) {
            solution_samples.at(i).emplace_back(element_prime.first);
            solution_sample_weights.at(i).emplace_back(element_prime.second);
        }
    }
    solutions1.clear();
    for (auto& element: solutions2) {  // Order of solutions not important
        i++;
        solutions.emplace_back(std::get<0>(element.first),std::get<1>(element.first));
        solution_samples.emplace_back(); solution_sample_weights.emplace_back();
        for (auto& element_prime: element.second) {
            solution_samples.at(i).emplace_back(element_prime.first);
            solution_sample_weights.at(i).emplace_back(element_prime.second);
        }
    }
}


void SteinerTree::build_sample_solutions(bool build_weights) {
    size_t length;
    int64_t i, j;
    int64_t sample_id;
    float weight;

    if (!sample_solutions.empty() && (!build_weights || !sample_solution_weights.empty())) return;
    sample_solutions.clear();
    sample_solutions.reserve(n_samples);
    if (build_weights) {
        sample_solution_weights.clear();
        sample_solution_weights.reserve(n_samples);
    }
    for (i=0; i<n_solutions; i++) {
        length=solution_samples.at(i).size();
        for (j=0; j<length; j++) {
            sample_id=solution_samples.at(i).at(j);
            weight=solution_sample_weights.at(i).at(j);
            if (sample_solutions.contains(sample_id)) {
                sample_solutions[sample_id].emplace_back(i);
                if (build_weights) sample_solution_weights[sample_id].emplace_back(weight);
            }
            else {
                vector<int64_t> new_vector;
                new_vector.emplace_back(i);
                sample_solutions[sample_id]=new_vector;
                if (build_weights) {
                    vector<float> new_vector_prime;
                    new_vector_prime.emplace_back(weight);
                    sample_solution_weights[sample_id]=new_vector_prime;
                }
            }
        }
    }
}


void SteinerTree::approximate_clear() {
    size_t i, j;
    size_t length;

    objective=0.0;
    selected_solutions.clear(); selected_solutions.reserve(n_solutions);
    for (i=0; i<n_solutions; i++) selected_solutions.emplace_back(false);
    selected_edges.clear(); selected_edges.reserve(n_solutions);
    for (i=0; i<n_solutions; i++) {
        selected_edges.emplace_back();
        length=solution_samples.at(i).size();
        for (j=0; j<length; j++) selected_edges.at(i).emplace_back(false);
    }
}


double SteinerTree::optimize_d(bool minimize, bool build_solution) {

    double min, max;

    build_sample_solutions(true);
    approximate_clear();
    for (auto& [sample_id, list]: sample_solutions) {
        min=INT32_MAX; max=0;
        for (auto& solution: list) {
            ----->
        }
    }




}


double SteinerTree::approximate_shortest_paths() {
    size_t length;
    int64_t i;
    int64_t solution_id, sample_id;
    double d, distance;

    vector<double> solution_distance;
    unordered_map<int64_t,int64_t> sample_previous;
    unordered_map<int64_t,double> sample_distance;

    // Initializing data structures
    solution_distance.reserve(n_solutions);
    for (i=0; i<n_solutions; i++) solution_distance.emplace_back(solutions.at(i).first==solutions.at(i).second?hap_cost:(2*hap_cost));
    sample_distance.reserve(n_samples);
    for (i=0; i<n_samples; i++) sample_distance[i]=INT32_MAX;
    sample_previous.reserve(n_samples);

    // Shortest path iteration
    for (solution_id=0; solution_id<n_solutions; solution_id++) {
        distance=solution_distance.at(solution_id);
        length=solution_samples.at(solution_id).size();
        for (i=0; i<length; i++) {
            sample_id=solution_samples.at(solution_id).at(i);
            d=distance+solution_sample_weights.at(solution_id).at(i)*edge_cost_multiplier;
            if (d<sample_distance[sample_id]) {
                sample_distance[sample_id]=d;
                sample_previous[sample_id]=solution_id;
            }
        }
    }

    // Outputting
    approximate_clear();
    for (auto& element: sample_previous) {
        sample_id=element.first; solution_id=element.second;
        selected_solutions.at(solution_id)=true;
        length=solution_samples.at(solution_id).size();
        for (i=0; i<length; i++) {
            if (solution_samples.at(solution_id).at(i)!=sample_id) continue;
            selected_edges.at(solution_id).at(i)=true;
            objective+=solution_sample_weights.at(solution_id).at(i)*edge_cost_multiplier;
            break;
        }
    }
    for (i=0; i<n_solutions; i++) {
        if (selected_solutions.at(i)) objective+=solution_distance.at(solution_id);
    }
    return objective;
}


template<typename T>
class removable_priority_queue: public std::priority_queue<T, std::vector<T>, std::less<T>> {
  public:
      bool remove(const T& value) {
          auto iter = std::find(this->c.begin(), this->c.end(), value);
          if (iter==this->c.end()) return false;
          else if (iter==this->c.begin()) this->pop();
          else this->c.erase(iter);
          return true;
     }
};


double SteinerTree::approximate_dense_bunch() {
    size_t length, n_covered;
    int64_t i;
    int64_t sample_id, solution_id;
    double d, numerator, denominator;
    vector<double> density;
    unordered_set<int64_t> covered, covered_new, solutions_to_update;
    auto comparator = [&](int64_t i, int64_t j) { return density.at(i)<density.at(j); };
    removable_priority_queue<int64_t, vector<int64_t>, decltype(comparator)> queue(comparator);

    // Initializing data structures
    covered.reserve(n_samples); covered_new.reserve(n_samples);
    build_sample_solutions(false);
    density.reserve(n_solutions);
    i=-1;
    for (i=0; i<n_solutions; i++) {
        i++;
        denominator=solutions.at(i).first==solutions.at(i).second?hap_cost:(2*hap_cost);
        d=0;
        for (auto& weight: solution_sample_weights.at(i)) d+=weight;
        denominator+=d*edge_cost_multiplier;
        density.emplace_back(((double)solution_samples.at(i).size())/denominator);
    }
    for (i=0; i<n_solutions; i++) queue.push(i);

    // Greedy iteration
    approximate_clear(); n_covered=0;
    while (true) {
        solution_id=queue.top(); queue.pop();
        // Updating the approximation
        selected_solutions.at(solution_id)=true;
        objective+=solutions.at(solution_id).first==solutions.at(solution_id).second?hap_cost:(2*hap_cost);
        covered_new.clear();
        d=0; length=solution_samples.at(solution_id).size();
        for (i=0; i<length; i++) {
            sample_id=solution_samples.at(solution_id).at(i);
            if (covered.contains(sample_id)) continue;
            covered_new.insert(sample_id);
            selected_edges.at(solution_id).at(i)=true;
            d+=solution_sample_weights.at(solution_id).at(i);
        }
        objective+=d*edge_cost_multiplier;
        n_covered+=covered_new.size();
        if (n_covered==n_samples) break;
        // Updating the priority queue
        solutions_to_update.clear();
        for (auto& sample: covered_new) {
            for (auto& solution: sample_solutions.at(sample)) {
                if (selected_solutions.at(solution)) continue;
                solutions_to_update.insert(solution);
                queue.remove(solution);
            }
        }
        for (auto& solution: solutions_to_update) {
            denominator=solutions.at(solution).first==solutions.at(solution).second?hap_cost:(2*hap_cost);
            length=solution_sample_weights.at(solution).size();
            d=0; numerator=0;
            for (i=0; i<length; i++) {
                sample_id=solution_samples.at(solution).at(i);
                if (!covered.contains(sample_id) && !covered_new.contains(sample_id)) {
                    numerator++;
                    d+=solution_sample_weights.at(solution).at(i);
                }
            }
            denominator+=d*edge_cost_multiplier;
            density.at(solution)=((double)numerator)/denominator;
            queue.push(solution);
        }
        covered.merge(covered_new);
    }
    return objective;
}


void SteinerTree::parse_approximation(TransMap& transmap, path output_dir) {
    size_t length;
    int64_t i, j;
    int64_t hap1, hap2, sample_id;
    string sample_name;
    ofstream file;
    set<pair<int64_t,int64_t>> to_be_kept, to_be_removed;

    // Writing the solution file
    if (!output_dir.empty()) {
        path out_path = output_dir/"solution.csv";
        file.open(out_path);
        if (!file.is_open() or !file.good()) throw runtime_error("ERROR: cannot write to file: " + out_path.string());
        file << "sample,read,path" << '\n';
    }
    for (i=0; i<n_solutions; i++) {
        if (!selected_solutions.at(i)) continue;
        hap1=solutions.at(i).first; hap2=solutions.at(i).second;
        length=solution_samples.at(i).size();
        for (j=0; j<length; j++) {
            if (!selected_edges.at(i).at(j)) continue;
            sample_id=solution_samples.at(i).at(j);
            sample_name=transmap.get_node(sample_id).name;
            for (const auto& [read_id,path_id]: solution_sample_reads.at(std::make_tuple(hap1,hap2,sample_id))) {
                if (!output_dir.empty()) file << sample_name << ',' << transmap.get_node(read_id).name << ',' << transmap.get_node(path_id).name << '\n';
                to_be_kept.emplace(read_id,path_id);
            }
        }
    }
    if (!output_dir.empty()) file.close();

    // Updating the transmap
    transmap.for_each_read_id([&](int64_t read_id) {
        transmap.for_each_path_of_read(read_id, [&](int64_t path_id){
            const auto key = std::make_pair(read_id,path_id);
            if (!to_be_kept.contains(key)) to_be_removed.emplace(key);
        });
    });
    to_be_kept.clear();
    for (const auto& [read_id,path_id]: to_be_removed) transmap.remove_edge(read_id,path_id);
}



}