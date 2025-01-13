#include <queue>

using std::priority_queue;

#include "SteinerTree.hpp"


namespace sv_merge {

SteinerTree::SteinerTree(TransMap& transmap):
    solutions(),
    solution_samples(),
    solution_sample_weights(),
    solution_sample_reads()
{
    size_t n_reads, n_solutions1, n_solutions2, n_haps_prime, n_solutions;
    int64_t i, j;
    int64_t read_id, hap_id, hap1, hap2;
    double weight1, weight2;
    vector<int64_t> read_ids;
    unordered_set<int64_t> tmp_set;
    vector<vector<int64_t>> neighbors;
    vector<vector<float>> weights;
    unordered_map<int64_t, unordered_set<int64_t>> hap_to_reads;  // hapID -> readID (for one sample)
    unordered_map<int64_t, double> solutions1_prime;  // hapID -> weight (for one sample)
    unordered_map<tuple<int64_t,int64_t>, double> solutions2_prime;  // (hapID1,hapID2) -> weight (for one sample)
    unordered_map<int64_t, vector<pair<int64_t,double>>> solutions1;  // hapID -> (sampleID,weight)
    unordered_map<tuple<int64_t,int64_t>, vector<pair<int64_t,double>>> solutions2;  // (hapID1,hapID2) -> (sampleID,weight)

    n_samples=transmap.get_n_samples();
    n_haps=transmap.get_n_paths();
    transmap.update_first_of_type();

    // Enumerating solution-sample edges
    hap_to_reads.reserve(n_haps);
    solutions1.reserve(n_haps); solutions2.reserve(n_haps);
    solutions1_prime.reserve(n_haps); solutions2_prime.reserve(n_haps);
    solution_sample_reads.reserve(n_haps);
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
            n_haps_prime=neighbors.at(i).size();
            for (j=0; j<n_haps_prime; j++) {
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
            n_haps_prime=neighbors.at(i).size();
            for (j=0; j<n_haps_prime; j++) {
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
                n_haps_prime=neighbors.at(i).size();
                weight1=-1; weight2=-1;
                for (j=0; j<n_haps_prime; j++) {
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
    cerr << "Found " << to_string(n_solutions1) << " one-hap candidate solutions (" << to_string((100.0*n_solutions1)/n_haps) << "% of haps) and " << to_string(n_solutions2) << " two-hap candidate solutions (" << to_string((100.0*n_solutions2)/(n_haps*n_haps)) << "% of hap pairs)\n";

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


void SteinerTree::build_sample_solutions() {
    size_t length;
    int64_t i, j;
    int64_t sample_id;
    double weight;

    if (!sample_solutions.empty()) return;
    sample_solutions.clear();
    sample_solutions.reserve(n_samples);
    for (i=0; i<n_solutions; i++) {
        length=solution_samples.at(i).size();
        for (j=0; j<length; j++) {
            sample_id=solution_samples.at(i).at(j);
            weight=solution_sample_weights.at(i).at(j);
            if (sample_solutions.contains(sample_id)) sample_solutions[sample_id].emplace_back(i,weight);
            else {
                vector<pair<int64_t,double>> new_vector;
                new_vector.emplace_back(i,weight);
                sample_solutions[sample_id]=new_vector;
            }
        }
    }
}


void SteinerTree::get_sample_solutions_histogram(vector<int64_t> out) {
    const size_t OUT_LENGTH = out.size();
    size_t length;

    build_sample_solutions();
    for (int64_t i=0; i<n_samples; i++) {
        length=sample_solutions.at(i).size();
        out[length>=OUT_LENGTH?OUT_LENGTH-1:length]++;
    }
}


void SteinerTree::build_hap_solutions() {
    int64_t i;
    int64_t hap_id;

    if (!hap_solutions.empty()) return;
    hap_solutions.clear();
    for (i=0; i<n_solutions; i++) {
        hap_id=solutions.at(i).first;
        if (hap_solutions.contains(hap_id)) hap_solutions.at(hap_id).emplace_back(i);
        else {
            vector<int64_t> new_vector;
            new_vector.emplace_back(i);
            hap_solutions[hap_id]=new_vector;
        }
        hap_id=solutions.at(i).second;
        if (hap_solutions.contains(hap_id)) hap_solutions.at(hap_id).emplace_back(i);
        else {
            vector<int64_t> new_vector;
            new_vector.emplace_back(i);
            hap_solutions[hap_id]=new_vector;
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
    selected_haps.clear();
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




// -------------------------------------------- GREEDY DENSE SOLUTION --------------------------------------------------

double SteinerTree::get_density(int64_t solution_id, unordered_set<int64_t>& covered_samples, float hap_cost, float edge_cost_multiplier) {
    const size_t length = solution_sample_weights.at(solution_id).size();
    const int64_t hap1 = solutions.at(solution_id).first;
    const int64_t hap2 = solutions.at(solution_id).second;
    int64_t i;
    int64_t sample_id;
    double d, denominator, numerator;

    denominator=0.0;
    if (!selected_haps.contains(hap1)) denominator+=hap_cost;
    if (hap2!=hap1 && !selected_haps.contains(hap2)) denominator+=hap_cost;
    d=0; numerator=0.0;
    for (i=0; i<length; i++) {
        sample_id=solution_samples.at(solution_id).at(i);
        if (!covered_samples.contains(sample_id)) {
            numerator+=1.0;
            d+=solution_sample_weights.at(solution_id).at(i);
        }
    }
    denominator+=d*edge_cost_multiplier;
    return numerator/denominator;
}


double SteinerTree::get_density_init(int64_t solution_id, float hap_cost, float edge_cost_multiplier) {
    double denominator = solutions.at(solution_id).first==solutions.at(solution_id).second?hap_cost:(2.0*hap_cost);
    double d = 0.0;
    for (auto& weight: solution_sample_weights.at(solution_id)) d+=weight;
    denominator+=d*edge_cost_multiplier;
    return ((double)solution_samples.at(solution_id).size())/denominator;
}


double SteinerTree::greedy_dense_solution(float hap_cost, float edge_cost_multiplier, int64_t mode, vector<int64_t>& n_samples_with_solutions) {
    int64_t i;
    vector<double> density;
    unordered_set<int64_t> covered_samples, solutions_to_update;
    priority_queue<pair<int64_t,double>, vector<pair<int64_t,double>>, Compare> queue;

    // Initializing reused data structures
    covered_samples.reserve(n_samples);
    build_sample_solutions();
    get_sample_solutions_histogram(n_samples_with_solutions);
    build_hap_solutions();
    density.reserve(n_solutions);
    for (i=0; i<n_solutions; i++) density.emplace_back(get_density_init(i,hap_cost,edge_cost_multiplier));

    // Greedy steps
    if (mode==0) {
        approximate_clear();
        greedy_dense_solution_impl(density,queue,covered_samples,solutions_to_update,hap_cost,edge_cost_multiplier);
    }
    else if (mode==1) {
        double objective_min;
        vector<bool> selected_solutions_min;
        vector<double> density_prime;
        vector<vector<bool>> selected_edges_min;
        unordered_set<int64_t> covered_samples_new, selected_haps_min;

        covered_samples_new.reserve(n_samples);
        objective_min=INT32_MAX;
        for (i=0; i<n_solutions; i++) {
            approximate_clear();
            density_prime=density;
            greedy_dense_solution_step(i,density_prime,covered_samples,covered_samples_new,solutions_to_update,false,queue,hap_cost,edge_cost_multiplier);
            greedy_dense_solution_impl(density_prime,queue,covered_samples,covered_samples_new,solutions_to_update,hap_cost,edge_cost_multiplier);
            if (objective<objective_min) {
                objective_min=objective;
                selected_haps_min=selected_haps;
                selected_solutions_min=selected_solutions;
                selected_edges_min=selected_edges;
            }
        }
        objective=objective_min;
        selected_haps=selected_haps_min;
        selected_solutions=selected_solutions_min;
        selected_edges=selected_edges_min;
    }
    return objective;
}


void SteinerTree::greedy_dense_solution_impl(vector<double>& density, priority_queue<pair<int64_t, double>, vector<pair<int64_t, double>>, Compare>& queue, unordered_set<int64_t>& covered_samples, unordered_set<int64_t>& covered_samples_new, unordered_set<int64_t>& solutions_to_update, float hap_cost, float edge_cost_multiplier) {
    int64_t i;
    int64_t solution_id;
    pair<int64_t,double> current_solution;

    covered_samples.clear(); covered_samples_new.clear();
    while (!queue.empty()) queue.pop();
    for (i=0; i<n_solutions; i++) {
        if (!selected_solutions.at(i)) queue.emplace(i,density.at(i));
    }
    while (covered_samples.size()<n_samples) {
        current_solution=queue.top(); queue.pop();
        solution_id=current_solution.first;
        if (current_solution.second!=density.at(solution_id)) continue;  // Skipping out-of-date elements
        greedy_dense_solution_step(solution_id,density,covered_samples,covered_samples_new,solutions_to_update,true,queue,hap_cost,edge_cost_multiplier);
    }
}


void SteinerTree::greedy_dense_solution_step(int64_t solution_id, vector<double>& density, unordered_set<int64_t>& covered_samples, unordered_set<int64_t>& covered_samples_new, unordered_set<int64_t>& solutions_to_update, bool update_queue, priority_queue<pair<int64_t,double>, vector<pair<int64_t,double>>, Compare> queue, float hap_cost, float edge_cost_multiplier) {
    bool new1, new2;
    size_t length;
    int64_t i;
    int64_t hap1, hap2, sample_id;
    double d, denominator, numerator;

    // Updating the approximation
    selected_solutions.at(solution_id)=true;
    hap1=solutions.at(solution_id).first; hap2=solutions.at(solution_id).second;
    if (!selected_haps.contains(hap1)) {
        objective+=hap_cost;
        selected_haps.insert(hap1);
        new1=true;
    }
    else new1=false;
    if (hap2!=hap1 && !selected_haps.contains(hap2)) {
        objective+=hap_cost;
        selected_haps.insert(hap2);
        new2=true;
    }
    else new2=false;
    covered_samples_new.clear();
    d=0; length=solution_samples.at(solution_id).size();
    for (i=0; i<length; i++) {
        sample_id=solution_samples.at(solution_id).at(i);
        if (covered_samples.contains(sample_id)) continue;
        covered_samples_new.insert(sample_id);
        selected_edges.at(solution_id).at(i)=true;
        d+=solution_sample_weights.at(solution_id).at(i);
    }
    objective+=d*edge_cost_multiplier;

    // Updating solution densities
    solutions_to_update.clear();
    for (auto& sample: covered_samples_new) {
        for (auto& solution: sample_solutions.at(sample)) {
            if (!selected_solutions.at(solution.first)) solutions_to_update.insert(solution.first);
        }
    }
    covered_samples.merge(covered_samples_new);
    covered_samples_new.clear();
    if (covered_samples.size()==n_samples) {
        // No need to update densities for the next iteration
        return;
    }
    if (new1) {
        for (auto& solution: hap_solutions.at(hap1)) {
            if (!selected_solutions.at(solution)) solutions_to_update.insert(solution);
        }
    }
    if (new2) {
        for (auto& solution: hap_solutions.at(hap2)) {
            if (!selected_solutions.at(solution)) solutions_to_update.insert(solution);
        }
    }
    for (auto& solution: solutions_to_update) {
        density.at(solution)=get_density(solution,covered_samples,hap_cost,edge_cost_multiplier);
        if (update_queue) queue.emplace(solution,density.at(solution));
    }
}




// ------------------------------------------------- SHORTEST PATHS ----------------------------------------------------

double SteinerTree::approximate_shortest_paths(bool minimize, float hap_cost, float edge_cost_multiplier, bool build_solution, vector<int64_t>& n_samples_with_solutions) {
    size_t length;
    int64_t i;
    int64_t sample_id, solution_id, solution_min, solution_max;
    double weight, weight_min, weight_max, distance, distance_min, distance_max;

    build_sample_solutions();
    get_sample_solutions_histogram(n_samples_with_solutions);
    approximate_clear();
    for (auto& element: sample_solutions) {
        sample_id=element.first;
        distance_min=INT32_MAX; solution_min=-1; weight_min=-1; distance_max=0; solution_max=-1; weight_max=-1;
        for (auto& solution: element.second) {
            solution_id=solution.first;
            weight=solution.second*edge_cost_multiplier;
            distance=weight+(solutions.at(solution_id).first==solutions.at(solution_id).second?hap_cost:(2.0*hap_cost));
            if (distance<distance_min) { distance_min=distance; weight_min=weight; solution_min=solution.first; }
            if (distance>distance_max) { distance_max=distance; weight_max=weight; solution_max=solution.first; }
        }
        if (minimize) { solution_id=solution_min; objective+=weight_min; }
        else { solution_id=solution_max; objective+=weight_max; }
        if (build_solution) {
            selected_solutions.at(solution_id)=true;
            length=solution_samples.at(solution_id).size();
            for (i=0; i<length; i++) {
                if (solution_samples.at(solution_id).at(i)==sample_id) {
                    selected_edges.at(solution_id).at(i)=true;
                    break;
                }
            }
        }
    }
    for (i=0; i<n_solutions; i++) {
        if (selected_solutions.at(i)) { selected_haps.insert(solutions.at(i).first); selected_haps.insert(solutions.at(i).second); }
    }
    objective+=((double)selected_haps.size())*hap_cost;
    return objective;
}



}