#include "SteinerTree.hpp"


namespace sv_merge {


SteinerTree::SteinerTree(TransMap& transmap, float hap_cost, float edge_cost_multiplier):
    solutions(),
    solution_samples(),
    solution_sample_weights()
{
    this->hap_cost=hap_cost;
    this->edge_cost_multiplier=edge_cost_multiplier;
    const int64_t N_HAPS = transmap.get_n_paths();

    size_t n_haps, n_reads;
    int64_t i, j, k;
    int64_t read_id, hap_id, hap1, hap2, n_solutions1, n_solutions2, n_solutions;
    float weight, weight1, weight2;
    vector<int64_t> read_ids;
    unordered_set<int64_t> tmp_set;
    vector<vector<int64_t>> neighbors;
    vector<vector<float>> weights;
    unordered_map<int64_t,float> solutions1_prime;
    unordered_map<int64_t,unordered_set<int64_t>> hap_to_reads;
    unordered_map<int64_t,vector<pair<int64_t,float>>> solutions1;
    unordered_map<int64_t,unordered_map<int64_t,float>> solutions2_prime;
    unordered_map<int64_t,unordered_map<int64_t,vector<pair<int64_t,float>>>> solutions2;

    transmap.update_first_of_type();

    // Enumerating solution-sample edges
    hap_to_reads.reserve(N_HAPS);
    solutions1.reserve(N_HAPS); solutions2.reserve(N_HAPS);
    solutions1_prime.reserve(N_HAPS); solutions2_prime.reserve(N_HAPS);
    n_solutions=0;
    transmap.for_each_sample([&](const string& sample_name, int64_t sample_id) {
        // Building read-hap and hap-read maps.
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
        read_ids.clear();
        // Enumerating one-hap solutions
        solutions1_prime.clear();
        for (auto& element: hap_to_reads) {
            if (element.second.size()==n_reads) solutions1_prime[element.first]=0;
        }
        for (i=0; i<n_reads; i++) {
            n_haps=neighbors.at(i).size();
            for (j=0; j<n_haps; j++) {
                hap_id=neighbors.at(i).at(j);
                if (solutions1_prime.contains(hap_id)) solutions1_prime.at(hap_id)+=weights.at(i).at(j);
            }
        }
        for (auto& element: solutions1_prime) {
            if (solutions1.contains(element.first)) solutions1.at(element.first).emplace_back(sample_id,element.second);
            else {
                vector<pair<int64_t,float>> new_vector;
                new_vector.emplace_back(sample_id,element.second);
                solutions1[element.first]=new_vector;
                n_solutions1++;
            }
        }
        solutions1_prime.clear();
        // Enumerating two-hap solutions
        solutions2_prime.clear();
        for (auto& element: hap_to_reads) {
            hap1=element.first;
            for (auto& element_prime: hap_to_reads) {
                hap2=element_prime.first;
                if (hap2<=hap1) continue;
                tmp_set.clear();
                for (auto& r: element.second) tmp_set.insert(r);
                for (auto& r: element_prime.second) tmp_set.insert(r);
                if (tmp_set.size()!=n_reads) continue;
                if (solutions2_prime.contains(hap1)) solutions2_prime[hap1][hap2]=0;
                else {
                    unordered_map<int64_t,float> new_map;
                    new_map[hap2]=0;
                    solutions2_prime[hap1]=new_map;
                }
            }
        }
        for (i=0; i<n_reads; i++) {
            n_haps=neighbors.at(i).size();
            for (j=0; j<n_haps; j++) {
                hap1=neighbors.at(i).at(j);
                if (!solutions2.contains(hap1)) continue;
                weight1=weights.at(i).at(j);
                for (k=j+1; k<n_haps; k++) {
                    hap2=neighbors.at(i).at(k);
                    if (!solutions2[hap1].contains(hap2)) continue;
                    weight2=weights.at(i).at(k);
                    solutions2_prime[hap1][hap2]+=weight1<weight2?weight1:weight2;
                }
            }
        }
        neighbors.clear(); weights.clear(); hap_to_reads.clear();
        for (auto& element1: solutions2_prime) {
            hap1=element1.first;
            if (solutions2.contains(hap1)) {
                for (auto& element2: element1.second) {
                    hap2=element2.first;
                    if (solutions2.at(hap1).contains(hap2)) solutions2[hap1][hap2].emplace_back(sample_id,element2.second);
                    else {
                        vector<pair<int64_t,float>> new_vector;
                        new_vector.emplace_back(sample_id,element2.second);
                        solutions2[hap1][hap2]=new_vector;
                        n_solutions2++;
                    }
                }
            }
            else {
                unordered_map<int64_t,vector<pair<int64_t,float>>> new_map;
                for (auto& element2: element1.second) {
                    hap2=element2.first;
                    vector<pair<int64_t,float>> new_vector;
                    new_vector.emplace_back(sample_id,element2.second);
                    new_map[hap2]=new_vector;
                    n_solutions2++;
                }
                solutions2[hap1]=new_map;
            }
        }
        solutions2_prime.clear();
    });
    n_solutions=n_solutions1+n_solutions2;
    cerr << "Found " << to_string(n_solutions1) << " one-hap candidate solutions (" << to_string((100.0*n_solutions1)/N_HAPS) << "% of haps) and " << to_string(n_solutions2) << " two-hap candidate solutions (" << to_string((100.0*n_solutions2)/(N_HAPS*N_HAPS)) << "% of hap pairs)\n";

    // Building output matrices
    solutions.clear(); solutions.reserve(n_solutions);
    solution_samples.clear(); solution_samples.reserve(n_solutions);
    solution_sample_weights.clear(); solution_sample_weights.reserve(n_solutions);
    i=-1;
    for (auto& element: solutions1) {
        i++;
        solutions.emplace_back(element.first,element.first);
        solution_samples.emplace_back(); solution_sample_weights.emplace_back();
        for (auto& element_prime: element.second) {
            solution_samples.at(i).emplace_back(element_prime.first);
            solution_sample_weights.at(i).emplace_back(element_prime.second);
        }
    }
    solutions1.clear();
    for (auto& element: solutions2) {
        for (auto& element_prime: element.second) {
            i++;
            solutions.emplace_back(element.first,element_prime.first);
            solution_samples.emplace_back(); solution_sample_weights.emplace_back();
            for (auto& element_prime_prime: element_prime.second) {
                solution_samples.at(i).emplace_back(element_prime_prime.first);
                solution_sample_weights.at(i).emplace_back(element_prime_prime.second);
            }
        }
    }
}


void SteinerTree::dense_bunch() {

}




}