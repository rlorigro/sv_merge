#include "TransitiveMap.hpp"
#include "VectorHeteroGraph.hpp"
#include "bdsg/include/bdsg/internal/hash_map.hpp"

using sv_merge::HeteroNode;
using sv_merge::TransMap;

#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <string>
#include <vector>
#include <set>
#include <iostream>
#include <stdexcept>

using std::runtime_error;
using std::cerr;
using std::unordered_map;
using std::unordered_map;
using std::unordered_set;
using std::vector;
using std::string;
using std::pair;
using std::set;

#include "ortools/base/logging.h"
#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model.pb.h"
#include "ortools/sat/cp_model_solver.h"
#include "ortools/util/sorted_interval_list.h"

using operations_research::sat::CpModelBuilder;
using operations_research::Domain;
using operations_research::sat::IntVar;
using operations_research::sat::BoolVar;
using operations_research::sat::CpSolverResponse;
using operations_research::sat::CpSolverStatus;
using operations_research::sat::LinearExpr;



/// Data class to contain the variables of the ILP, int64_t IDs intended to be reused from the TransMap
class Variables{
public:
    // Edges where the left node is the parent
    unordered_map <pair<int64_t,int64_t>, BoolVar> edges;

    // Cost of all edges assigned
    LinearExpr cost_d;

    // Cost of adding parents
    LinearExpr cost_n;
};


int main(){
    TransMap transmap;

    transmap.add_sample("HG001");

    transmap.add_read("read_01");
    transmap.add_read("read_02");
    transmap.add_read("read_03");
    transmap.add_read("read_04");
//    transmap.add_read("read_05");

    transmap.add_edge("read_01", "HG001");
    transmap.add_edge("read_02", "HG001");
    transmap.add_edge("read_03", "HG001");
    transmap.add_edge("read_04", "HG001");
//    transmap.add_edge("read_05", "HG001");

    // All vs All read-read edges
    transmap.add_edge("read_01", "read_02", 2);
    transmap.add_edge("read_01", "read_03", 2);
    transmap.add_edge("read_01", "read_04", 1);
    transmap.add_edge("read_02", "read_03", 1);
    transmap.add_edge("read_02", "read_04", 2);
    transmap.add_edge("read_03", "read_04", 2);

    cerr << '\n';

    CpModelBuilder model;
    Variables vars;
    string sample_name = "HG001";

    vector<int64_t> read_ids;

    transmap.for_each_read_of_sample(sample_name, [&](const string& name, int64_t id) {
        read_ids.emplace_back(id);
    });

    // Establish indicator variables for all pairs and build a linear expression for the total cost of pairs in clusters
    for (size_t i=0; i<read_ids.size(); i++){
        for (size_t j=i+1; j<read_ids.size(); j++) {
            auto a = read_ids[i];
            auto b = read_ids[j];

            auto [success,w] = transmap.try_get_edge_weight(a, b);

            if (not success){
                cerr << "edge not found: " << a << ',' << b << '\n';
                continue;
            }

            pair<int64_t,int64_t> p = {a,b};
            auto& e = vars.edges.emplace(p,model.NewBoolVar()).first->second;

            vars.cost_d += int64_t(w)*e;
            vars.cost_n -= e;
        }
    }

    for (size_t i=0; i<read_ids.size(); i++){
        for (size_t j=i+1; j<read_ids.size(); j++) {
            for (size_t k=j+1; k<read_ids.size(); k++) {
                auto a = read_ids[i];
                auto b = read_ids[j];
                auto c = read_ids[k];

                auto& ab = vars.edges.at({a,b});
                auto& bc = vars.edges.at({b,c});
                auto& ac = vars.edges.at({a,c});

                model.AddLessOrEqual(-ab + bc + ac, 1);
                model.AddLessOrEqual( ab - bc + ac, 1);
                model.AddLessOrEqual( ab + bc - ac, 1);
            }
        }
    }

    int64_t d_weight = 1;
    int64_t n_weight = 1;

    // First find one extreme of the pareto set
    model.Minimize(vars.cost_d);

    const CpSolverResponse response_d = Solve(model.Build());

    int64_t n_max = SolutionIntegerValue(response_d, vars.cost_n);
    int64_t d_min = SolutionIntegerValue(response_d, vars.cost_d);

    // Then find the other extreme of the pareto set
    model.Minimize(vars.cost_n);

    const CpSolverResponse response_n = Solve(model.Build());

    int64_t n_min = SolutionIntegerValue(response_n, vars.cost_n);
    int64_t d_max = SolutionIntegerValue(response_n, vars.cost_d);

    // Use pareto extremes to normalize the ranges of each objective and then jointly minimize distance from (0,0)
    auto n_range = n_max - n_min;
    n_range *= n_range;
    auto d_range = d_max - d_min;
    d_range *= d_range;

    auto d_norm = model.NewIntVar({0, d_max*d_max*n_range*2*d_weight*n_weight});
    auto n_norm = model.NewIntVar({0, d_max*d_max*n_range*2*d_weight*n_weight});

    model.AddEquality(d_norm, (vars.cost_d-d_min));
    model.AddEquality(n_norm, (vars.cost_n-n_min));

    auto d_square = model.NewIntVar({0, d_max*d_max*n_range*2*d_weight*n_weight});
    model.AddMultiplicationEquality(d_square,{d_norm,d_norm});

    auto n_square = model.NewIntVar({0, d_max*d_max*n_range*2*d_weight*n_weight});
    model.AddMultiplicationEquality(n_square,{n_norm,n_norm});

    model.Minimize(n_square*d_range*n_weight + d_square*n_range*d_weight);

//    SatParameters parameters;
//    parameters.set_log_search_progress(true);
//    const CpSolverResponse response_n_d = SolveWithParameters(model.Build(), parameters);

    const CpSolverResponse response = Solve(model.Build());

    if (response.status() == CpSolverStatus::OPTIMAL || response.status() == CpSolverStatus::FEASIBLE) {
        cerr << "Maximum of objective function: " << response.objective_value() << '\n';

        for (auto& [edge, var]: vars.edges){
            const auto& a = transmap.get_node(edge.first).name;
            const auto& b = transmap.get_node(edge.second).name;

            bool is_assigned = SolutionIntegerValue(response, var);
            cerr << a << ',' << b << ' ' << SolutionIntegerValue(response, var) << '\n';

            if (not is_assigned){
                transmap.remove_edge(edge.first,edge.second);
            }
        }
    }
    else {
        cerr << "No solution found." << '\n';
    }

    cerr << '\n';
    cerr << "Statistics" << '\n';
    cerr << CpSolverResponseStats(response) << '\n';

    unordered_set<int64_t> unvisited;

    for (const auto& item: read_ids){
        unvisited.emplace(item);
    }

    auto criteria = [](const HeteroNode& node){
        return (node.type == 'R');
    };

    vector <vector<int64_t> > clusters;

    while (unvisited.size() > 0){
        auto start_name = transmap.get_node(*unvisited.begin()).name;

        // Initialize new empty cluster
        clusters.emplace_back();

        // Update cluster with each connected component from the solution
        transmap.for_node_in_bfs(start_name, 0, criteria, [&](const HeteroNode& node, int64_t id){
            unvisited.erase(id);
            clusters.back().emplace_back(id);
        });
    }

    int c = 0;
    for (const auto& cluster: clusters){
        cerr << c++ << '\n';
        for (const auto& id: cluster){
            cerr << '\t' << id << ' ' << transmap.get_node(id).name << '\n';
        }
    }

    return 0;
}
