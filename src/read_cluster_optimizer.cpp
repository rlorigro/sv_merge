#include "read_cluster_optimizer.hpp"

#include <cmath>

using std::max;

namespace sv_merge{

void construct_cluster_model(
        const TransMap& transmap,
        int64_t sample_id,
        CpModelBuilder& model,
        ReadClusterVariables& vars
        ){

    vector<int64_t> read_ids;

    transmap.for_each_read_of_sample(sample_id, [&](int64_t id) {
        read_ids.emplace_back(id);
    });

    // Establish indicator variables for all pairs and build a linear expression for the total cost of pairs in clusters
    for (size_t i=0; i<read_ids.size(); i++){
        for (size_t j=i+1; j<read_ids.size(); j++) {
            auto a = read_ids[i];
            auto b = read_ids[j];

            auto [success,w] = transmap.try_get_edge_weight(a, b);

            if (not success){
//                w = max(transmap.get_sequence(a).size(), transmap.get_sequence(b).size());
                cerr << "not found" << '\n';
                continue;
            }

            pair<int64_t,int64_t> p = {a,b};
            auto& e = vars.edges.emplace(p,model.NewBoolVar()).first->second;

            cerr << "adding edge to model: " << transmap.get_node(a).name << ',' << transmap.get_node(b).name << ' ' << w << '\n';

            vars.cost_d += int64_t(w)*e + int64_t(-w)*e.Not();

//            vars.cost_n += e.Not();
        }
    }

    for (size_t i=0; i<read_ids.size(); i++){
        for (size_t j=i+1; j<read_ids.size(); j++) {
            for (size_t k=j+1; k<read_ids.size(); k++) {
                auto a = read_ids[i];
                auto b = read_ids[j];
                auto c = read_ids[k];

                auto result_ab = vars.edges.find({a,b});
                if (result_ab == vars.edges.end()){
                    continue;
                }

                auto result_bc = vars.edges.find({b,c});
                if (result_bc == vars.edges.end()){
                    continue;
                }

                auto result_ac = vars.edges.find({a,c});
                if (result_ac == vars.edges.end()){
                    continue;
                }

                auto& ab = result_ab->second;
                auto& bc = result_bc->second;
                auto& ac = result_ac->second;

                model.AddLessOrEqual(-ab + bc + ac, 1);
                model.AddLessOrEqual( ab - bc + ac, 1);
                model.AddLessOrEqual( ab + bc - ac, 1);
            }
        }
    }
}


void parse_read_model_solution(
        const CpSolverResponse& response,
        const ReadClusterVariables& vars,
        int64_t sample_id,
        TransMap& transmap,
        vector <vector<int64_t> >& clusters
        ){

    cerr << '\n';

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

    transmap.for_each_read_of_sample(sample_id, [&](int64_t id) {
        unvisited.emplace(id);
    });

    auto criteria = [](const HeteroNode& node){
        return (node.type == 'R');
    };

    clusters.clear();

    while (not unvisited.empty()){
        auto start_name = transmap.get_node(*unvisited.begin()).name;

        // Initialize new empty cluster
        clusters.emplace_back();

        // Update cluster with each connected component from the solution
        transmap.for_node_in_bfs(start_name, 0, criteria, [&](const HeteroNode& node, int64_t id){
            unvisited.erase(id);
            clusters.back().emplace_back(id);
        });
    }
}


void cluster_reads(TransMap& transmap, int64_t sample_id, int64_t d_weight, int64_t n_weight, vector <vector<int64_t> >& clusters){
    CpModelBuilder model;
    ReadClusterVariables vars;

    construct_cluster_model(transmap, sample_id, model, vars);

    model.Maximize(vars.cost_d);

//    // First find one extreme of the pareto set
//    model.Minimize(vars.cost_d);
//
//    const CpSolverResponse response_d = Solve(model.Build());
//
//    int64_t n_max = SolutionIntegerValue(response_d, vars.cost_n);
//    int64_t d_min = SolutionIntegerValue(response_d, vars.cost_d);
//
//    // Then find the other extreme of the pareto set
//    model.Minimize(vars.cost_n);
//
//    const CpSolverResponse response_n = Solve(model.Build());
//
//    int64_t n_min = SolutionIntegerValue(response_n, vars.cost_n);
//    int64_t d_max = SolutionIntegerValue(response_n, vars.cost_d);
//
//    cerr << '\n';
//    cerr << "n_max: " << n_max << '\n';
//    cerr << "n_min: " << n_min << '\n';
//    cerr << "d_max: " << d_max << '\n';
//    cerr << "d_min: " << d_min << '\n';
//
//    // Use pareto extremes to normalize the ranges of each objective and then jointly minimize distance from (0,0)
//    auto n_range = n_max - n_min;
//    n_range *= n_range;
//    auto d_range = d_max - d_min;
//    d_range *= d_range;
//
//    cerr << "d_max*d_max*n_range*2*d_weight*n_weight: " << d_max*d_max*n_range*2*d_weight*n_weight << '\n';
//    cerr << "int64: " << numeric_limits<int64_t>::max() << '\n';
//    cerr << "int32: " << numeric_limits<int32_t>::max() << '\n';
//
//    auto d_norm = model.NewIntVar({0, d_max*d_max*n_range*2*d_weight*n_weight});
//    auto n_norm = model.NewIntVar({0, d_max*d_max*n_range*2*d_weight*n_weight});
//
//    model.AddEquality(d_norm, (vars.cost_d-d_min));
//    model.AddEquality(n_norm, (vars.cost_n-n_min));
//
//    auto d_square = model.NewIntVar({0, d_max*d_max*n_range*2*d_weight*n_weight});
//    model.AddMultiplicationEquality(d_square,{d_norm,d_norm});
//
//    auto n_square = model.NewIntVar({0, d_max*d_max*n_range*2*d_weight*n_weight});
//    model.AddMultiplicationEquality(n_square,{n_norm,n_norm});
//
//    model.Minimize(n_square*d_range*n_weight + d_square*n_range*d_weight);
//
//    SatParameters parameters;
//    parameters.set_log_search_progress(true);
//    const CpSolverResponse response_n_d = SolveWithParameters(model.Build(), parameters);

    const CpSolverResponse response_n_d = Solve(model.Build());

    parse_read_model_solution(response_n_d, vars, sample_id, transmap, clusters);
//    throw runtime_error("DEBUG EXIT EARLY");
}



}