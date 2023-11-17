#include "read_optimizer.hpp"


namespace sv_merge{

void construct_read_model(
        const TransMap& transmap,
        int64_t sample_id,
        CpModelBuilder& model,
        ReadVariables& vars,
        vector<int64_t>& representatives
        ){

    transmap.for_each_read_of_sample(sample_id, [&](int64_t id) {
        // Construct indicator booleans which represent parent status of each node
        auto& p = vars.is_parent.emplace(id, model.NewBoolVar()).first->second;

        // Keep track of the total number of parents assigned
        // TODO: currently unused
        vars.cost_n += p;
    });

    unordered_map<int64_t, LinearExpr> child_edge_sums;
    unordered_map<int64_t, LinearExpr> parent_edge_sums;

    int64_t n_reads = 0;
    int64_t n_edges = 0;

    // TODO use BFS/connected components to enumerate edges instead, some edges won't be present in practice
    transmap.for_each_read_of_sample(sample_id, [&](int64_t id_a){
//        cerr << id_a << ' ' << transmap.get_node(id_a).name << '\n';

        transmap.for_each_read_of_sample(sample_id, [&](int64_t id_b){
            if (id_a == id_b){
                return;
            }

            auto [success,w] = transmap.try_get_edge_weight(id_a, id_b);

            if (not success){
                return;
            }

            // Choosing any edge starting with id_a implies that id_a is parent and id_b is child
            pair<int64_t,int64_t> p = {id_a, id_b};
            auto& e = vars.edges.emplace(p, model.NewBoolVar()).first->second;

//            cerr << '\t' <<  transmap.get_node(id_a).name << "->" << transmap.get_node(id_b).name << ' ' << w << '\n';

            // Keep track of which edges would make id_a a parent
            parent_edge_sums[id_a] += e;

            // Keep track of which edges would make id_b a child
            child_edge_sums[id_b] += e;

            // Cost of selecting an edge
            vars.cost_d += int64_t(w)*e;
            n_edges++;
        });

        n_reads++;
    });

    transmap.for_each_read_of_sample(sample_id, [&](int64_t id){
        auto& p = vars.is_parent.at(id);

        auto p_result = parent_edge_sums.find(id);

        if (p_result == parent_edge_sums.end()){
            representatives.push_back(id);
            return;
        }

        auto& p_sum = parent_edge_sums.at(id);
        auto& c_sum = child_edge_sums.at(id);

        // If there is at least one edge from id to another node, it is a parent node
        model.AddGreaterThan(p_sum, 0).OnlyEnforceIf(p);
        model.AddLessOrEqual(p_sum, 0).OnlyEnforceIf(Not(p));

        // Constraint that being a parent node excludes possibility of incoming child edges
        auto c = model.NewBoolVar();
        model.AddEquality(c_sum, 1).OnlyEnforceIf(c);
        model.AddEquality(c_sum, 0).OnlyEnforceIf(Not(c));
        model.AddImplication(p, Not(c));

        // Constraint that multiple edges cannot point to a child
        model.AddLessOrEqual(c_sum, 1);

        // Constraint that a node must be covered by an edge
        model.AddGreaterThan(c_sum + p_sum, 0);
    });

    cerr << "Number of reads: " << n_reads << '\n';
    cerr << "Number of edges: " << n_edges << '\n';
}


void parse_read_model_solution(
        const CpSolverResponse& response,
        const ReadVariables& vars,
        TransMap& transmap,
        vector<int64_t>& representatives
        ){

    cerr << '\n';

    if (response.status() == CpSolverStatus::OPTIMAL || response.status() == CpSolverStatus::FEASIBLE) {
        cerr << "Maximum of objective function: " << response.objective_value() << '\n';

        set<string> results;

        // Iterate the parent assignment variables and print them
        for (const auto& [id,var]: vars.is_parent){
            auto name = transmap.get_node(id).name;

            bool is_assigned = (SolutionIntegerValue(response, var) == 1);
            if (is_assigned){
                results.insert(name + ": parent");
                representatives.push_back(id);
            }
            else{
                results.insert(name + ": child");
            }
        }

        // Iterate the edges and remove any that are not part of the solution
        for (const auto& [e,var]: vars.edges){
            auto [a,b] = e;

            // Model uses a directed graph, but transmap uses undirected. We remove the edge in transmap if neither
            // direction was used by the solution
            if (a < b){
                auto var_reverse = vars.edges.at({b,a});
                bool f_is_assigned = (SolutionIntegerValue(response, var) == 1);
                bool r_is_assigned = (SolutionIntegerValue(response, var_reverse) == 1);

                if (not (f_is_assigned or r_is_assigned)){
                    transmap.remove_edge(a,b);
                }
            }
        }

        for (auto& item: results){
            cerr << item << '\n';
        }

    }
    else {
        cerr << "No solution found." << '\n';
    }

    cerr << '\n';
    cerr << "Statistics" << '\n';
    cerr << CpSolverResponseStats(response) << '\n';

}


void optimize_reads_with_d(TransMap& transmap, int64_t sample_id, vector<int64_t>& representatives){
    CpModelBuilder model;
    ReadVariables vars;

    construct_read_model(transmap, sample_id, model, vars, representatives);

    model.Minimize(vars.cost_d);

    cerr << '\n';

    const CpSolverResponse response = Solve(model.Build());

    parse_read_model_solution(response, vars, transmap, representatives);
}


void optimize_reads_with_d_and_n(TransMap& transmap, int64_t sample_id, vector<int64_t>& representatives){
    CpModelBuilder model;
    ReadVariables vars;

    construct_read_model(transmap, sample_id, model, vars, representatives);

    // First find one extreme of the pareto set
    model.Minimize(vars.cost_d);

    const CpSolverResponse response_d = Solve(model.Build());

    auto n_max = SolutionIntegerValue(response_d, vars.cost_n);
    auto d_min = SolutionIntegerValue(response_d, vars.cost_d);

    // Then find the other extreme of the pareto set
    model.Minimize(vars.cost_n);

    const CpSolverResponse response_n = Solve(model.Build());

    auto n_min = SolutionIntegerValue(response_n, vars.cost_n);
    auto d_max = SolutionIntegerValue(response_n, vars.cost_d);

    // Use pareto extremes to normalize the ranges of each objective and then jointly minimize distance from (0,0)
    auto d_norm = LinearExpr(vars.cost_d);
    d_norm *= n_min;
    auto n_norm = LinearExpr(vars.cost_n);
    n_norm *= d_min;

    auto d_square = model.NewIntVar({d_min*d_min*n_min, d_max*d_max*n_min});
    model.AddMultiplicationEquality(d_square,{d_norm,d_norm});

    auto n_square = model.NewIntVar({n_min*n_min*d_min, n_max*n_max*d_min});
    model.AddMultiplicationEquality(n_square,{n_norm,n_norm});

    model.Minimize(n_square + d_square);

    const CpSolverResponse response_n_d = Solve(model.Build());

    parse_read_model_solution(response_n_d, vars, transmap, representatives);
}



}