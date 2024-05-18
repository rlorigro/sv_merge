#include "path_optimizer.hpp"
#include <fstream>

using std::ofstream;


namespace sv_merge{

/**
 * Construct a model such that each read must be assigned to exactly one path and each sample must have at most 2 paths.
 * For the sake of modularity, no explicit call to CpModelBuilder::Minimize(vars.cost_d) is made here, it must be made
 * after constructing the model.
 * @param transmap - TransMap representing the relationship between samples, paths, and reads
 * @param model - model to be constructed
 * @param vars - CPSAT BoolVar and LinearExpr objects to be filled in
 */
void construct_d_model(const TransMap& transmap, CpModelBuilder& model, PathVariables& vars){
    // TODO: reserve Variables map data structures based on how many edges there are known to be in the transmap?

    // Read->path boolean assignment indicators (the basis for the rest of the indicators)
    transmap.for_each_read_to_path_edge([&](int64_t read_id, int64_t path_id, float weight){
        pair<int64_t,int64_t> p = {path_id, read_id};

        // Update the boolean var map, and keep a reference to the inserted BoolVar
        const auto result = vars.path_to_read.emplace(p,model.NewBoolVar());
        const auto& v = result.first->second;

        auto w = int64_t(weight);
        vars.cost_d += v*w;
    });

    transmap.for_each_sample([&](const string& name, int64_t sample_id){
        unordered_map<int64_t,LinearExpr> path_to_sample_sums;

        // Sample->read boolean indicators
        transmap.for_each_read_of_sample(sample_id, [&](int64_t read_id){
            LinearExpr readwise_vars;

            // Read->path boolean indicators
            transmap.for_each_path_of_read(read_id, [&](int64_t path_id){
                const BoolVar& var = vars.path_to_read.at({path_id, read_id});
                path_to_sample_sums[path_id] += var;
                readwise_vars += var;
            });

            // Enforce (sum(r) == 1). i.e. read must be assigned to exactly one path
            model.AddEquality(readwise_vars, 1);
        });

        LinearExpr s;

        for (const auto& [path_id, sum]: path_to_sample_sums){
            pair<int64_t,int64_t> p = {path_id, sample_id};
            const auto result = vars.path_to_sample.emplace(p, model.NewBoolVar());
            const auto& b = result.first->second;

            // Add boolean indicator to represent sample->path edge is active
            model.AddGreaterThan(sum, 0).OnlyEnforceIf(b);
            model.AddLessOrEqual(sum, 0).OnlyEnforceIf(Not(b));

            // Increment a sum which indicates the total number of paths used in this sample
            s += b;
        }

        // Enforce (s <= 2). i.e. "ploidy constraint": sample must be assigned at most 2 paths
        model.AddLessOrEqual(s, 2);
    });
}


/**
 * Construct a model such that each read must be assigned to exactly one path and each sample must have at most 2 paths.
 * For the sake of modularity, no explicit call to CpModelBuilder::Minimize(vars.cost_n + vars.cost_d) is made here, it
 * must be made after constructing the model.
 * @param transmap - TransMap representing the relationship between samples, paths, and reads
 * @param model - model to be constructed
 * @param vars - CPSAT BoolVar and LinearExpr objects to be filled in
 */
void construct_joint_n_d_model(const TransMap& transmap, CpModelBuilder& model, PathVariables& vars){
    // Start with the base model, only tracking total cost `d` and read/ploidy constraints
    construct_d_model(transmap, model, vars);

    // Keep track of whether each path is used
    transmap.for_each_path([&](const string& name, int64_t path_id){
        LinearExpr pathwise_reads;
        transmap.for_each_read_of_path(path_id, [&](int64_t read_id){
            pathwise_reads += vars.path_to_read.at({path_id,read_id});
        });

        const auto result = vars.path_indicators.emplace(path_id,model.NewBoolVar());
        const auto& p = result.first->second;

        // Implement p == (sum(r) > 0). i.e. indicator for whether path is used
        model.AddGreaterThan(pathwise_reads, 0).OnlyEnforceIf(p);
        model.AddLessOrEqual(pathwise_reads, 0).OnlyEnforceIf(Not(p));
    });

    for (const auto& [path_id,p]: vars.path_indicators){
        vars.cost_n += p;
    }
}


void optimize_d(const TransMap& transmap){
    CpModelBuilder model;
    PathVariables vars;

    construct_d_model(transmap, model, vars);

    model.Minimize(vars.cost_d);

    const CpSolverResponse response = Solve(model.Build());

    if (response.status() == CpSolverStatus::OPTIMAL || response.status() == CpSolverStatus::FEASIBLE) {
        cerr << "Maximum of objective function: " << response.objective_value() << '\n';

        set<string> results;

        // Iterate the path->sample assignment variables and print them
        for (const auto& [p,var]: vars.path_to_sample){
            auto [path_id, sample_id] = p;
            auto path_name = transmap.get_node(path_id).name;
            auto sample_name = transmap.get_node(sample_id).name;

            bool is_assigned = (SolutionIntegerValue(response, var) == 1);
            if (is_assigned){
                results.insert(sample_name + ": " + path_name);
            }
        }

        for (auto& item: results){
            cerr << item << '\n';
        }

    }
    else {
        cerr << "No solution found." << '\n';
    }

    cerr << "Statistics" << '\n';
    cerr << CpSolverResponseStats(response) << '\n';
}


void optimize_d_plus_n(const TransMap& transmap, int64_t d_coeff, int64_t n_coeff){
    CpModelBuilder model;
    PathVariables vars;

    construct_joint_n_d_model(transmap, model, vars);

    model.Minimize(d_coeff*vars.cost_d + n_coeff*vars.cost_n);

    const CpSolverResponse response = Solve(model.Build());

    if (response.status() == CpSolverStatus::OPTIMAL || response.status() == CpSolverStatus::FEASIBLE) {
        cerr << "Maximum of objective function: " << response.objective_value() << '\n';

        set<string> results;

        // Iterate the path->sample assignment variables and print them
        for (const auto& [p,var]: vars.path_to_sample){
            auto [path_id, sample_id] = p;
            auto path_name = transmap.get_node(path_id).name;
            auto sample_name = transmap.get_node(sample_id).name;

            bool is_assigned = (SolutionIntegerValue(response, var) == 1);
            if (is_assigned){
                results.insert(sample_name + ": " + path_name);
            }
        }

        for (auto& item: results){
            cerr << item << '\n';
        }

    }
    else {
        cerr << "No solution found." << '\n';
    }

    cerr << "Statistics" << '\n';
    cerr << CpSolverResponseStats(response) << '\n';
}


void parse_read_model_solution(const CpSolverResponse& response_n_d, const PathVariables& vars, const TransMap& transmap, path output_dir){
    // Open a file
    path out_path = output_dir/"solution.csv";
    ofstream file(out_path);

    if (not file.is_open() or not file.good()){
        throw runtime_error("ERROR: cannot write to file: " + out_path.string());
        return;
    }

    // Write header
    file << "sample,read,path" << '\n';

    // Print the results of the ILP by iterating all samples, all reads of each sample, and all read/path edges in the transmap
    if (response_n_d.status() == CpSolverStatus::OPTIMAL || response_n_d.status() == CpSolverStatus::FEASIBLE) {
        transmap.for_each_sample([&](const string& sample_name, int64_t sample_id){
            transmap.for_each_read_of_sample(sample_id, [&](int64_t read_id){
                transmap.for_each_path_of_read(read_id, [&](int64_t path_id){
                    const BoolVar& var = vars.path_to_read.at({path_id, read_id});
                    bool is_assigned = SolutionIntegerValue(response_n_d, var);
                    if (is_assigned){
                        file << sample_name << ',' << transmap.get_node(read_id).name << ',' << transmap.get_node(path_id).name << '\n';
                    }
                });
            });
        });
    }
}


void optimize_reads_with_d_and_n(TransMap& transmap, int64_t d_weight, int64_t n_weight, size_t n_threads, path output_dir){
    CpModelBuilder model;
    PathVariables vars;

    construct_joint_n_d_model(transmap, model, vars);

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

    SatParameters parameters;
//    parameters.set_log_search_progress(true);
    parameters.set_num_workers(n_threads);
    const CpSolverResponse response_n_d = SolveWithParameters(model.Build(), parameters);

//    const CpSolverResponse response_n_d = Solve(model.Build());

    parse_read_model_solution(response_n_d, vars, transmap, output_dir);

    // Write a log contatining the solutioninfo and responsestats
    path out_path = output_dir/"log_optimizer.txt";
    ofstream file(out_path);

    if (not file.is_open() or not file.good()){
        throw runtime_error("ERROR: cannot write to file: " + out_path.string());
        return;
    }

    file << "n_path_to_read_vars: " << vars.path_to_read.size() << '\n';
    file << CpSolverResponseStats(response_n_d) << '\n';
}



}