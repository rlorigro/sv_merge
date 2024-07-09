#include "path_optimizer_mathopt.hpp"
#include "bdsg/include/bdsg/internal/hash_map.hpp"
#include <fstream>

#include <map>

using std::map;

using std::ofstream;


namespace sv_merge{

/**
 * Construct a model such that each read must be assigned to exactly one path and each sample must have at most 2 paths.
 * For the sake of modularity, no explicit call to CpModelBuilder::Minimize(vars.cost_n + vars.cost_d) is made here, it
 * must be made after constructing the model.
 * @param transmap - TransMap representing the relationship between samples, paths, and reads
 * @param model - model to be constructed
 * @param vars - container to hold ORTools objects which are filled in and queried later after solving
 */
void construct_joint_n_d_model(const TransMap& transmap, Model& model, PathVariables& vars, bool integral, bool use_ploidy_constraint){
    // DEFINE: hap vars
    transmap.for_each_path([&](const string& hap_name, int64_t hap_id){\
        string name = "h" + std::to_string(hap_id);
        vars.haps.emplace(hap_id, model.AddVariable(0,1,integral,name));
    });

    transmap.for_each_sample([&](const string& sample_name, int64_t sample_id){
        transmap.for_each_read_of_sample(sample_id, [&](int64_t read_id){
            transmap.for_each_path_of_read(read_id, [&](int64_t hap_id) {
                string hap_name = transmap.get_node(hap_id).name;
                string read_name = transmap.get_node(read_id).name;

                // DEFINE: read-hap vars
                string r_h_name = "r" + std::to_string(read_id) + "h" + std::to_string(hap_id);
                auto result = vars.read_hap.emplace(std::make_pair(read_id, hap_id), model.AddVariable(0,1,integral,r_h_name));
                auto& r_h = result.first->second;

                // DEFINE: flow
                vars.read_flow[read_id] += r_h;

                auto [success, w] = transmap.try_get_edge_weight(read_id, hap_id);
                if (not success){
                    throw runtime_error("ERROR: edge weight not found for read-hap: " + std::to_string(read_id) + ", " + std::to_string(hap_id));
                }

                // OBJECTIVE: accumulate d cost sum
                vars.cost_d += w*r_h;

                // Do only once for each unique pair of sample-hap
                if (vars.sample_hap.find({sample_id, hap_id}) == vars.sample_hap.end()){
                    // DEFINE: sample-hap vars
                    string s_h_name = "s" + std::to_string(sample_id) + "h" + std::to_string(hap_id);
                    auto result2 = vars.sample_hap.emplace(std::make_pair(sample_id, hap_id),model.AddVariable(0,1,integral,s_h_name));
                    auto& s_h = result2.first->second;

                    // DEFINE: ploidy
                    vars.ploidy[sample_id] += s_h;

                    // CONSTRAINT: vsh <= vh (indicator for usage of hap w.r.t. sample-hap)
                    model.AddLinearConstraint(s_h <= vars.haps.at(hap_id));
                }

                // CONSTRAINT: vrh <= vsh (indicator for usage of read-hap, w.r.t. sample-hap)
                model.AddLinearConstraint(r_h <= vars.sample_hap.at({sample_id, hap_id}));
            });
        });
    });

    // CONSTRAINT: read assignment (flow)
    for (const auto& [read_id,f]: vars.read_flow){
        model.AddLinearConstraint(f == 1);
    }

    if (use_ploidy_constraint){
        // CONSTRAINT: ploidy
        for (const auto& [sample_id,p]: vars.ploidy){
            model.AddLinearConstraint(p <= 2);
        }
    }

    // OBJECTIVE: accumulate n cost sum
    for (const auto& [hap_id,h]: vars.haps){
        vars.cost_n += h;
    }

}


void parse_read_model_solution(const SolveResult& result_n_d, const PathVariables& vars, TransMap& transmap, path output_dir){
    // Open a file
    path out_path = output_dir/"solution.csv";
    ofstream file(out_path);

    if (not file.is_open() or not file.good()){
        throw runtime_error("ERROR: cannot write to file: " + out_path.string());
    }

    // Write header
    file << "sample,read,path" << '\n';

    unordered_set <pair <int64_t, int64_t> > to_be_removed;

    // Print the results of the ILP by iterating all samples, all reads of each sample, and all read/path edges in the transmap
    if (result_n_d.termination.reason == TerminationReason::kOptimal) {
        transmap.for_each_sample([&](const string& sample_name, int64_t sample_id){
            transmap.for_each_read_of_sample(sample_id, [&](int64_t read_id){
                transmap.for_each_path_of_read(read_id, [&](int64_t path_id){
                    const Variable& var = vars.read_hap.at({read_id, path_id});

                    if (var.is_integer()){
                        auto is_assigned = bool(int64_t(round(result_n_d.variable_values().at(var))));

                        if (is_assigned){
                            file << sample_name << ',' << transmap.get_node(read_id).name << ',' << transmap.get_node(path_id).name << '\n';
                        }
                        else{
                            // Delete all the edges that are not assigned (to simplify iteration later)
                            to_be_removed.emplace(read_id, path_id);
                        }
                    }
                    else{
                        throw runtime_error("ERROR: solution parsing not implemented for non-integer variables");
                    }
                });
            });
        });

        for (const auto& [read_id, path_id]: to_be_removed){
            transmap.remove_edge(read_id, path_id);
        }
    }
    else{
        cerr << "WARNING: cannot update transmap for non-optimal solution" << '\n';
        transmap = {};
    }
}


/**
 * Optimize the assignment of reads to a given number of paths
 * @param model
 * @param vars
 * @param transmap
 * @param solver_type
 * @param result_read_path_edges The set of read-path edges that are assigned in the solution (result)
 * @param n the given number of paths
 * @param n_threads
 * @param output_dir
 * @return
 */
double optimize_d_given_n(
        Model& model,
        PathVariables& vars,
        const SolverType& solver_type,
        unordered_set <pair <int64_t,int64_t> >& result_read_path_edges,
        size_t n,
        size_t n_threads,
        path output_dir
        ){

    cerr << "Optimizing d for n: " << n << '\n';

    result_read_path_edges.clear();

    SolveArguments args;
    args.parameters.threads = n_threads;

    auto constraint = model.AddLinearConstraint(vars.cost_n == n);
    model.Minimize(vars.cost_d);

    const absl::StatusOr<SolveResult> response = Solve(model, solver_type, args);

    // Immediately undo the constraint after solving
    model.DeleteLinearConstraint(constraint);

    const auto result = response.value();

    // Check if the first solution is feasible
    if (result.termination.reason != TerminationReason::kOptimal){
        cerr << "WARNING: no solution for d_min found: " << output_dir << '\n';
        return -1;
    }

    // Extract the d value
    double d = vars.cost_d.Evaluate(result.variable_values());

    // Print the n and d values of the solution
    cerr << "n: " << n << "\td: " << d << '\n';

    // Parse the solution
    for (const auto& [edge,var]: vars.read_hap){
        if (var.is_integer()){
            // Round the value to the nearest integer (0 or 1)
            auto is_assigned = bool(int64_t(round(result.variable_values().at(var))));
            if (is_assigned){
                result_read_path_edges.emplace(edge);
            }
        }
        else{
            throw runtime_error("ERROR: solution parsing not implemented for non-integer variables");
        }
    }

    // Append log line to output file which contains the result of each optimization
    path out_path = output_dir/"log_optimizer.txt";
    ofstream file(out_path, std::ios_base::app);

    if (not file.is_open() or not file.good()){
        throw runtime_error("ERROR: cannot write to file: " + out_path.string());
    }

    file << n << ',' << d << '\n';

    return d;
}


double get_normalized_distance(double d, double n, double n_min, double n_max, double d_min, double d_max){
    double n_range = n_max - n_min;
    double d_range = d_max - d_min;

    if (n_range == 0){
        n_range = 1;
    }

    if (d_range == 0){
        d_range = 1;
    }

    // Compute normalized distance from utopia point as "cost"
    double d_norm = (d - d_min)/(d_range);
    double n_norm = (n - n_min)/(n_range);

    return sqrt(d_norm*d_norm + n_norm*n_norm);
}


/**
 * Given the range of n values, use binary search to find the optimal n and d values
 * @param model
 * @param vars
 * @param solver_type
 * @param result_read_path_edges
 * @param n_min Result of optimization for n_min (NOT normalized!)
 * @param d_min Result of optimization for d_min (NOT normalized!)
 * @param n_max Result of optimizing for d_min, assuming convexity (NOT normalized!)
 * @param d_max Result of optimizing for n_min, assuming convexity (NOT normalized!)
 * @param n_threads
 * @param output_dir
 */
void solve_with_golden_search(
        Model& model,
        PathVariables& vars,
        const SolverType& solver_type,
        TransMap& transmap,
        int64_t n_min,
        int64_t n_max,
        double d_min,
        double d_max,
        size_t n_threads,
        path output_dir
){
    // Within this function the d_cost is now referred to as simply "distance" and the d_cost for a given n value is
    // referred to as "n_distance". 'd' is no loger short for "distance" but is used to refer to the d point which
    // is generated during the sectioning step of the golden section search algorithm, which yields a,b and interior
    // points c,d.
    int64_t n;
    int64_t n_distance;

    // With the range of n determined, we can now use binary search to find the optimal n and d values
    unordered_map<int64_t, double> results;
    unordered_map<int64_t, unordered_set <pair <int64_t,int64_t> > > result_edges_cache;

    // Add the extreme values to the results map
    results.emplace(n_min, get_normalized_distance(d_max, n_min, n_min, n_max, d_min, d_max));
    results.emplace(n_max, get_normalized_distance(d_min, n_max, n_min, n_max, d_min, d_max));

    // Set arbitrary limit on maximum iterations
    int64_t max_iter = 20;

    int64_t i = 0;
    auto a_i = n_min;
    auto b_i = n_max;

    double phi_inverse = (sqrt(5) - 1)/2;

    // Minimize d for each given n until it can be proven that the d value is optimal (left and right values are larger)
    while (i < max_iter){
        ///        c = b - (b - a) * invphi
        ///        d = a + (b - a) * invphi
        ///        if f(c) < f(d):
        ///            b = d
        ///        else:  # f(c) > f(d) to find the maximum
        ///            a = c
        cerr << "---- Iteration: " << i << '\n';

        int64_t c_i = floor(double(b_i) - (double(b_i - a_i) * phi_inverse));
        int64_t d_i = floor(double(a_i) + (double(b_i - a_i) * phi_inverse));

        double c_distance = -1;
        double d_distance = -1;

        auto c_result = results.find(c_i);
        if (c_result == results.end()){
            double c = optimize_d_given_n(model, vars, solver_type, result_edges_cache[c_i], c_i, n_threads, output_dir);
            c_distance = get_normalized_distance(c, c_i, n_min, n_max, d_min, d_max);
            results.emplace(c_i, c_distance);
        }
        else{
            c_distance = c_result->second;
        }

        auto d_result = results.find(d_i);
        if (d_result == results.end()){
            double d = optimize_d_given_n(model, vars, solver_type, result_edges_cache[d_i], d_i, n_threads, output_dir);
            d_distance = get_normalized_distance(d, d_i, n_min, n_max, d_min, d_max);
            results.emplace(d_i, d_distance);

        }
        else{
            d_distance = d_result->second;
        }

        cerr << "a: " << a_i << "\tc: " << c_i << "\td: " << d_i << "\tb: " << b_i << '\n';
        cerr << "a: " << results.at(a_i) << "\tc: " << c_distance << "\td: " << d_distance << "\tb: " << results.at(b_i) << '\n';

        if (c_distance < d_distance){
            // Prevent recomputing the same interval
            if (b_i == d_i){
                break;
            }

            b_i = d_i;
        }
        else {
            // Prevent recomputing the same interval
            if (a_i == c_i){
                break;
            }

            a_i = c_i;
        }

        i++;
    }

    int64_t c = floor(double(b_i) - (double(b_i - a_i) * phi_inverse));
    int64_t d = floor(double(a_i) + (double(b_i - a_i) * phi_inverse));

    double c_distance = results.at(c);
    double d_distance = results.at(d);

    // Compute the final minimum
    if (c_distance < d_distance){
        n = c;
        n_distance = c_distance;
    }
    else {
        n = d;
        n_distance = d_distance;
    }

    // Filter the transmap using edges of the solution
    transmap.for_each_sample([&](const string& sample_name, int64_t sample_id){
        transmap.for_each_read_of_sample(sample_id, [&](int64_t read_id){
            transmap.for_each_path_of_read(read_id, [&](int64_t path_id){
                if (result_edges_cache[n].find({read_id, path_id}) == result_edges_cache[n].end()){
                    transmap.remove_edge(read_id, path_id);
                }
            });
        });
    });

    cerr << "n: " << n << "\td: " << n_distance << '\n';
}


/** Instead of giving the joint/quadratic distance objective to the solver, assume that the joint distance as a function
 * of n and d is convex and use binary search across n to find the optimal n and d values.
 * @param transmap
 * @param n
 * @param n_threads
 * @param output_dir
 * @param solver_type
 * @param use_ploidy_constraint
 */
void optimize_reads_with_d_and_n_using_golden_search(
        TransMap& transmap,
        double d_weight,
        double n_weight,
        size_t n_threads,
        path output_dir,
        const SolverType& solver_type,
        bool use_ploidy_constraint
        ){
    Model model;
    SolveArguments args;
    PathVariables vars;

    args.parameters.threads = n_threads;

    double n_max = -1;
    double d_min = -1;
    double n_min = -1;
    double d_max = -1;
    double n = -1;
    double d = -1;

    bool integral = true;
    if (solver_type == SolverType::kPdlp or solver_type == SolverType::kGlop){
        integral = false;
    }

    construct_joint_n_d_model(transmap, model, vars, integral, use_ploidy_constraint);

    // First find one extreme of the pareto set (using a tie-breaker cost for the other objective)
    model.Minimize(vars.cost_d + vars.cost_n*1e-6);

    const absl::StatusOr<SolveResult> response_d = Solve(model, solver_type, args);
    const auto result_d = response_d.value();

    // Check if the first solution is feasible
    if (result_d.termination.reason != TerminationReason::kOptimal){
        cerr << "WARNING: no solution for d_min found: " << output_dir << '\n';
        transmap = {};
        return;
    }

    // Infer the n and d values of the D_MIN solution (ignoring tie-breaker cost)
    n_max = vars.cost_n.Evaluate(result_d.variable_values());
    d_min = vars.cost_d.Evaluate(result_d.variable_values());

    cerr << "n_max: " << n_max << "\td_min: " << d_min << '\n';

    // Then find the other extreme of the pareto set (using a tie-breaker cost for the other objective)
    model.Minimize(vars.cost_n + vars.cost_d*1e-6);

    const absl::StatusOr<SolveResult> response_n = Solve(model, solver_type, args);
    const auto result_n = response_n.value();

    // Check if the second solution is feasible
    if (result_n.termination.reason != TerminationReason::kOptimal){
        cerr << "WARNING: no solution for n_min found: " << output_dir << '\n';
        transmap = {};
        return;
    }

    // Infer the n and d values of the N_MIN solution (ignoring tie-breaker cost)
    n_min = vars.cost_n.Evaluate(result_n.variable_values());
    d_max = vars.cost_d.Evaluate(result_n.variable_values());

    cerr << "n_min: " << n_min << "\td_max: " << d_max << '\n';

    // Solve with golden search
    solve_with_golden_search(model, vars, solver_type, transmap, n_min, n_max, d_min, d_max, n_threads, output_dir);

    return;
}


void optimize_reads_with_d_and_n(
        TransMap& transmap,
        double d_weight,
        double n_weight,
        size_t n_threads,
        path output_dir,
        const SolverType& solver_type,
        bool use_ploidy_constraint
        ){

    Model model;
    SolveArguments args;
    PathVariables vars;

    args.parameters.threads = n_threads;

    double n_max = -1;
    double d_min = -1;
    double n_min = -1;
    double d_max = -1;
    double n = -1;
    double d = -1;

    bool integral = true;
    if (solver_type == SolverType::kPdlp or solver_type == SolverType::kGlop){
        integral = false;
    }

    construct_joint_n_d_model(transmap, model, vars, integral, use_ploidy_constraint);

    // First find one extreme of the pareto set (using a tie-breaker cost for the other objective)
    model.Minimize(vars.cost_d + vars.cost_n*1e-6);

    const absl::StatusOr<SolveResult> response_d = Solve(model, solver_type, args);
    const auto result_d = response_d.value();

    // Check if the first solution is feasible
    if (result_d.termination.reason != TerminationReason::kOptimal){
        cerr << "WARNING: no solution for d_min found: " << output_dir << '\n';
        transmap = {};
        return;
    }

    // Infer the n and d values of the D_MIN solution (ignoring tie-breaker cost)
    n_max = vars.cost_n.Evaluate(result_d.variable_values());
    d_min = vars.cost_d.Evaluate(result_d.variable_values());

    cerr << "n_max: " << n_max << "\td_min: " << d_min << '\n';

    // Then find the other extreme of the pareto set (using a tie-breaker cost for the other objective)
    model.Minimize(vars.cost_n + vars.cost_d*1e-6);

    const absl::StatusOr<SolveResult> response_n = Solve(model, solver_type, args);
    const auto result_n = response_n.value();

    // Check if the second solution is feasible
    if (result_n.termination.reason != TerminationReason::kOptimal){
        cerr << "WARNING: no solution for n_min found: " << output_dir << '\n';
        transmap = {};
        return;
    }

    // Infer the n and d values of the N_MIN solution (ignoring tie-breaker cost)
    n_min = vars.cost_n.Evaluate(result_n.variable_values());
    d_max = vars.cost_d.Evaluate(result_n.variable_values());

    cerr << "n_min: " << n_min << "\td_max: " << d_max << '\n';

    // Now that we have the range of n and d values, we can normalize the costs and construct the quadratic objective
    auto n_range = n_max - n_min;
    auto d_range = d_max - d_min;

    // Avoid division by zero
    n_range = (n_range == 0 ? 1 : n_range);
    d_range = (d_range == 0 ? 1 : d_range);

    Variable d_norm = model.AddContinuousVariable(0,1,"d");
    Variable n_norm = model.AddContinuousVariable(0,1,"n");

    model.AddLinearConstraint(d_norm == (vars.cost_d - d_min)/d_range);
    model.AddLinearConstraint(n_norm == (vars.cost_n - n_min)/n_range);

    model.Minimize(d_norm*d_norm*d_weight + n_norm*n_norm*n_weight);

    const absl::StatusOr<SolveResult> response_n_d = Solve(model, solver_type, args);
    const auto result_n_d = response_n_d.value();

    // Write a log containing the solutioninfo and responsestats
    path out_path = output_dir/"log_optimizer.txt";
    ofstream file(out_path);

    if (not file.is_open() or not file.good()){
        throw runtime_error("ERROR: cannot write to file: " + out_path.string());
        return;
    }

    file << "n_read_to_hap_vars: " << vars.read_hap.size() << '\n';
    file << result_n_d << '\n';

    // Check if the final solution is feasible
    if (result_n_d.termination.reason != TerminationReason::kOptimal){
        cerr << "WARNING: no solution for joint optimization found: " << output_dir << '\n';
        transmap = {};
        return;
    }

    // Infer the optimal n and d values of the joint solution
    n = vars.cost_n.Evaluate(result_n_d.variable_values());
    d = vars.cost_d.Evaluate(result_n_d.variable_values());

    cerr << "n: " << n << "\td: " << d << '\n';

    if (integral) {
        parse_read_model_solution(result_n_d, vars, transmap, output_dir);
    }
    else{
        cerr << "WARNING: solution parsing not implemented for non-integer variables" << '\n';
        transmap = {};
    }
}


}
