#include "bdsg/include/bdsg/internal/hash_map.hpp"
#include "path_optimizer_mathopt.hpp"
//#include "gurobi_manager.hpp"
#include "Timer.hpp"

#include <fstream>
#include <thread>
#include <map>

using std::map;
using std::ofstream;


namespace sv_merge{


string termination_reason_to_string(const TerminationReason& reason){
    if (reason ==  TerminationReason::kOptimal){
        return "Optimal";
    }
    else if (reason ==  TerminationReason::kInfeasible){
        return "Infeasible";
    }
    else if (reason ==  TerminationReason::kUnbounded){
        return "Unbounded";
    }
    else if (reason ==  TerminationReason::kInfeasibleOrUnbounded){
        return "InfeasibleOrUnbounded";
    }
    else if (reason ==  TerminationReason::kImprecise){
        return "Imprecise";
    }
    else if (reason ==  TerminationReason::kFeasible){
        return "Feasible";
    }
    else if (reason ==  TerminationReason::kNoSolutionFound){
        return "NoSolutionFound";
    }
    else if (reason ==  TerminationReason::kNumericalError){
        return "NumericalError";
    }
    else if (reason ==  TerminationReason::kOtherError){
        return "OtherError";
    }
    else{
        throw runtime_error("ERROR: unrecongnized TerminationReason");
    }
}


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


/**
 * Construct a model such that each read must be assigned to exactly one path and each sample must have at most 2 paths.
 * For the sake of modularity, no explicit call to CpModelBuilder::Minimize(vars.cost_n + vars.cost_d) is made here, it
 * must be made after constructing the model.
 * @param transmap - TransMap representing the relationship between samples, paths, and reads
 * @param model - model to be constructed
 * @param vars - container to hold ORTools objects which are filled in and queried later after solving
 */
void construct_r_model(const TransMap& transmap, Model& model, ReadVariables& vars, bool integral){
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

                // Do only once for each unique pair of sample-hap
                if (vars.sample_hap.find({sample_id, hap_id}) == vars.sample_hap.end()){
                    // DEFINE: sample-hap vars
                    string s_h_name = "s" + std::to_string(sample_id) + "h" + std::to_string(hap_id);
                    auto result2 = vars.sample_hap.emplace(std::make_pair(sample_id, hap_id),model.AddVariable(0,1,integral,s_h_name));
                    auto& s_h = result2.first->second;

                    // DEFINE: ploidy
                    vars.ploidy[sample_id] += s_h;
                }

                // CONSTRAINT: vrh <= vsh (indicator for usage of read-hap, w.r.t. sample-hap)
                model.AddLinearConstraint(r_h <= vars.sample_hap.at({sample_id, hap_id}));
            });
        });
    });

    // CONSTRAINT: read assignment (flow)
    for (const auto& [read_id,f]: vars.read_flow){
        auto result = vars.reads.emplace(read_id, model.AddVariable(0,1,integral,"r" + std::to_string(read_id) + "_retained"));
        auto& r = result.first->second;

        // If this read is "on" then the sum of its read-hap variables must be 1, otherwise they are 0
        model.AddLinearConstraint(f == r);

        // OBJECTIVE: accumulate r cost sum
        // This will be maximized to retain the most reads while maintaining feasibility
        vars.cost_r += r;
    }

    // CONSTRAINT: ploidy (the point of this model)
    for (const auto& [sample_id,p]: vars.ploidy){
        model.AddLinearConstraint(p <= 2);
    }
}


void parse_read_model_solution(
    const TerminationReason& termination_reason,
    const VariableMap<double>& result_var_map,
    const PathVariables& vars,
    TransMap& transmap,
    path output_dir)
{

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
    if (termination_reason == TerminationReason::kOptimal) {
        transmap.for_each_sample([&](const string& sample_name, int64_t sample_id){
            transmap.for_each_read_of_sample(sample_id, [&](int64_t read_id){
                transmap.for_each_path_of_read(read_id, [&](int64_t path_id){
                    const Variable& var = vars.read_hap.at({read_id, path_id});

                    if (var.is_integer()){
                        auto is_assigned = bool(int64_t(round(result_var_map.at(var))));

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


void parse_read_model_solution(
    const SolveResult& result_n_d,
    const PathVariables& vars,
    TransMap& transmap,
    path output_dir)
{
    parse_read_model_solution(
        result_n_d.termination.reason,
        result_n_d.variable_values(),
        vars,
        transmap,
        output_dir
    );
}


void parse_read_feasibility_solution(
    const TerminationReason& termination_reason,
    const VariableMap<double>& result_var_map,
    const ReadVariables& vars,
    TransMap& transmap,
    path output_dir)
{

    path out_path = output_dir/"solution.csv";

    // Check if output path exists already
    bool log_exists = std::filesystem::exists(out_path);

    // Open a file
    ofstream file(out_path, std::ios_base::app);

    if (not file.is_open() or not file.good()){
        throw runtime_error("ERROR: cannot write to file: " + out_path.string());
    }

    // Write header if the file is new
    if (not log_exists){
        file << "sample,read,path" << '\n';
    }

    unordered_set <int64_t> to_be_removed;

    // Print the results of the ILP by iterating all samples, all reads of each sample, and all read/path edges in the transmap
    if (termination_reason == TerminationReason::kOptimal) {
        transmap.for_each_read([&](const string& read_name, int64_t read_id){
            const Variable& var = vars.reads.at(read_id);

            if (var.is_integer()){
                auto is_assigned = bool(int64_t(round(result_var_map.at(var))));

                if (not is_assigned){
                    // Delete all the edges that are not assigned (to simplify iteration later)
                    to_be_removed.emplace(read_id);
                }
            }
            else{
                throw runtime_error("ERROR: solution parsing not implemented for non-integer variables");
            }
        });

        for (const auto& read_id: to_be_removed){
            transmap.remove_node(read_id);
        }
    }
    else{
        cerr << "WARNING: cannot update transmap for non-optimal solution" << '\n';
        transmap = {};
    }
}


void parse_read_feasibility_solution(
    const SolveResult& result,
    const ReadVariables& vars,
    TransMap& transmap,
    path output_dir)
{
    parse_read_feasibility_solution(
        result.termination.reason,
        result.variable_values(),
        vars,
        transmap,
        output_dir
    );
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
        const SolveArguments& args,
        TerminationReason& termination_reason,
    	std::chrono::milliseconds& result_duration,
        unordered_set <pair <int64_t,int64_t> >& result_read_path_edges,
        int64_t n
        ){

    cerr << "Optimizing d for n: " << n << '\n';

    result_read_path_edges.clear();

    auto constraint = model.AddLinearConstraint(vars.cost_n <= n + 0.5);
    model.Minimize(vars.cost_d);

    const absl::StatusOr<SolveResult> response = Solve(model, solver_type, args);

    const auto result = response.value();
    termination_reason = result.termination.reason;

    // Parse Duration
    result_duration = std::chrono::milliseconds(result.solve_stats.solve_time / absl::Milliseconds(1));

    // Immediately undo the constraint after solving
    model.DeleteLinearConstraint(constraint);

    // Check if the first solution is feasible
    if (result.termination.reason != TerminationReason::kOptimal){
        return -1;
    }

    // Extract the d value
    double d = vars.cost_d.Evaluate(result.variable_values());

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

    return d;
}


/**
 * Optimize the assignment of reads to a given number of paths
 * @param model
 * @param vars
 * @param transmap
 * @param solver_type
 * @param result_read_path_edges The set of read-path edges that are assigned in the solution (result)
 * @param d the given distance cost to constrain to
 * @param n_threads
 * @param output_dir
 * @return
 */
double optimize_n_given_d(
        Model& model,
        PathVariables& vars,
        const SolverType& solver_type,
        const SolveArguments& args,
        TerminationReason& termination_reason,
    	std::chrono::milliseconds& result_duration,
        unordered_set <pair <int64_t,int64_t> >& result_read_path_edges,
        double d
        ){

    cerr << "Optimizing n for d: " << d << '\n';

    result_read_path_edges.clear();

    auto constraint = model.AddLinearConstraint(vars.cost_d <= d + 0.5);
    model.Minimize(vars.cost_n);

    const absl::StatusOr<SolveResult> response = Solve(model, solver_type, args);

    const auto result = response.value();
    termination_reason = result.termination.reason;

    // Parse Duration
    result_duration = std::chrono::milliseconds(result.solve_stats.solve_time / absl::Milliseconds(1));

    // Immediately undo the constraint after solving
    model.DeleteLinearConstraint(constraint);

    // Check if the first solution is feasible
    if (result.termination.reason != TerminationReason::kOptimal){
        return -1;
    }

    // Extract the n value
    double n = vars.cost_n.Evaluate(result.variable_values());

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

    return n;
}


/**
 * Optimize the assignment of reads to a given number of paths
 * @param model
 * @param vars
 * @param transmap
 * @param solver_type
 * @param n the given number of paths
 * @param n_threads
 * @return
 */
double optimize_d_given_n(
        Model& model,
        PathVariables& vars,
        const SolverType& solver_type,
        const SolveArguments& args,
        TerminationReason& termination_reason,
    	std::chrono::milliseconds& result_duration,
        int64_t n
        ){

    auto constraint = model.AddLinearConstraint(vars.cost_n <= n + 0.5);
    model.Minimize(vars.cost_d);

    const absl::StatusOr<SolveResult> response = Solve(model, solver_type, args);

    const auto result = response.value();
    termination_reason = result.termination.reason;
    result_duration = std::chrono::milliseconds(result.solve_stats.solve_time / absl::Milliseconds(1));

    // Immediately undo the constraint after solving
    model.DeleteLinearConstraint(constraint);

    // Check if the first solution is feasible
    if (result.termination.reason != TerminationReason::kOptimal){
        return -1;
    }

    // Extract the d value
    double d = vars.cost_d.Evaluate(result.variable_values());

    return d;
}


/**
 * Optimize the assignment of reads to a given number of paths
 * @param model
 * @param vars
 * @param transmap
 * @param solver_type
 * @param termination_reason This is updated so that it is clear why a failure happened
 * @param result_varmap This is updated so that it can be used for the solution in a special cases where n=1
 * @param d the given cost to constrain to
 * @param n_threads
 * @return
 */
double optimize_n_given_d(
        Model& model,
        PathVariables& vars,
        const SolverType& solver_type,
        const SolveArguments& args,
        TerminationReason& termination_reason,
        std::chrono::milliseconds& result_duration,
        int64_t d
        ){

    auto constraint = model.AddLinearConstraint(vars.cost_d <= d + 0.5);
    model.Minimize(vars.cost_n);

    const absl::StatusOr<SolveResult> response = Solve(model, solver_type, args);

    const auto result = response.value();
	termination_reason = result.termination.reason;
    result_duration = std::chrono::milliseconds(result.solve_stats.solve_time / absl::Milliseconds(1));

    // Immediately undo the constraint after solving
    model.DeleteLinearConstraint(constraint);

    // Check if the first solution is feasible
    if (result.termination.reason != TerminationReason::kOptimal){
        return -1;
    }

    // Extract the n value
    double n = vars.cost_n.Evaluate(result.variable_values());

    return n;
}


/**
 * Optimize the assignment of reads to a given number of paths
 * @param model
 * @param vars
 * @param transmap
 * @param solver_type
 * @param d the given cost to constrain to
 * @param n_threads
 * @return
 */
double optimize_n(
        Model& model,
        PathVariables& vars,
        const SolverType& solver_type,
        const SolveArguments& args,
        TerminationReason& termination_reason,
    	std::chrono::milliseconds& result_duration
        ){

    model.Minimize(vars.cost_n);

    const absl::StatusOr<SolveResult> response = Solve(model, solver_type, args);

    const auto result = response.value();
	termination_reason = result.termination.reason;
	result_duration = std::chrono::milliseconds(result.solve_stats.solve_time / absl::Milliseconds(1));

    // Check if the first solution is feasible
    if (result.termination.reason != TerminationReason::kOptimal){
        return -1;
    }

    // Extract the n value
    double n = vars.cost_n.Evaluate(result.variable_values());

    return n;
}


double optimize_d(
        Model& model,
        PathVariables& vars,
        const SolverType& solver_type,
        const SolveArguments& args,
        TerminationReason& termination_reason,
        std::chrono::milliseconds& result_duration
        ){

    model.Minimize(vars.cost_d);

    const absl::StatusOr<SolveResult> response = Solve(model, solver_type, args);

    const auto result = response.value();
	termination_reason = result.termination.reason;
    result_duration = std::chrono::milliseconds(result.solve_stats.solve_time / absl::Milliseconds(1));

    // Check if the first solution is feasible
    if (result.termination.reason != TerminationReason::kOptimal){
        return -1;
    }

    // Extract the d value
    double d = vars.cost_d.Evaluate(result.variable_values());

    return d;
}


double get_normalized_distance(double d, double n, double n_min, double n_max, double d_min, double d_max, double n_weight=1, double d_weight=1){
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

    return sqrt(d_norm*d_norm*d_weight + n_norm*n_norm*n_weight);
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
void optimize_with_golden_search(
        Model& model,
        PathVariables& vars,
        const SolverType& solver_type,
        TransMap& transmap,
        double n_weight,
        double d_weight,
        size_t n_threads,
        path output_dir
){
    TerminationReason tr;
    std::chrono::milliseconds duration;
    VariableMap<double> result_varmap;

    SolveArguments args;
    args.parameters.threads = n_threads;

    unordered_map<int64_t, double> results;
    unordered_map<int64_t, unordered_set <pair <int64_t,int64_t> > > result_edges_cache;

    double n_max = -1;
    double d_min = -1;
    double n_min = -1;
    double d_max = -1;

    Timer t;

    // First find one extreme of the pareto set (D_MIN)
    unordered_set <pair <int64_t,int64_t> > n_max_result_edges;
    d_min = round(optimize_d(model, vars, solver_type, args, tr, duration));

    cerr << t << '\n';
    t.reset();

    n_max = round(optimize_n_given_d(model, vars, solver_type, args, tr, duration, n_max_result_edges, d_min));

    // Need to manually add the result edges to the cache (chicken and egg problem)
    result_edges_cache.emplace(n_max, n_max_result_edges);

    cerr << t << '\n';
    t.reset();

    cerr << "n_max: " << n_max << "\td_min: " << d_min << '\n';

    // Then find the other extreme of the pareto set (N_MIN)
    n_min = round(optimize_n(model, vars, solver_type, args, tr, duration));

    cerr << t << '\n';
    t.reset();

    d_max = round(optimize_d_given_n(model, vars, solver_type, args, tr, duration, result_edges_cache[n_min], n_min));

    cerr << t << '\n';
    t.reset();

    cerr << "n_min: " << n_min << "\td_max: " << d_max << '\n';

    // Within this function the d_cost is now referred to as simply "distance" and the d_cost for a given n value is
    // referred to as "n_distance". 'd' is no loger short for "distance" but is used to refer to the d point which
    // is generated during the sectioning step of the golden section search algorithm, which yields a,b and interior
    // points c,d.
    int64_t n;
    int64_t n_distance = 0;

    // With the range of n determined, we can now use golden section search to find the optimal n and d values
    // Add the extreme values to the results map
    results.emplace(n_min, get_normalized_distance(d_max, n_min, n_min, n_max, d_min, d_max, n_weight, d_weight));
    results.emplace(n_max, get_normalized_distance(d_min, n_max, n_min, n_max, d_min, d_max, n_weight, d_weight));

    // Set arbitrary limit on maximum iterations
    int64_t max_iter = 20;

    int64_t i = 0;
    auto a_i = n_min;
    auto b_i = n_max;

    double phi_inverse = (sqrt(5) - 1)/2;

    // Minimize d for each given n until it can be proven that the d value is optimal (left and right values are larger)
    while (i < max_iter){
        ///        c = b - (b - a) * phi_inv
        ///        d = a + (b - a) * phi_inv
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
            double c = optimize_d_given_n(model, vars, solver_type, args, tr, duration, result_edges_cache[c_i], c_i);
            c_distance = get_normalized_distance(c, c_i, n_min, n_max, d_min, d_max, n_weight, d_weight);
            results.emplace(c_i, c_distance);
            cerr << t << '\n';
            t.reset();
        }
        else{
            c_distance = c_result->second;
        }

        auto d_result = results.find(d_i);
        if (d_result == results.end()){
            double d = optimize_d_given_n(model, vars, solver_type, args, tr, duration, result_edges_cache[d_i], d_i);
            d_distance = get_normalized_distance(d, d_i, n_min, n_max, d_min, d_max, n_weight, d_weight);
            results.emplace(d_i, d_distance);
            cerr << t << '\n';
            t.reset();
        }
        else{
            d_distance = d_result->second;
        }

//        cerr << "a: " << a_i << "\tc: " << c_i << "\td: " << d_i << "\tb: " << b_i << '\n';
//        cerr << "a: " << results.at(a_i) << "\tc: " << c_distance << "\td: " << d_distance << "\tb: " << results.at(b_i) << '\n';

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
    }
    else {
        n = d;
    }

    // Check if the minimum is at the edge of the range (not assessed by the golden section search)
    if (results.at(n_min) < results.at(n)){
        n = n_min;
    }

    // Check if the minimum is at the edge of the range (not assessed by the golden section search)
    if (results.at(n_max) < results.at(n)){
        n = n_max;
    }

    cerr << "Final n: " << n << '\n';

    const auto& optimal_result = result_edges_cache.at(n);

    vector <pair <int64_t,int64_t> > to_be_removed;

    // Filter the transmap using edges of the solution and re-infer the raw distance
    transmap.for_each_sample([&](const string& sample_name, int64_t sample_id){
        transmap.for_each_read_of_sample(sample_id, [&](int64_t read_id){
            transmap.for_each_path_of_read(read_id, [&](int64_t path_id){
                if (optimal_result.find({read_id, path_id}) == optimal_result.end()){
                    to_be_removed.emplace_back(read_id, path_id);
                }
                else{
                    auto [success, w] = transmap.try_get_edge_weight(read_id, path_id);
                    n_distance += round(w);
                }
            });
        });
    });

    for (const auto& [read_id, path_id]: to_be_removed){
        transmap.remove_edge(read_id, path_id);
    }
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
    double n = -1;
    double d = -1;

    bool integral = true;
    if (solver_type == SolverType::kPdlp or solver_type == SolverType::kGlop){
        integral = false;
    }

    construct_joint_n_d_model(transmap, model, vars, integral, use_ploidy_constraint);

    // Solve with golden search
    optimize_with_golden_search(model, vars, solver_type, transmap, n_weight, d_weight, n_threads, output_dir);

    return;
}


void write_optimization_log(
    const TerminationReason& termination_reason,
    const std::chrono::milliseconds& duration,
    const TransMap& transmap,
    const string& name,
    path output_dir
){
	string time_csv;
    duration_to_csv(duration, time_csv);

    size_t n_reads = 0;
    size_t n_paths = 0;
    size_t n_edges = 0;

    transmap.for_each_read([&](const string& name, int64_t id){
    	n_reads++;
    });

    transmap.for_each_path([&](const string& name, int64_t id){
    	n_paths++;
    });

    transmap.for_each_read_to_path_edge([&](int64_t read_id, int64_t path_id, float weight){
    	n_edges++;
    });

    string notes;
    notes += termination_reason_to_string(termination_reason);
    notes += ";";
    notes += "n_reads=" + to_string(n_reads);
    notes += ";";
    notes += "n_paths=" + to_string(n_paths);
    notes += ";";
    notes += "n_edges=" + to_string(n_edges);

    bool success = termination_reason == TerminationReason::kOptimal;

	// Append log line to output file which contains the result of each optimization
    write_time_log(output_dir, name, time_csv, success, notes);

}


TerminationReason optimize_reads_with_d_and_n(
        TransMap& transmap,
        double d_weight,
        double n_weight,
        size_t n_threads,
        size_t time_limit_seconds,
        path output_dir,
        const SolverType& solver_type,
        bool use_ploidy_constraint
        ){

    TerminationReason termination_reason;
    std::chrono::milliseconds duration;
    Model model;
    SolveArguments args;
    PathVariables vars;

    args.parameters.threads = n_threads;

    if (time_limit_seconds > 0){
        args.parameters.time_limit = absl::Seconds(time_limit_seconds);
    }
    else{
        args.parameters.time_limit = absl::InfiniteDuration();
    }

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

    // First find one extreme of the pareto set (D_MIN)
    d_min = round(optimize_d(model, vars, solver_type, args, termination_reason, duration));
    write_optimization_log(termination_reason, duration, transmap, "optimize_d", output_dir);

    if (termination_reason != TerminationReason::kOptimal){
        return termination_reason;
    }

    n_max = round(optimize_n_given_d(model, vars, solver_type, args, termination_reason, duration, d_min));
    write_optimization_log(termination_reason, duration, transmap, "optimize_n_given_d", output_dir);

    if (termination_reason != TerminationReason::kOptimal){
        return termination_reason;
    }

    // Then find the other extreme of the pareto set (N_MIN)
    n_min = round(optimize_n(model, vars, solver_type, args, termination_reason, duration));
    write_optimization_log(termination_reason, duration, transmap, "optimize_n", output_dir);

    if (termination_reason != TerminationReason::kOptimal){
        return termination_reason;
    }

    d_max = round(optimize_d_given_n(model, vars, solver_type, args, termination_reason, duration, n_min));
    write_optimization_log(termination_reason, duration, transmap, "optimize_d_given_n", output_dir);

    if (termination_reason != TerminationReason::kOptimal){
        return termination_reason;
    }

    // Now that we have the range of n and d values, we can normalize the costs and construct the quadratic objective
    auto n_range = n_max - n_min;
    auto d_range = d_max - d_min;

    // Avoid division by zero
    n_range = (n_range == 0 ? 1 : n_range);
    d_range = (d_range == 0 ? 1 : d_range);

    Variable d_norm = model.AddContinuousVariable(0,2,"d");
    Variable n_norm = model.AddContinuousVariable(0,2,"n");

    // Normalize the costs and add 1 to ensure that squaring the normalized values does not make them smaller in the minimization
    model.AddLinearConstraint(d_norm == 1 + (vars.cost_d - d_min)/d_range);
    model.AddLinearConstraint(n_norm == 1 + (vars.cost_n - n_min)/n_range);

    model.Minimize(d_norm*d_norm*d_weight + n_norm*n_norm*n_weight);

    const absl::StatusOr<SolveResult> response_n_d = Solve(model, solver_type, args);
    const auto& result_n_d = response_n_d.value();
    termination_reason = result_n_d.termination.reason;

    // Write a log
    duration = std::chrono::milliseconds(result_n_d.solve_stats.solve_time / absl::Milliseconds(1));
	write_optimization_log(termination_reason, duration, transmap, "optimize_n_d_quadratic", output_dir);

    if (termination_reason != TerminationReason::kOptimal){
        return termination_reason;
    }

    // Infer the optimal n and d values of the joint solution
    n = vars.cost_n.Evaluate(result_n_d.variable_values());
    d = vars.cost_d.Evaluate(result_n_d.variable_values());

    if (integral) {
        parse_read_model_solution(result_n_d, vars, transmap, output_dir);
    }
    else{
        throw runtime_error("ERROR: solution parsing not implemented for non-integer variables");
    }

    return result_n_d.termination.reason;
}


/**
 * This reproduces the function from hapslap
 */
TerminationReason optimize_reads_with_d_plus_n(
        TransMap& transmap,
        double d_weight,
        double n_weight,
        size_t n_threads,
        size_t time_limit_seconds,
        path output_dir,
        const SolverType& solver_type,
        bool use_ploidy_constraint
        ){

    Model model;
    SolveArguments args;
    PathVariables vars;
    TerminationReason termination_reason;
    std::chrono::milliseconds duration;

    args.parameters.threads = n_threads;

    if (time_limit_seconds > 0){
        args.parameters.time_limit = absl::Seconds(time_limit_seconds);
    }
    else{
        args.parameters.time_limit = absl::InfiniteDuration();
    }

    double n_max = -1;
    double d_min = -1;
    double n = -1;
    double d = -1;

    bool integral = true;
    if (solver_type == SolverType::kPdlp or solver_type == SolverType::kGlop){
        integral = false;
    }

    construct_joint_n_d_model(transmap, model, vars, integral, use_ploidy_constraint);

    // First find one extreme of the pareto set (D_MIN)
    d_min = round(optimize_d(model, vars, solver_type, args, termination_reason, duration));
    write_optimization_log(termination_reason, duration, transmap, "optimize_d", output_dir);

    if (termination_reason != TerminationReason::kOptimal){
		return termination_reason;
    }

    n_max = round(optimize_n_given_d(model, vars, solver_type, args, termination_reason, duration, d_min));
    write_optimization_log(termination_reason, duration, transmap, "optimize_n_given_d", output_dir);

    if (termination_reason != TerminationReason::kOptimal){
		return termination_reason;
    }

    // Playing it safe with the variable domains. We actually don't know how much worse the d_max value could be, so
    // using an arbitrary factor of 32.
    Variable d_norm = model.AddContinuousVariable(0,32,"d");
    Variable n_norm = model.AddContinuousVariable(0,n_max,"n");

    // Normalize the costs
    model.AddLinearConstraint(d_norm == vars.cost_d/d_min);
    model.AddLinearConstraint(n_norm == vars.cost_n/n_max);

    model.Minimize(d_norm*d_weight + n_norm*n_weight);

    const absl::StatusOr<SolveResult> response_n_d = Solve(model, solver_type, args);
    const auto& result_n_d = response_n_d.value();

    // Write a log
    duration = std::chrono::milliseconds(result_n_d.solve_stats.solve_time / absl::Milliseconds(1));
    termination_reason = result_n_d.termination.reason;

    write_optimization_log(termination_reason, duration, transmap, "optimize_d_plus_n", output_dir);

    // Check if the final solution is optimal
    if (termination_reason != TerminationReason::kOptimal){
        return termination_reason;
    }

    // Infer the optimal n and d values of the joint solution
    n = vars.cost_n.Evaluate(result_n_d.variable_values());
    d = vars.cost_d.Evaluate(result_n_d.variable_values());

    if (integral) {
        parse_read_model_solution(result_n_d, vars, transmap, output_dir);
    }
    else{
        throw runtime_error("ERROR: solution parsing not implemented for non-integer variables");
    }

    return termination_reason;
}


/**
 * This reproduces the function from hapslap
 */
TerminationReason optimize_read_feasibility(
        TransMap& transmap,
        size_t n_threads,
        size_t time_limit_seconds,
        path output_dir,
        const SolverType& solver_type
        ){

    Model model;
    SolveArguments args;
    ReadVariables vars;

    int64_t r_in = 0;
    int64_t r_out = 0;

    // Count the reads before optimizing feasibility
    transmap.for_each_read([&](const string& name, int64_t id){
      r_in++;
    });

    args.parameters.threads = n_threads;

    if (time_limit_seconds > 0){
        args.parameters.time_limit = absl::Seconds(time_limit_seconds);
    }
    else{
        args.parameters.time_limit = absl::InfiniteDuration();
    }

    double r = -1;

    bool integral = true;
    if (solver_type == SolverType::kPdlp or solver_type == SolverType::kGlop){
        integral = false;
    }

    construct_r_model(transmap, model, vars, integral);

    model.Maximize(vars.cost_r);

    const absl::StatusOr<SolveResult> response = Solve(model, solver_type, args);
    const auto& result = response.value();

    // Convert the Google time object to a usable value
    auto duration = std::chrono::milliseconds(result.solve_stats.solve_time / absl::Milliseconds(1));
    string time_csv;
    duration_to_csv(duration, time_csv);

    bool success = result.termination.reason == TerminationReason::kOptimal;

    // Check if the final solution is feasible
    if (result.termination.reason != TerminationReason::kOptimal){
        // Append log line to output file which contains the result of each optimization
        // If it failed, then r_out is 0 by default
        string notes = "n_read_hap_vars=" + to_string(vars.read_hap.size()) + ";r_in=" + to_string(r_in) + ";r_out=0"  + ";" + termination_reason_to_string(result.termination.reason);
        write_time_log(output_dir, "feasibility", time_csv, success, notes);

        return result.termination.reason;
    }

    // Infer the optimal r value of the feasibility solution
    r = vars.cost_r.Evaluate(result.variable_values());

    if (integral) {
        parse_read_feasibility_solution(result, vars, transmap, output_dir);
    }
    else{
        throw runtime_error("ERROR: solution parsing not implemented for non-integer variables");
    }

    // Count the remaining reads
    transmap.for_each_read([&](const string& name, int64_t id){
      r_out++;
    });

    // Append log line to output file which contains the result of each optimization
    string notes = termination_reason_to_string(result.termination.reason) + ";n_read_hap_vars=" + to_string(vars.read_hap.size()) + ";r_in=" + to_string(r_in) + ";r_out=" + to_string(r_out);
    write_time_log(output_dir, "feasibility", time_csv, success, notes);

    return result.termination.reason;
}


}
