#include "bdsg/include/bdsg/internal/hash_map.hpp"
#include "path_optimizer_mathopt.hpp"
//#include "gurobi_manager.hpp"
#include "Timer.hpp"

#include <fstream>
#include <thread>
#include <map>

using std::ofstream;
using std::tuple;
using std::map;


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
        throw runtime_error("ERROR: unrecognized TerminationReason");
    }
}


void write_optimization_log(
    const TerminationReason& termination_reason,
    const milliseconds& duration,
    const TransMap& transmap,
    const string& name,
    path output_dir
){
    string time_csv;
    duration_to_csv(duration, time_csv);

    size_t n_reads = transmap.get_read_count();
    size_t n_paths = transmap.get_path_count();
    size_t n_edges = 0;

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


/**
 * Construct a model such that each read must be assigned to exactly one path and each sample must have at most 2 paths.
 * For the sake of modularity, no explicit call to CpModelBuilder::Minimize(vars.cost_n + vars.cost_d) is made here, it
 * must be made after constructing the model.
 * @param transmap - TransMap representing the relationship between samples, paths, and reads
 * @param model - model to be constructed
 * @param vars - container to hold ORTools objects which are filled in and queried later after solving
 */
void construct_joint_n_d_model(const TransMap& transmap, Model& model, PathVariables& vars, bool integral, bool use_ploidy_constraint){
    int64_t n_paths;
    unordered_set<int64_t> mandatory_haps;

    // DEFINE: hap vars
    transmap.for_each_path([&](const string& hap_name, int64_t hap_id){\
        string name = "h" + std::to_string(hap_id);
        auto result = vars.haps.emplace(hap_id, model.AddVariable(0,1,integral,name));

        if (transmap.present_haps.contains(hap_id)) {
            auto& h = result.first->second;
            model.AddLinearConstraint(h == 1);
        }
    });

    transmap.for_each_sample([&](const string& sample_name, int64_t sample_id){
        transmap.for_each_read_of_sample(sample_id, [&](int64_t read_id){
            transmap.for_each_path_of_read(read_id, [&](int64_t hap_id) {
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

                // Transmap compression may have already identified edges that are always 1 in the solution
                if (transmap.present_edges.contains(std::make_pair(read_id,hap_id))) {
                    model.AddLinearConstraint(r_h == 1);
                }

                // Do only once for each unique pair of sample-hap
                if (not vars.sample_hap.contains({sample_id, hap_id})){
                    // DEFINE: sample-hap vars
                    string s_h_name = "s" + std::to_string(sample_id) + "h" + std::to_string(hap_id);
                    auto result2 = vars.sample_hap.emplace(std::make_pair(sample_id, hap_id),model.AddVariable(0,1,integral,s_h_name));
                    auto& s_h = result2.first->second;

                    // DEFINE: ploidy
                    vars.ploidy[sample_id] += s_h;

                    // CONSTRAINT: vsh <= vh (indicator for usage of hap w.r.t. sample-hap)
                    model.AddLinearConstraint(s_h <= vars.haps.at(hap_id));

                    // Transmap compression may have identified edges that are always 1 in the solution.
                    // This constraint is technically redundant in the transitive chain of read-hap --> sample-hap
                    // implications, but we add it anyway
                    if (transmap.present_edges.contains(std::make_pair(read_id,hap_id))) {
                        model.AddLinearConstraint(s_h == 1);
                    }
                }

                // CONSTRAINT: vrh <= vsh (indicator for usage of read-hap, w.r.t. sample-hap)
                model.AddLinearConstraint(vars.read_hap.at({read_id, hap_id}) <= vars.sample_hap.at({sample_id, hap_id}));
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
void construct_joint_n_d_model_with_sum_constraints(const TransMap& transmap, Model& model, PathVariables& vars, bool integral, bool use_ploidy_constraint){
    // DEFINE: hap vars
    transmap.for_each_path([&](const string& hap_name, int64_t hap_id){\
        string name = "h" + std::to_string(hap_id);
        auto result = vars.haps.emplace(hap_id, model.AddVariable(0,1,integral,name));

        if (transmap.present_haps.contains(hap_id)) {
            auto& h = result.first->second;
            model.AddLinearConstraint(h == 1);
        }
    });

    unordered_map <pair <int64_t,int64_t>, LinearExpression> sample_hap_read_sums;
    unordered_map <int64_t, LinearExpression> hap_vsh_sums;

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

                // Transmap compression may have already identified edges that are always 1 in the solution
                if (transmap.present_edges.contains(std::make_pair(read_id,hap_id))) {
                    model.AddLinearConstraint(r_h == 1);
                }

                // Do only once for each unique pair of sample-hap
                if (not vars.sample_hap.contains({sample_id, hap_id})){
                    // DEFINE: sample-hap vars
                    string s_h_name = "s" + std::to_string(sample_id) + "h" + std::to_string(hap_id);
                    auto result2 = vars.sample_hap.emplace(std::make_pair(sample_id, hap_id),model.AddVariable(0,1,integral,s_h_name));
                    auto& s_h = result2.first->second;

                    // DEFINE: ploidy
                    vars.ploidy[sample_id] += s_h;

                    // Accumulated sum for use later: sum(vsh) <= vh (indicator for usage of hap w.r.t. sample-hap)
                    hap_vsh_sums[hap_id] += s_h;

                    // Transmap compression may have identified edges that are always 1 in the solution.
                    // This constraint is technically redundant in the transitive chain of read-hap --> sample-hap
                    // implications, but we add it anyway
                    if (transmap.present_edges.contains(std::make_pair(read_id,hap_id))) {
                        model.AddLinearConstraint(s_h == 1);
                    }
                }

                sample_hap_read_sums[{sample_id, hap_id}] += r_h;
            });
        });
    });

    for (auto [sample_hap, read_sum]: sample_hap_read_sums) {
        auto& s_h = vars.sample_hap.at(sample_hap);

        // CONSTRAINT: sum(vrh) <= vsh (indicator for usage of read-hap, w.r.t. sample-hap)
        // model.AddIndicatorConstraint(s_h, read_sum >= 1);
        // model.AddIndicatorConstraint(s_h, read_sum <= 0, true);

        model.AddLinearConstraint(read_sum >= s_h);
        model.AddLinearConstraint(read_sum <= read_sum.terms().size()*s_h);
    }

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

    for (const auto& [hap_id,vsh_sum]: hap_vsh_sums){
        auto& h = vars.haps.at(hap_id);

        // CONSTRAINT: sum(vsh) <= vh (indicator for usage of hap w.r.t. sample-hap)
        // model.AddIndicatorConstraint(h, vsh_sum >= 1);
        // model.AddIndicatorConstraint(h, vsh_sum <= 0, true);

        model.AddLinearConstraint(vsh_sum >= h);
        model.AddLinearConstraint(vsh_sum <= vsh_sum.terms().size()*h);

        // OBJECTIVE: accumulate n cost sum
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
    path output_dir
    )
{
    ofstream file;

    if (not output_dir.empty()){
        // Open a file
        path out_path = output_dir/"solution.csv";
        file.open(out_path);

        if (not file.is_open() or not file.good()){
            throw runtime_error("ERROR: cannot write to file: " + out_path.string());
        }

        // Write header
        file << "sample,read,path" << '\n';
    }

    unordered_set <pair <int64_t, int64_t> > to_be_removed;

    // Print the results of the ILP by iterating all samples, all reads of each sample, and all read/path edges in the transmap
    if (termination_reason == TerminationReason::kOptimal) {
        transmap.for_each_sample([&](const string& sample_name, int64_t sample_id){
            transmap.for_each_read_of_sample(sample_id, [&](int64_t read_id){
                transmap.for_each_path_of_read(read_id, [&](int64_t path_id){
                    const Variable& var = vars.read_hap.at({read_id, path_id});

                    if (var.is_integer()){
                        auto is_assigned = bool(int64_t(round(result_var_map.at(var))));

                        if (is_assigned and not output_dir.empty()){
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


TerminationReason optimize_d(
        const TransMap& transmap,
        const OptimizerConfig& config,
        path output_dir,
        size_t n_threads,
        double& d_min,
        double& n_max
        ){
    TerminationReason termination_reason;
    SolveArguments args;
    Model model;
    PathVariables vars;

    Timer t;

    // Construct the model
    if (config.use_sum_constraints) {
        construct_joint_n_d_model_with_sum_constraints(transmap, model, vars, true, config.use_ploidy_constraint);
    }
    else {
        construct_joint_n_d_model(transmap, model, vars, true, config.use_ploidy_constraint);
    }
    write_time_log(output_dir, "optimize_d_plus_n_construct", t, true);
    t.reset();

    args.parameters.threads = n_threads;

    if (config.timeout_sec > 0){
        args.parameters.time_limit = absl::Seconds(config.timeout_sec);
    }
    else{
        args.parameters.time_limit = absl::InfiniteDuration();
    }

    model.Minimize(vars.cost_d);

    const absl::StatusOr<SolveResult> response = Solve(model, config.solver_type, args);

    const auto result = response.value();
    termination_reason = result.termination.reason;
    auto duration = std::chrono::milliseconds(result.solve_stats.solve_time / absl::Milliseconds(1));

    auto initialize_time = t.elapsed_milliseconds() - duration;
    write_time_log(output_dir, "optimize_d_init", initialize_time, true);
    write_optimization_log(termination_reason, duration, transmap, "optimize_d", output_dir);

    // Check if the first solution is feasible/optimal
    if (termination_reason != TerminationReason::kOptimal){
        return termination_reason;
    }

    // Ideally we would minimize( n | d=d_min ) but we are choosing to be lazy here because it is a costly step
    d_min = round(vars.cost_d.Evaluate(result.variable_values()));
    n_max = round(vars.cost_n.Evaluate(result.variable_values()));

    return termination_reason;
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


TerminationReason optimize_reads_with_d_and_n(
        TransMap& transmap,
        double d_weight,
        double n_weight,
        size_t n_threads,
        size_t time_limit_seconds,
        path output_dir,
        const SolverType& solver_type,
        bool use_ploidy_constraint,
        bool write_solution
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
        if (not write_solution) {
            output_dir.clear();
        }

        parse_read_model_solution(result_n_d, vars, transmap, output_dir);
    }
    else{
        throw runtime_error("ERROR: solution parsing not implemented for non-integer variables");
    }

    return result_n_d.termination.reason;
}


TerminationReason prune_paths_with_d_min(
        TransMap& transmap,
        size_t n_threads,
        path output_dir,
        const OptimizerConfig& config,
        double& d_min_result,
        double& n_max_result
        ){
    Timer t;

    Model model;
    SolveArguments args;
    PathVariables vars;
    TerminationReason termination_reason;
    std::chrono::milliseconds duration;

    args.parameters.threads = n_threads;

    if (config.timeout_sec > 0){
        args.parameters.time_limit = absl::Seconds(config.timeout_sec);
    }
    else{
        args.parameters.time_limit = absl::InfiniteDuration();
    }

    bool integral = true;
    if (config.solver_type == SolverType::kPdlp or config.solver_type == SolverType::kGlop){
        integral = false;
    }

    // Construct the model using the Transmap
    construct_joint_n_d_model(transmap, model, vars, integral, config.use_ploidy_constraint);
    write_time_log(output_dir, "optimize_d_prune_construct", t, true);
    t.reset();

    // First find one extreme of the pareto set (D_MIN)
    model.Minimize(vars.cost_d);

    const absl::StatusOr<SolveResult> response = Solve(model, config.solver_type, args);

    const auto result = response.value();
    termination_reason = result.termination.reason;

    // Record duration reported by MathOpt and the total time outside of that duration
    duration = std::chrono::milliseconds(result.solve_stats.solve_time / absl::Milliseconds(1));
    auto initialize_time = t.elapsed_milliseconds() - duration;

    write_time_log(output_dir, "optimize_d_prune_init", initialize_time, true);
    write_optimization_log(termination_reason, duration, transmap, "optimize_d_prune", output_dir);
    t.reset();

    // EXIT EARLY
    if (termination_reason != TerminationReason::kOptimal){
        return termination_reason;
    }

    if (not integral) {
        throw runtime_error("ERROR: solution parsing not implemented for non-integer variables");
    }

    // Assign the resulting n_max and d_min so they can be reused later
    d_min_result = round(vars.cost_d.Evaluate(result.variable_values()));
    n_max_result = round(vars.cost_n.Evaluate(result.variable_values()));

    vector<int64_t> to_be_removed;

    transmap.sort_adjacency_lists();
    transmap.update_first_of_type();

    transmap.for_each_path([&](const string& path_name, int64_t path_id){
        size_t n_reads = 0;
        bool is_covered = false;

        transmap.for_each_read_of_path(path_id, [&](int64_t read_id){
            const Variable& var = vars.read_hap.at({read_id, path_id});

            if (var.is_integer()){
                auto is_assigned = bool(int64_t(round(result.variable_values().at(var))));

                if (is_assigned){
                    is_covered = true;
                }
            }
            else{
                throw runtime_error("ERROR: solution parsing not implemented for non-integer variables");
            }

            n_reads++;
        });

        // If this path is not an island, and no assigned read-hap edge connects to it, it should be removed
        if (n_reads > 0 and not is_covered) {
            to_be_removed.push_back(path_id);
        }
    });

    for (auto path_id: to_be_removed){
        transmap.remove_node(path_id);
    }

    write_time_log(output_dir, "optimize_d_prune_parse", t, true);

    return termination_reason;

}


/**
 * Perform joint optimization on edit distance + total haps used
 * @param transmap Filled in Transmap which may or may not be compressed
 * @param model Model resulting from calling construct_*_model
 * @param vars Vars corresponding to construct_*_model
 * @param n_threads number of threads to use
 * @param output_dir where to write/append logs
 * @param config configuration for optimization, such as use_ploidy, use_sum_constraints, etc
 * @param d_min used for normalizing the objective terms
 * @param n_max used for normalizing the objective terms
 * @param write_solution if true, write a CSV containing the assigned variables
 * @return result of attempting optimization with ORTools
 */
TerminationReason optimize_reads_with_d_plus_n(
        TransMap& transmap,
        size_t n_threads,
        path output_dir,
        const OptimizerConfig& config,
        double d_min,
        double n_max,
        bool write_solution
        ){

    TerminationReason termination_reason;
    SolveArguments args;
    Model model;
    PathVariables vars;

    Timer t;
    milliseconds duration;

    // Construct the model
    if (config.use_sum_constraints) {
        construct_joint_n_d_model_with_sum_constraints(transmap, model, vars, true, config.use_ploidy_constraint);
    }
    else {
        construct_joint_n_d_model(transmap, model, vars, true, config.use_ploidy_constraint);
    }
    write_time_log(output_dir, "optimize_d_plus_n_construct", t, true);
    t.reset();

    bool integral = true;

    if (config.solver_type == SolverType::kPdlp or config.solver_type == SolverType::kGlop){
        integral = false;
    }

    args.parameters.threads = n_threads;

    if (config.timeout_sec > 0){
        args.parameters.time_limit = absl::Seconds(config.timeout_sec);
    }
    else{
        args.parameters.time_limit = absl::InfiniteDuration();
    }

    // ------------------- Normalize -----------------------------

    // In rare cases, all the edges in the graph are pruned, which indicates that none of the candidates are viable,
    // and therefore the d_min and n_max are 0, resulting in a NaN for the norm step. Here we simply set them to 1
    // so that the solver exits normally and the solution is parsed as given: no read-hap edges.
    //
    // An example of a case where this may happen is in a window where none of the reads span. Generally this is the
    // result of a Graphaligner issue. In some cases a perfect duplication can "capture" all the reads and prevent
    // them from spanning the flanks.
    if (d_min == 0){
        d_min = 1;
    }
    if (n_max == 0){
        n_max = 1;
    }

    // Normalize the costs
    const double d_multiplier = config.d_weight/d_min;
    const double n_multiplier = config.n_weight/n_max;

    // --------------- Minimize joint model -------------------------

    model.Minimize(vars.cost_d*d_multiplier + vars.cost_n*n_multiplier);
    t.reset();

    const absl::StatusOr<SolveResult> response_n_d = Solve(model, config.solver_type, args);
    const auto& result_n_d = response_n_d.value();

    // Write a log
    duration = std::chrono::milliseconds(result_n_d.solve_stats.solve_time / absl::Milliseconds(1));
    termination_reason = result_n_d.termination.reason;

    auto initialize_time = t.elapsed_milliseconds() - duration;
    write_time_log(output_dir, "optimize_d_plus_n_init", initialize_time, true);
    write_optimization_log(termination_reason, duration, transmap, "optimize_d_plus_n", output_dir);
    t.reset();

    // Check if the final solution is optimal
    if (termination_reason != TerminationReason::kOptimal){
        return termination_reason;
    }

    // Infer the optimal n and d values of the joint solution
    double n = vars.cost_n.Evaluate(result_n_d.variable_values());
    double d = vars.cost_d.Evaluate(result_n_d.variable_values());

    if (integral) {
        if (not write_solution) {
            parse_read_model_solution(result_n_d, vars, transmap, "");
        }
        else {
            parse_read_model_solution(result_n_d, vars, transmap, output_dir);
        }
    }
    else{
        throw runtime_error("ERROR: solution parsing not implemented for non-integer variables");
    }

    write_time_log(output_dir, "optimize_d_plus_n_parse", t, true);


    cerr << "-----------------" << '\n' <<
            output_dir << '\n' <<
            "Objective=" << to_string(result_n_d.objective_value()) << '\n' <<
            "d_min=" << to_string(d_min) << '\n' <<
            "n_max=" << to_string(n_max) << '\n' <<
            "n=" << to_string(n) << '\n' <<
            "d=" << to_string(d) << '\n' <<
            "-----------------" << '\n';

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
    Timer t;

    Model model;
    SolveArguments args;
    ReadVariables vars;

    int64_t r_in = int64_t(transmap.get_read_count());
    int64_t r_out = 0;

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
    write_time_log(output_dir, "feasibility_construct", t, true);
    t.reset();

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

    auto initialize_time = t.elapsed_milliseconds() - duration;
    write_time_log(output_dir, "feasibility_init", initialize_time, true);
    t.reset();

    // Infer the optimal r value of the feasibility solution
    r = vars.cost_r.Evaluate(result.variable_values());

    if (integral) {
        parse_read_feasibility_solution(result, vars, transmap, output_dir);
    }
    else{
        throw runtime_error("ERROR: solution parsing not implemented for non-integer variables");
    }

    // Count the remaining reads
    r_out = int64_t(transmap.get_read_count());

    // Append log line to output file which contains the result of each optimization
    string notes = termination_reason_to_string(result.termination.reason) + ";n_read_hap_vars=" + to_string(vars.read_hap.size()) + ";r_in=" + to_string(r_in) + ";r_out=" + to_string(r_out);
    write_time_log(output_dir, "feasibility", time_csv, success, notes);

    return result.termination.reason;
}


TerminationReason optimize(TransMap& transmap, const OptimizerConfig& config, path subdir, bool write_solution){
    TerminationReason termination_reason;

    double d_min = -1;
    double n_max = -1;

    // First resolve any samples that break ploidy feasibility by removing the minimum # of reads
    termination_reason = optimize_read_feasibility(transmap, 1, config.timeout_sec, subdir, config.solver_type);

    if (termination_reason != TerminationReason::kOptimal) {
        return termination_reason;
    }

    // Find d_min and n_max and also prune the haps from the graph based on the d_min solution
    if (config.prune_with_d_min) {
        termination_reason = prune_paths_with_d_min(transmap, 1, subdir, config, d_min, n_max);
    }
    else{
        termination_reason = optimize_d(transmap, config, subdir, 1, d_min, n_max);
    }

    if (termination_reason != TerminationReason::kOptimal) {
        return termination_reason;
    }

    termination_reason = optimize_reads_with_d_plus_n(transmap, 1, subdir, config, d_min, n_max, write_solution);

    return termination_reason;
}


TerminationReason optimize_compressed(TransMap& transmap, const OptimizerConfig& config, path subdir, bool write_solution){
    Timer t;

    TerminationReason termination_reason;

    double d_min = -1;
    double n_max = -1;

    // First resolve any samples that break ploidy feasibility by removing the minimum # of reads
    termination_reason = optimize_read_feasibility(transmap, 1, config.timeout_sec, subdir, config.solver_type);

    if (termination_reason != TerminationReason::kOptimal) {
        return termination_reason;
    }

    TransMap clone;

    // TODO: cannot currently use read compression for the read_feasibility step, unless we store # of reads collapsed
    // t.reset();
    // transmap.sort_adjacency_lists();
    // transmap.update_first_of_type();
    // transmap.get_mandatory_haplotypes();
    // transmap.compress_haplotypes_global(0);
    // write_time_log(subdir, "compress_transmap_initial", t, true);
    // t.reset();

    // Branch off from the parent Transmap after common compression steps
    clone = transmap;

    // For d_min step, n_weight is 0 and d_weight is 1 because n is not considered in the objective
    clone.sort_adjacency_lists();
    clone.update_first_of_type();
    clone.get_mandatory_haplotypes();
    clone.compress_haplotypes_global(0);
    clone.compress_haplotypes_local(0,1,0);
    clone.solve_easy_samples(0,1,0);
    clone.compress_reads(0);
    clone.compress_samples(0);

    // Find d_min and n_max
    if (config.prune_with_d_min) {
        // Also prune the haps from the graph based on the d_min solution, if requested
        termination_reason = prune_paths_with_d_min(clone, 1, subdir, config, d_min, n_max);

        // Since we are working from a clone, need to also apply the changes to the parent transmap
        vector<int64_t> to_be_deleted;
        transmap.for_each_path([&](const string& name, int64_t id) {
            if (not clone.has_node(name)) {
                to_be_deleted.emplace_back(id);
            }
        });

        // Must not modify while iterating Transmap
        for (auto id: to_be_deleted) {
            transmap.remove_node(id);
        }
    }
    else{
        termination_reason = optimize_d(clone, config, subdir, 1, d_min, n_max);
    }

    if (termination_reason != TerminationReason::kOptimal) {
        return termination_reason;
    }

    // Branch off from the parent Transmap
    clone = transmap;

    // For the final joint optimization we need to tell the compression methods what constants are used in the objective
    // for both of the n and d terms
    auto c_n = float(config.n_weight/n_max);
    auto c_d = float(config.d_weight/d_min);

    clone.sort_adjacency_lists();
    clone.update_first_of_type();
    clone.get_mandatory_haplotypes();
    clone.compress_haplotypes_global(0);
    clone.compress_haplotypes_local(c_n,c_d,0);
    clone.solve_easy_samples(c_n,c_d,0);
    clone.compress_reads(0);
    clone.compress_samples(0);

    termination_reason = optimize_reads_with_d_plus_n(clone, 1, subdir, config, d_min, n_max, write_solution);

    // Finally decompress
    clone.decompress_samples();

    // Overwrite the original transmap
    transmap = clone;

    return termination_reason;
}


TerminationReason optimize_samplewise(TransMap& transmap, const OptimizerConfig& config, path subdir, bool write_solution) {

    vector <pair<int64_t,int64_t> > edges_to_remove;
    TerminationReason termination_reason = TerminationReason::kOptimal;

    transmap.for_each_sample([&](const string& sample_name, int64_t sample_id) {
        TransMap submap;
        transmap.extract_sample_as_transmap(sample_name, submap);

        termination_reason = optimize(transmap, config, subdir, write_solution);

        // cerr << sample_name << ' ' << termination_reason_to_string(termination_reason) << '\n';

        // TODO handle this more smarter.. what about kFeasible?
        if (termination_reason != TerminationReason::kOptimal) {
            return;
        }

        // Find edges in transmap that weren't part of the samplewise solution before continuing to next sample
        transmap.for_each_read_of_sample(sample_name, [&](const string& read_name, int64_t read_id) {
            transmap.for_each_path_of_read(read_id, [&](const string& path_name, int64_t path_id) {
                // Convert names to submap IDs
                auto a = submap.get_id(read_name);
                auto b = submap.get_id(path_name);

                if (not submap.has_edge(a, b)) {
                    edges_to_remove.emplace_back(a,b);
                }
            });
        });
    });

    for (const auto& [a,b]: edges_to_remove) {
        transmap.remove_edge(a, b);
    }

    return termination_reason;
}


/// Rescale edge weights for each read as quadratic function of distance from best weight
void rescale_weights_as_quadratic_best_diff(TransMap& transmap, float domain_min, float domain_max){
    unordered_map<int64_t,float> best_weights;

    // First find the best weight for each read
    transmap.for_each_read_to_path_edge([&](int64_t read_id, int64_t path_id, float weight) {
        auto result = best_weights.find(read_id);

        if (result == best_weights.end()) {
            best_weights[read_id] = weight;
        }
        else {
            best_weights[read_id] = min(result->second, weight);
        }
    });

    vector<tuple<int64_t,int64_t,float> > edges_to_add;

    // Rescale the weights
    transmap.for_each_read_to_path_edge([&](int64_t read_id, int64_t path_id, float weight) {
        float best = best_weights[read_id];
        float w = round(pow((weight - best), 1.2) + 1);

        if (w < domain_min) {
            throw runtime_error("ERROR: weight rescaling resulted in weight below minimum: " + to_string(w));
        }

        // if exceeds max, then just clip it (for the sanity of building the model without int overflow)
        if (w > domain_max) {
            w = domain_max;
        }

//        cerr << "rescaling: " << read_id << ',' << path_id << ' ' << weight << ' ' << best << ' ' << w << '\n';

        edges_to_add.emplace_back(read_id, path_id, w);
    });

    // Update the DS
    for (const auto& [read_id, path_id, weight]: edges_to_add){
        // Will overwrite the edge if it already exists
        transmap.add_edge(read_id, path_id, weight);
    }
}





}
