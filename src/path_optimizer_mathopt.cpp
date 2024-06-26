#include "path_optimizer_mathopt.hpp"
#include "bdsg/include/bdsg/internal/hash_map.hpp"
#include <fstream>

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
void construct_joint_n_d_model(const TransMap& transmap, Model& model, PathVariables& vars){
    // DEFINE: hap vars
    transmap.for_each_path([&](const string& hap_name, int64_t hap_id){\
        string name = "h" + std::to_string(hap_id);
        vars.haps.emplace(hap_id, model.AddBinaryVariable(name));
    });

    transmap.for_each_sample([&](const string& name, int64_t sample_id){
        transmap.for_each_read_of_sample(sample_id, [&](int64_t read_id){
            transmap.for_each_path_of_read(read_id, [&](int64_t hap_id) {
                // DEFINE: read-hap vars
                string r_h_name = "r" + std::to_string(read_id) + "h" + std::to_string(hap_id);
                auto result = vars.read_hap.emplace(std::make_pair(read_id, hap_id), model.AddBinaryVariable(r_h_name));
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
                    auto result2 = vars.sample_hap.emplace(std::make_pair(sample_id, hap_id),model.AddBinaryVariable(s_h_name));
                    auto& s_h = result2.first->second;

                    // DEFINE: ploidy
                    vars.ploidy[hap_id] += s_h;

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

    // CONSTRAINT: ploidy
    for (const auto& [hap_id,p]: vars.ploidy){
        model.AddLinearConstraint(p <= 2);
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
        return;
    }

    // Write header
    file << "sample,read,path" << '\n';

    // Print the results of the ILP by iterating all samples, all reads of each sample, and all read/path edges in the transmap
    if (result_n_d.termination.reason == TerminationReason::kOptimal) {
        transmap.for_each_sample([&](const string& sample_name, int64_t sample_id){
            transmap.for_each_read_of_sample(sample_id, [&](int64_t read_id){
                transmap.for_each_path_of_read(read_id, [&](int64_t path_id){
                    const Variable& var = vars.read_hap.at({read_id, path_id});
                    bool is_assigned = result_n_d.variable_values().at(var);
                    if (is_assigned){
                        file << sample_name << ',' << transmap.get_node(read_id).name << ',' << transmap.get_node(path_id).name << '\n';
                    }
                    else{
                        // Delete all the edges that are not assigned (to simplify iteration later)
                        transmap.remove_edge(read_id, path_id);
                    }
                });
            });
        });
    }
}


void optimize_reads_with_d_and_n(TransMap& transmap, double d_weight, double n_weight, size_t n_threads, path output_dir, const SolverType& solver_type){
    cerr << "solver_type: " << solver_type << '\n';

    Model model;
    PathVariables vars;

    construct_joint_n_d_model(transmap, model, vars);

    // First find one extreme of the pareto set (using a tie-breaker cost for the other objective)
    model.Minimize(vars.cost_d + vars.cost_n*1e-6);

    const absl::StatusOr<SolveResult> response_d = Solve(model, solver_type);
    const auto result_d = response_d.value();

    // Check if the first solution is feasible
    if (result_d.termination.reason != TerminationReason::kOptimal){
        cerr << "WARNING: no solution for d_min found: " << output_dir << '\n';
        transmap = {};
        return;
    }

    // Infer the n and d values of the D_MIN solution (ignoring tie-breaker cost)
    double n_max = vars.cost_n.Evaluate(result_d.variable_values());
    double d_min = vars.cost_d.Evaluate(result_d.variable_values());

    // Then find the other extreme of the pareto set (using a tie-breaker cost for the other objective)
    model.Minimize(vars.cost_n + vars.cost_d*1e-6);

    const absl::StatusOr<SolveResult> response_n = Solve(model, solver_type);
    const auto result_n = response_n.value();

    // Check if the second solution is feasible
    if (result_n.termination.reason != TerminationReason::kOptimal){
        cerr << "WARNING: no solution for n_min found: " << output_dir << '\n';
        transmap = {};
        return;
    }

    // Infer the n and d values of the N_MIN solution (ignoring tie-breaker cost)
    double n_min = vars.cost_n.Evaluate(result_n.variable_values());
    double d_max = vars.cost_d.Evaluate(result_n.variable_values());

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

    cerr << "solving joint objective" << '\n';
    const absl::StatusOr<SolveResult> response_n_d = Solve(model, solver_type);
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

    // Infer the n and d values of the N_MIN solution (ignoring tie-breaker cost)
    double n = vars.cost_n.Evaluate(result_n.variable_values());
    double d = vars.cost_d.Evaluate(result_n.variable_values());

    cerr << "n: " << n << "\td: " << d << '\n';

    parse_read_model_solution(result_n_d, vars, transmap, output_dir);
}



}