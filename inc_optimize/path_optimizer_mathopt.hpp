#pragma once
#include "TransitiveMap.hpp"
#include "bdsg/include/bdsg/internal/hash_map.hpp"

#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <string>
#include <vector>
#include <set>

using std::unordered_map;
using std::unordered_set;
using std::vector;
using std::string;
using std::pair;
using std::set;

#include "ortools/base/logging.h"
#include "ortools/math_opt/cpp/math_opt.h"
#include "ortools/util/sorted_interval_list.h"

using operations_research::math_opt::Model;
using operations_research::math_opt::Variable;
using operations_research::math_opt::LinearExpression;
using operations_research::math_opt::BoundedLinearExpression;
using operations_research::math_opt::QuadraticExpression;
using operations_research::math_opt::Sum;
using operations_research::math_opt::Solve;
using operations_research::math_opt::SolveResult;
using operations_research::math_opt::SolverType;
using operations_research::math_opt::TerminationReason;
using operations_research::math_opt::SolveArguments;
using operations_research::math_opt::Termination;
using operations_research::math_opt::VariableMap;


namespace sv_merge {


class OptimizerConfig {
public:
        OptimizerConfig() = default;

        size_t timeout_sec = 30*60;
        float min_read_hap_identity = 0.5;
        float d_weight = 1.0;
        float n_weight = 1.0;
        bool skip_solve = false;
        bool rescale_weights = false;
        bool use_quadratic_objective = false;
        bool use_golden_search = false;
        bool use_sum_constraints = true;
        bool samplewise = false;
        bool prune_with_d_min = false;
        bool use_ploidy_constraint = true;
        size_t max_reads_per_sample = 0;
        SolverType solver_type = SolverType::kGscip;
};


string termination_reason_to_string(const TerminationReason& reason);


/// Data class to contain the variables of the ILP, int64_t IDs intended to be reused from the TransMap
class PathVariables{
public:
    unordered_map <pair<int64_t,int64_t>, Variable> sample_hap;
    unordered_map <pair<int64_t,int64_t>, Variable> read_hap;

    unordered_map <int64_t, LinearExpression> read_flow;
    unordered_map <int64_t, LinearExpression> ploidy;
    unordered_map <int64_t, Variable> haps;

    LinearExpression cost_d;
    LinearExpression cost_n;
};


/// Data class to contain the variables of the ILP, int64_t IDs intended to be reused from the TransMap
class ReadVariables{
public:
    unordered_map <pair<int64_t,int64_t>, Variable> sample_hap;
    unordered_map <pair<int64_t,int64_t>, Variable> read_hap;

    unordered_map <int64_t, LinearExpression> read_flow;
    unordered_map <int64_t, LinearExpression> ploidy;
    unordered_map <int64_t, Variable> haps;
    unordered_map <int64_t, Variable> reads;

    LinearExpression cost_r;
};


/**
 * Construct a model such that each read must be assigned to exactly one path and each sample must have at most 2 paths.
 * For the sake of modularity, no explicit call to CpModelBuilder::Minimize(vars.cost_d) is made here, it must be made
 * after constructing the model.
 * @param transmap - TransMap representing the relationship between samples, paths, and reads
 * @param model - model to be constructed
 * @param vars - CPSAT BoolVar and LinearExpr objects to be filled in
 */
void construct_d_model(const TransMap& transmap, Model& model, PathVariables& vars);

/**
 * Construct a model such that each read must be assigned to exactly one path and each sample must have at most 2 paths.
 * For the sake of modularity, no explicit call to CpModelBuilder::Minimize(vars.cost_n + vars.cost_d) is made here, it
 * must be made after constructing the model.
 * @param transmap - TransMap representing the relationship between samples, paths, and reads
 * @param model - model to be constructed
 * @param vars - CPSAT BoolVar and LinearExpr objects to be filled in
 */
void construct_joint_n_d_model(
        const TransMap& transmap,
        Model& model,
        PathVariables& vars,
        bool integral = true,
        bool use_ploidy_constraint = true
        );

TerminationReason optimize_reads_with_d_and_n(
        TransMap& transmap,
        double d_weight,
        double n_weight,
        size_t n_threads,
        size_t time_limit_seconds,
        path output_dir,
        const SolverType& solver_type,
        bool use_ploidy_constraint = true,
        bool write_solution = false
        );

TerminationReason optimize_reads_with_d_plus_n(
        TransMap& transmap,
        size_t n_threads,
        path output_dir,
        const OptimizerConfig& config,
        double d_min = -1,
        double n_max = -1,
        bool write_solution = false
        );

TerminationReason prune_paths_with_d_min(
        TransMap& transmap,
        size_t n_threads,
        path output_dir,
        const OptimizerConfig& config,
        double& d_min_result,
        double& n_max_result
);

void optimize_reads_with_d_and_n_using_golden_search(
        TransMap& transmap,
        double d_weight,
        double n_weight,
        size_t n_threads,
        path output_dir,
        const SolverType& solver_type,
        bool use_ploidy_constraint = true
        );

TerminationReason optimize_read_feasibility(
        TransMap& transmap,
        size_t n_threads,
        size_t time_limit_seconds,
        path output_dir,
        const SolverType& solver_type
);

TerminationReason optimize(TransMap& transmap, const OptimizerConfig& config, path subdir, bool write_solution = false);

TerminationReason optimize_samplewise(TransMap& transmap, const OptimizerConfig& config, path subdir, bool write_solution = false);

void rescale_weights_as_quadratic_best_diff(TransMap& transmap, float domain_min, float domain_max);

}
