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
using operations_research::sat::SatParameters;
using operations_research::sat::LinearExpr;


namespace sv_merge {


/// Data class to contain the variables of the ILP, int64_t IDs intended to be reused from the TransMap
class PathVariables{
public:
    unordered_map <pair<int64_t,int64_t>, BoolVar> path_to_read;
    unordered_map <pair<int64_t,int64_t>, BoolVar> path_to_sample;

    // Will be left empty in the unary objective d_model
    unordered_map <int64_t, BoolVar> path_indicators;

    LinearExpr cost_d;

    // Will be left default initialized in the unary objective d_model
    LinearExpr cost_n;
};

/**
 * Construct a model such that each read must be assigned to exactly one path and each sample must have at most 2 paths.
 * For the sake of modularity, no explicit call to CpModelBuilder::Minimize(vars.cost_d) is made here, it must be made
 * after constructing the model.
 * @param transmap - TransMap representing the relationship between samples, paths, and reads
 * @param model - model to be constructed
 * @param vars - CPSAT BoolVar and LinearExpr objects to be filled in
 */
void construct_d_model(const TransMap& transmap, CpModelBuilder& model, PathVariables& vars);

/**
 * Construct a model such that each read must be assigned to exactly one path and each sample must have at most 2 paths.
 * For the sake of modularity, no explicit call to CpModelBuilder::Minimize(vars.cost_n + vars.cost_d) is made here, it
 * must be made after constructing the model.
 * @param transmap - TransMap representing the relationship between samples, paths, and reads
 * @param model - model to be constructed
 * @param vars - CPSAT BoolVar and LinearExpr objects to be filled in
 */
void construct_joint_n_d_model(const TransMap& transmap, CpModelBuilder& model, PathVariables& vars);

void optimize_d(const TransMap& transmap);

void optimize_d_plus_n(const TransMap& transmap, int64_t d_coeff, int64_t n_coeff);

void optimize_reads_with_d_and_n(TransMap& transmap, int64_t d_weight, int64_t n_weight, size_t n_threads, path output_dir);

}
