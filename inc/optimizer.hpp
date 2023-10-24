#pragma once
#include "TransitiveMap.hpp"
#include "pair_hash.hpp"

#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <string>
#include <vector>

using std::unordered_map;
using std::unordered_set;
using std::vector;
using std::string;
using std::pair;

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


namespace sv_merge {


/// Data class to contain the variables of the ILP, int64_t IDs intended to be reused from the TransMap
class Variables{
public:
    unordered_map <pair<int64_t,int64_t>, BoolVar> path_to_read;

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
void construct_d_model(const TransMap& transmap, CpModelBuilder& model, Variables& vars){
    // TODO: reserve Variables map data structures based on how many edges there are known to be in the transmap?

    // Read->path boolean indicators
    transmap.for_each_read_to_path_edge([&](int64_t read_id, int64_t path_id, float weight){
        pair<int64_t,int64_t> p = {path_id, read_id};

        // Update the boolean var map, and keep a reference to the inserted BoolVar
        const auto result = vars.path_to_read.emplace(p,model.NewBoolVar());
        const auto& v = result.first->second;

        auto w = int64_t(weight);
        vars.cost_d += v*w;
    });

    transmap.for_each_sample([&](const string& name, int64_t sample_id){
        LinearExpr samplewise_vars;

        // Sample->read boolean indicators
        transmap.for_each_read_of_sample(sample_id, [&](int64_t read_id){
            LinearExpr readwise_vars;

            // Read->path boolean indicators
            transmap.for_each_path_of_read(read_id, [&](int64_t path_id){
                const BoolVar& var = vars.path_to_read.at({path_id, read_id});
                samplewise_vars += var;
                readwise_vars += var;
            });

            // Enforce (sum(r) == 1). i.e. read must be assigned to exactly one path
            model.AddEquality(readwise_vars, 1);
        });

        // Enforce (sum(s) <= 2). i.e. sample must be assigned at most 2 paths
        model.AddLessOrEqual(samplewise_vars, 2);
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
void construct_joint_n_d_model(const TransMap& transmap, CpModelBuilder& model, Variables& vars){
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


void optimize_d(TransMap& transmap){
    CpModelBuilder model;
    Variables vars;

    construct_d_model(transmap, model, vars);

    model.Minimize(vars.cost_d);

    const CpSolverResponse response = Solve(model.Build());

    if (response.status() == CpSolverStatus::OPTIMAL || response.status() == CpSolverStatus::FEASIBLE) {

        cerr << "Maximum of objective function: " << response.objective_value() << '\n';

        // Iterate the read assignment variables and print them
        transmap.for_each_read_to_path_edge([&](int64_t path_id, int64_t read_id, float weight){
            cerr << path_id << ',' << read_id << '\n';

            const auto& path_name = transmap.get_node(path_id).name;
            cerr << path_name << '\n' << std::flush;

            const auto& read_name = transmap.get_node(read_id).name;
            cerr << read_name << '\n' << std::flush;

            string sample_name;
            transmap.get_read_sample(read_id, sample_name);
            cerr << sample_name << '\n' << std::flush;

            const auto& var = vars.path_to_read.at({path_id, read_id});
            cerr << var << '\n' << std::flush;

            cerr << sample_name << ',' << path_name << ',' << read_name << " = " << SolutionIntegerValue(response, var) << '\n';
        });
    }
    else {
        cerr << "No solution found." << '\n';
    }

    cerr << "Statistics" << '\n';
    cerr << CpSolverResponseStats(response) << '\n';
}

}
