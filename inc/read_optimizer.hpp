#pragma once
#include "TransitiveMap.hpp"
#include "pair_hash.hpp"

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
using operations_research::sat::LinearExpr;


namespace sv_merge {


/// Data class to contain the variables of the ILP, int64_t IDs intended to be reused from the TransMap
class ReadVariables{
public:
    // Edges where the left node is the parent
    unordered_map <pair<int64_t,int64_t>, BoolVar> edges;
    unordered_map <int64_t, BoolVar> is_parent;

    // Cost of all edges assigned
    LinearExpr cost_d;

    // Cost of adding parents
    LinearExpr cost_n;
};


void construct_read_model(
        const TransMap& transmap,
        int64_t sample_id,
        CpModelBuilder& model,
        ReadVariables& vars,
        vector<int64_t>& representatives);

void optimize_reads(TransMap& transmap, int64_t sample_id, vector<int64_t>& representatives);


}
