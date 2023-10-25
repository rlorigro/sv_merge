#include <stdint.h>
#include <stdlib.h>

#include "ortools/base/logging.h"
#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model.pb.h"
#include "ortools/sat/cp_model_solver.h"
#include "ortools/util/sorted_interval_list.h"
// [END import]

using operations_research::sat::CpModelBuilder;
using operations_research::Domain;
using operations_research::sat::IntVar;
using operations_research::sat::CpSolverResponse;
using operations_research::sat::CpSolverStatus;


/// Taken from: https://github.com/google/or-tools/blob/stable/ortools/sat/samples/cp_sat_example.cc
void CpSatExample() {
    // [START model]
    CpModelBuilder cp_model;
    // [END model]

    // [START variables]
    int64_t var_upper_bound = std::max({50, 45, 37});
    const Domain domain(0, var_upper_bound);
    const IntVar x = cp_model.NewIntVar(domain).WithName("x");
    const IntVar y = cp_model.NewIntVar(domain).WithName("y");
    const IntVar z = cp_model.NewIntVar(domain).WithName("z");
    // [END variables]

    // [START constraints]
    cp_model.AddLessOrEqual(2 * x + 7 * y + 3 * z, 50);
    cp_model.AddLessOrEqual(3 * x - 5 * y + 7 * z, 45);
    cp_model.AddLessOrEqual(5 * x + 2 * y - 6 * z, 37);
    // [END constraints]

    // [START objective]
    cp_model.Maximize(2 * x + 2 * y + 3 * z);
    // [END objective]

    // Solving part.
    // [START solve]
    const CpSolverResponse response = Solve(cp_model.Build());
    // [END solve]

    // [START print_solution]
    if (response.status() == CpSolverStatus::OPTIMAL ||
        response.status() == CpSolverStatus::FEASIBLE) {
        // Get the value of x in the solution.
        LOG(INFO) << "Maximum of objective function: "
                  << response.objective_value();
        LOG(INFO) << "x = " << SolutionIntegerValue(response, x);
        LOG(INFO) << "y = " << SolutionIntegerValue(response, y);
        LOG(INFO) << "z = " << SolutionIntegerValue(response, z);
    } else {
        LOG(INFO) << "No solution found.";
    }
    // [END print_solution]

    // Statistics.
    // [START statistics]
    LOG(INFO) << "Statistics";
    LOG(INFO) << CpSolverResponseStats(response);
    // [END statistics]
}

int main() {
    CpSatExample();
    return EXIT_SUCCESS;
}