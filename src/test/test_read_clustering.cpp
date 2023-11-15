#include "TransitiveMap.hpp"
#include "pair_hash.hpp"

using sv_merge::TransMap;

#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <string>
#include <vector>
#include <set>
#include <iostream>
#include <stdexcept>

using std::runtime_error;
using std::cerr;
using std::unordered_map;
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



/// Data class to contain the variables of the ILP, int64_t IDs intended to be reused from the TransMap
class Variables{
public:
    // Edges where the left node is the parent
    unordered_map <pair<int64_t,int64_t>, BoolVar> edges;
    unordered_map <int64_t, BoolVar> is_parent;
    unordered_map <int64_t, BoolVar> is_child;

    LinearExpr cost_d;

    // Will be left default initialized in the unary objective d_model
    LinearExpr cost_n;
};


int main(){
    TransMap transmap;

    transmap.add_sample("HG001");

    transmap.add_read("read_01");
    transmap.add_read("read_02");
    transmap.add_read("read_03");
    transmap.add_read("read_04");

    transmap.add_edge("read_01", "HG001");
    transmap.add_edge("read_02", "HG001");
    transmap.add_edge("read_03", "HG001");
    transmap.add_edge("read_04", "HG001");

    // All vs All read-read edges
    transmap.add_edge("read_01", "read_02", 2);
    transmap.add_edge("read_01", "read_03", 2);
    transmap.add_edge("read_01", "read_04", 1);
    transmap.add_edge("read_02", "read_03", 1);
    transmap.add_edge("read_02", "read_04", 2);
    transmap.add_edge("read_03", "read_04", 2);

    cerr << '\n';

    transmap.for_each_sample([&](const string& sample_name, int64_t sample_id){
        CpModelBuilder model;
        Variables vars;

        transmap.for_each_read_of_sample(sample_id, [&](int64_t id) {
            auto& p = vars.is_parent.emplace(id, model.NewBoolVar()).first->second;
            auto& c = vars.is_child.emplace(id, model.NewBoolVar()).first->second;

            // Any node can only be a parent or a child, one implies NOT the other
            // TODO: Does this need its reciprocal??
            model.AddEquality(p,Not(c));

            vars.cost_n += p;
        });

        // TODO use BFS/connected components to enumerate edges instead, some edges won't be present
        transmap.for_each_read_of_sample(sample_id, [&](int64_t id_a){
            transmap.for_each_read_of_sample(sample_id, [&](int64_t id_b){
                auto [success,w] = transmap.try_get_edge_weight(id_a, id_b);

                if (not success){
                    return;
                }

                pair<int64_t,int64_t> p = {id_a, id_b};
                auto& result = vars.edges.emplace(p, model.NewBoolVar()).first->second;

                // Choosing this edge implies that the left is parent and the right is child
                auto& p_a = vars.is_parent.at(id_a);
                auto& c_b = vars.is_child.at(id_b);

                // TODO: FIX THIS, need to change to at least one, not a strict equality
                model.AddEquality(result, p_a);
                model.AddEquality(result, c_b);

                vars.cost_d += w*p_a;
            });
        });

        model.AddGreaterThan(vars.cost_n, 0);
        model.Minimize(vars.cost_d);


        const CpSolverResponse response = Solve(model.Build());

        if (response.status() == CpSolverStatus::OPTIMAL || response.status() == CpSolverStatus::FEASIBLE) {
            cerr << "Maximum of objective function: " << response.objective_value() << '\n';

            set<string> results;

            // Iterate the parent assignment variables and print them
            for (const auto& [id,var]: vars.is_parent){
                auto name = transmap.get_node(id).name;

                bool is_assigned = (SolutionIntegerValue(response, var) == 1);
                if (is_assigned){
                    results.insert(name + ": parent");
                }
            }

            // Iterate the child assignment variables and print them
            for (const auto& [id,var]: vars.is_child){
                auto name = transmap.get_node(id).name;

                bool is_assigned = (SolutionIntegerValue(response, var) == 1);
                if (is_assigned){
                    results.insert(name + ": child");
                }
            }

            // Iterate the child assignment variables and print them
            for (const auto& [e,var]: vars.edges){
                auto [a,b] = e;
                auto name_a = transmap.get_node(a).name;
                auto name_b = transmap.get_node(b).name;

                bool is_assigned = (SolutionIntegerValue(response, var) == 1);
                cerr << is_assigned << ' ' << name_a << "->" << name_b << '\n';
            }

            for (auto& item: results){
                cerr << item << '\n';
            }

        }
        else {
            cerr << "No solution found." << '\n';
        }

        cerr << '\n';
        cerr << "Statistics" << '\n';
        cerr << CpSolverResponseStats(response) << '\n';

    });


    return 0;
}
