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

    // Cost of all edges assigned
    LinearExpr cost_d;

    // Cost of adding parents
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
            // Construct indicator booleans which represent parent status of each node
            auto& p = vars.is_parent.emplace(id, model.NewBoolVar()).first->second;

            // Keep track of the total number of parents assigned
            // TODO: currently unused
            vars.cost_n += p;
        });

        unordered_map<int64_t, LinearExpr> child_edge_sums;
        unordered_map<int64_t, LinearExpr> parent_edge_sums;

        int64_t n_reads = 0;
        int64_t n_edges = 0;

        // TODO use BFS/connected components to enumerate edges instead, some edges won't be present in practice
        transmap.for_each_read_of_sample(sample_id, [&](int64_t id_a){
            cerr << id_a << ' ' << transmap.get_node(id_a).name << '\n';

            transmap.for_each_read_of_sample(sample_id, [&](int64_t id_b){
                if (id_a == id_b){
                    return;
                }

                auto [success,w] = transmap.try_get_edge_weight(id_a, id_b);

                if (not success){
                    return;
                }

                // Choosing any edge starting with id_a implies that id_a is parent and id_b is child
                pair<int64_t,int64_t> p = {id_a, id_b};
                auto& e = vars.edges.emplace(p, model.NewBoolVar()).first->second;

                cerr << '\t' <<  transmap.get_node(id_a).name << "->" << transmap.get_node(id_b).name << ' ' << w << '\n';

                // Keep track of which edges would make id_a a parent
                parent_edge_sums[id_a] += e;

                // Keep track of which edges would make id_b a child
                child_edge_sums[id_b] += e;

                // Cost of selecting an edge
                vars.cost_d += int64_t(w)*e;
                n_edges++;
            });

            n_reads++;
        });

        transmap.for_each_read_of_sample(sample_id, [&](int64_t id){
            auto& p = vars.is_parent.at(id);

            auto& p_sum = parent_edge_sums.at(id);
            auto& c_sum = child_edge_sums.at(id);

            // If there is at least one edge from id to another node, it is a parent node
            model.AddGreaterThan(p_sum, 0).OnlyEnforceIf(p);
            model.AddLessOrEqual(p_sum, 0).OnlyEnforceIf(Not(p));

            // Constraint that being a parent node excludes possibility of incoming child edges
            auto c = model.NewBoolVar();
            model.AddEquality(c_sum, 1).OnlyEnforceIf(c);
            model.AddEquality(c_sum, 0).OnlyEnforceIf(Not(c));
            model.AddImplication(p, Not(c));

            // Constraint that multiple edges cannot point to a child
            model.AddLessOrEqual(c_sum, 1);

            // Constraint that a node must be covered by an edge
            model.AddGreaterThan(c_sum + p_sum, 0);
        });

        cerr << "Number of reads: " << n_reads << '\n';
        cerr << "Number of edges: " << n_edges << '\n';

        model.Minimize(vars.cost_d);

        cerr << '\n';

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
                else{
                    results.insert(name + ": child");
                }
            }

            // Iterate the child assignment variables and print them
            for (const auto& [e,var]: vars.edges){
                auto [a,b] = e;
                auto name_a = transmap.get_node(a).name;
                auto name_b = transmap.get_node(b).name;
                auto [success, w] = transmap.try_get_edge_weight(a,b);

                bool is_assigned = (SolutionIntegerValue(response, var) == 1);
                cerr << is_assigned << ' ' << name_a << "->" << name_b << ' ' << w << '\n';
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
