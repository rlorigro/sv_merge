#pragma once

#include "HeteroGraph.hpp"

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


class TransMap {
    /// Attributes
    HeteroGraph<HeteroNode> graph;

    // Pains me to add yet another map but here it is
    unordered_map<int64_t,string> sequences;

    const string sample_node_name;
    const string read_node_name;
    const string path_node_name;

public:
    TransMap();

    /// Building
    void add_sample(const string& name);
    void add_read(const string& name);
    void add_read(const string& name, const string& sequence);
    void add_path(const string& name);
    void add_edge(const string& a, const string& b);
    void add_edge(const string& a, const string& b, float weight);
    void construct_optimizer();

    /// Accessing
    const HeteroNode& get_node(int64_t id) const;

    void for_each_sample(const function<void(const string& name, int64_t id)>& f) const;
    void for_each_read(const function<void(const string& name, int64_t id)>& f) const;
    void for_each_path(const function<void(const string& name, int64_t id)>& f) const;

    void get_read_sample(const string& read_name, string& result) const;
    void get_read_sample(int64_t read_id, string& result) const;

    void for_each_read_of_sample(const string& sample_name, const function<void(const string& name, int64_t id)>& f) const;
    void for_each_read_of_sample(int64_t sample_id, const function<void(int64_t read_id)>& f) const;
    void for_each_read_of_path(int64_t path_id, const function<void(int64_t id)>& f) const;

    void for_each_sample_of_read(const string& read_name, const function<void(const string& name, int64_t id)>& f) const;
    void for_each_sample_of_path(const string& path_name, const function<void(const string& name, int64_t id)>& f) const;

    void for_each_path_of_sample(const string& sample_name, const function<void(const string& name, int64_t id)>& f) const;

    void for_each_path_of_read(int64_t read_id, const function<void(int64_t path_id)>& f) const;

    void for_each_read_to_path_edge(const function<void(int64_t read_id, int64_t path_id, float weight)>& f) const;
};


TransMap::TransMap():
    sample_node_name("sample_node"),
    read_node_name("read_node"),
    path_node_name("path_node")
{
    graph.add_node(sample_node_name, 'S');
    graph.add_node(read_node_name, 'R');
    graph.add_node(path_node_name, 'P');
}


const HeteroNode& TransMap::get_node(int64_t id) const{
    return graph.get_node(id);
}


void TransMap::add_read(const string& name){
    if (name == read_node_name or name == sample_node_name or name == path_node_name){
        throw runtime_error("ERROR: cannot add node with preset node name: " + name);
    }
    graph.add_node(name, 'R');
    graph.add_edge(read_node_name, name, 0);
}


void TransMap::add_read(const string& name, const string& sequence){
    if (name == read_node_name or name == sample_node_name or name == path_node_name){
        throw runtime_error("ERROR: cannot add node with preset node name: " + name);
    }
    graph.add_node(name, 'R');
    graph.add_edge(read_node_name, name, 0);

    sequences.emplace(graph.name_to_id(name), sequence);
}


void TransMap::add_sample(const string& name){
    if (name == read_node_name or name == sample_node_name or name == path_node_name){
        throw runtime_error("ERROR: cannot add node with preset node name: " + name);
    }
    graph.add_node(name, 'S');
    graph.add_edge(sample_node_name, name, 0);
}


void TransMap::add_path(const string& name){
    if (name == read_node_name or name == sample_node_name or name == path_node_name){
        throw runtime_error("ERROR: cannot add node with preset node name: " + name);
    }
    graph.add_node(name, 'P');
    graph.add_edge(path_node_name, name, 0);
}


void TransMap::add_edge(const string& a, const string& b){
    graph.add_edge(a,b,0);
}


void TransMap::add_edge(const string& a, const string& b, float weight){
    graph.add_edge(a,b,weight);
}


void TransMap::get_read_sample(const string& read_name, string& result) const{
    result.clear();

    // Using the HeteroGraph back end means that we iterate, even for a 1:1 mapping.
    graph.for_each_neighbor_of_type(read_name, 'S', [&](const HeteroNode& neighbor, int64_t id){
        // Check for cases that should be impossible
        if (not result.empty()){
            throw runtime_error("ERROR: multiple samples found for read: " + read_name);
        }

        result = neighbor.name;
    });

    if (result.empty()){
        throw runtime_error("ERROR: no sample found for read: " + read_name);
    }
}


void TransMap::get_read_sample(int64_t read_id, string& result) const{
    result.clear();

    // Using the HeteroGraph back end means that we iterate, even for a 1:1 mapping.
    graph.for_each_neighbor_of_type(read_id, 'S', [&](int64_t id){
        // Check for cases that should be impossible
        if (not result.empty()){
            throw runtime_error("ERROR: multiple samples found for read id: " + read_id);
        }

        result = graph.get_node(id).name;
    });

    if (result.empty()){
        throw runtime_error("ERROR: no sample found for read id: " + read_id);
    }
}


void TransMap::for_each_read(const function<void(const string& name, int64_t id)>& f) const{
    graph.for_each_neighbor_of_type(read_node_name, 'R', [&](const HeteroNode& neighbor, int64_t id){
        f(neighbor.name, id);
    });
}


void TransMap::for_each_sample(const function<void(const string& name, int64_t id)>& f) const{
    graph.for_each_neighbor_of_type(sample_node_name, 'S', [&](const HeteroNode& neighbor, int64_t id){
        f(neighbor.name, id);
    });
}


void TransMap::for_each_path(const function<void(const string& name, int64_t id)>& f) const{
    graph.for_each_neighbor_of_type(path_node_name, 'P', [&](const HeteroNode& neighbor, int64_t id){
        f(neighbor.name, id);
    });
}


void TransMap::for_each_read_of_sample(const string& sample_name, const function<void(const string& name, int64_t id)>& f) const{
    graph.for_each_neighbor_of_type(sample_name, 'R', [&](const HeteroNode& neighbor, int64_t id){
        f(neighbor.name, id);
    });
}


void TransMap::for_each_read_of_sample(int64_t sample_id, const function<void(int64_t read_id)>& f) const{
    graph.for_each_neighbor_of_type(sample_id, 'R', [&](int64_t id){
        f(id);
    });
}


void TransMap::for_each_path_of_read(int64_t read_id, const function<void(int64_t path_id)>& f) const{
    graph.for_each_neighbor_of_type(read_id, 'P', [&](int64_t id){
        f(id);
    });
}


void TransMap::for_each_read_of_path(int64_t path_id, const function<void(int64_t id)>& f) const{
    graph.for_each_neighbor_of_type(path_id, 'R', [&](int64_t r){
        f(r);
    });
}


void TransMap::for_each_sample_of_read(const string& read_name, const function<void(const string& name, int64_t id)>& f) const{
    graph.for_each_neighbor_of_type(read_name, 'S', [&](const HeteroNode& neighbor, int64_t id){
        f(neighbor.name, id);
    });
}


void TransMap::for_each_sample_of_path(const string& path_name, const function<void(const string& name, int64_t id)>& f) const{
    auto id = graph.name_to_id(path_name);
    graph.for_each_neighbor_of_type(id, 'R', [&](int64_t r){
        graph.for_each_neighbor_of_type(r, 'S', [&](int64_t s){
            f(graph.get_node(s).name, s);
        });
    });
}


void TransMap::for_each_path_of_sample(const string& sample_name, const function<void(const string& name, int64_t id)>& f) const{
    auto id = graph.name_to_id(sample_name);
    graph.for_each_neighbor_of_type(id, 'R', [&](int64_t r){
        graph.for_each_neighbor_of_type(r, 'P', [&](int64_t p){
            f(graph.get_node(p).name, p);
        });
    });
}


void TransMap::for_each_read_to_path_edge(const function<void(int64_t read_id, int64_t path_id, float weight)>& f) const{
    // Starting from the source node which connects to all reads, find all neighbors (read nodes)
    graph.for_each_neighbor_of_type(read_node_name, 'R', [&](const HeteroNode& n, int64_t read_id){
        // Find all path-type neighbors
        graph.for_each_neighbor_of_type(read_id, 'P', [&](int64_t path_id, float w){
            f(read_id, path_id, w);
        });
    });
}


//void TransMap::construct_optimizer(){
//    CpModelBuilder model;
//
//    // TODO: reserve these data structures based on how many edges there are known to be?
//    unordered_map <pair<int64_t,int64_t>, BoolVar> path_to_read_variables;
//    unordered_map <int64_t, BoolVar> path_indicators;
//
//    LinearExpr cost_d;
//    LinearExpr cost_n;
//
//    // Read->path boolean indicators
//    for_each_read_to_path_edge([&](int64_t read_id, int64_t path_id, float weight){
//        pair<int64_t,int64_t> p = {path_id, read_id};
//
//        // Update the boolean var map, and keep a reference to the inserted BoolVar
//        const auto result = path_to_read_variables.emplace(p,model.NewBoolVar());
//        const auto& v = result.first->second;
//
//        auto w = int64_t(weight);
//        cost_d += v*w;
//    });
//
//    for_each_sample([&](const string& name, int64_t sample_id){
//        LinearExpr samplewise_vars;
//
//        // Sample->read boolean indicators
//        for_each_read_of_sample(sample_id, [&](int64_t read_id){
//            LinearExpr readwise_vars;
//
//            // Read->path boolean indicators
//            for_each_path_of_read(read_id, [&](int64_t path_id){
//                const BoolVar& var = path_to_read_variables.at({path_id, read_id});
//                samplewise_vars += var;
//                readwise_vars += var;
//            });
//
//            // Enforce (r == 1).
//            model.AddEquality(readwise_vars, 1);
//        });
//
//        // Enforce (s <= 2).
//        model.AddLessOrEqual(samplewise_vars, 2);
//    });
//
//    // Keep track of whether each path is used
//    for_each_path([&](const string& name, int64_t path_id){
//        LinearExpr pathwise_reads;
//        for_each_read_of_path(path_id, [&](int64_t read_id){
//            pathwise_reads += path_to_read_variables.at({path_id,read_id});
//        });
//
//        const auto result = path_indicators.emplace(path_id,model.NewBoolVar());
//        const auto& p = result.first->second;
//
//        // Implement p == (sum(r) > 0).
//        model.AddGreaterThan(pathwise_reads, 0).OnlyEnforceIf(p);
//        model.AddLessOrEqual(pathwise_reads, 0).OnlyEnforceIf(Not(p));
//    });
//
//    for (const auto& [path_id,p]: path_indicators){
//        cost_n += p;
//    }
//
//    model.Minimize(cost_d);
//
//    // TODO: Move elsewhere
//    const CpSolverResponse response = Solve(model.Build());
//
//    if (response.status() == CpSolverStatus::OPTIMAL ||
//        response.status() == CpSolverStatus::FEASIBLE) {
//
//        cerr << "Maximum of objective function: " << response.objective_value() << '\n';
//
//        // Iterate the read assignment variables and print them
//        for (auto& [item, var]: path_to_read_variables){
//            int64_t path_id = item.first;
//            int64_t read_id = item.second;
//
//            auto path_name = graph.get_node(path_id).name;
//            auto read_name = graph.get_node(read_id).name;
//            string sample_name;
//            get_read_sample(read_id, sample_name);
//
//            cerr << sample_name << ',' << path_name << ',' << read_name << " = " << SolutionIntegerValue(response, var) << '\n';
//        }
//    }
//    else {
//        cerr << "No solution found." << '\n';
//    }
//
//    cerr << "Statistics" << '\n';
//    cerr << CpSolverResponseStats(response) << '\n';
//}


}
