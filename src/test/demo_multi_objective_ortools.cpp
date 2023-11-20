#include <unordered_map>
#include <unordered_set>
#include <stdexcept>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <utility>
#include <random>
#include <string>
#include <vector>
#include <cmath>
#include <array>
#include <set>

using std::uniform_int_distribution;
using std::normal_distribution;
using std::runtime_error;
using std::unordered_map;
using std::unordered_set;
using std::ofstream;
using std::vector;
using std::string;
using std::tuple;
using std::array;
using std::hash;
using std::pair;
using std::cerr;
using std::pow;
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

using color_t = array<double,3>;


void generate_data(
        int64_t v_min,
        int64_t v_max,
        vector <pair <double,double> >& coords,
        vector <pair <double,double> >& offset_coords,
        vector <int>& labels
        ){

    std::default_random_engine generator;
    uniform_int_distribution offsets(v_min+1,v_max-1);

    vector<double> stdevs = {1.0, 1.0, 1.0};

    int64_t n = 8;

    int l = 0;
    for (auto s: stdevs){
        // Pick a random place to center the distribution
        auto x_offset = offsets(generator);
        auto y_offset = offsets(generator);

        offset_coords.emplace_back(x_offset, y_offset);

        // Construct the distribution with one of the provided sigmas
        normal_distribution<double> distribution_a(0,s);

        cerr << "Sampling distribution: " << '\n';
        cerr << "x:\t" << x_offset << '\n';
        cerr << "y:\t" << y_offset << '\n';
        cerr << "s:\t" << s << '\n';
        cerr << "n:\t" << n << '\n';
        cerr << '\n';

        for (int i=0; i<n; i++) {
            labels.emplace_back(l);

            double x = -1;
            double y = -1;

            // Draw x and y coords from the distribution, don't allow them to exist outside the bounds
            while (x<v_min or y>v_max){
                x = x_offset + distribution_a(generator);
            }
            while (y<v_min or y>v_max){
                y = y_offset + distribution_a(generator);
            }

            coords.emplace_back(x, y);
        }

        l++;
    }
}


double get_distance(const pair<double,double>& a, const pair<double,double>& b){
    auto t1 = pow(a.first - b.first, 2);
    auto t2 = pow(a.second - b.second, 2);
    return sqrt(t1 + t2);
}


/// Data class to contain the variables of the ILP, int64_t IDs intended to be reused from the TransMap
class Variables{
public:
    // Edges where the left node is the parent
    unordered_map <int64_t, unordered_map <int64_t,int64_t> > costs;
    unordered_map <int64_t, unordered_map <int64_t,BoolVar> > edges;

    // Nodes
    unordered_map <int64_t, BoolVar> is_parent;

    // Cost of all edges assigned
    LinearExpr cost_d;

    // Cost of adding parents
    LinearExpr cost_n;
};


void construct_model(const vector <pair <double,double> >& coords, CpModelBuilder& model, Variables& vars){
    for (int64_t i=0; i<coords.size(); i++) {
        vars.is_parent[i] = model.NewBoolVar();
        vars.cost_n += vars.is_parent[i];
    }

    unordered_map<int64_t, LinearExpr> child_edge_sums;
    unordered_map<int64_t, LinearExpr> parent_edge_sums;

    int64_t n_nodes = coords.size();
    int64_t n_edges = 0;

    for (int64_t id_a=0; id_a<coords.size(); id_a++){
        const auto& a = coords[id_a];

        for (int64_t id_b=0; id_b<coords.size(); id_b++){
            const auto& b = coords[id_b];

            if (a == b){
                continue;
            }

            // The cost of this edge is the literal euclidian distance, with 1 decimals of precision, as an integer
            auto d = get_distance(a,b);
            auto cost = int64_t(round(d*10));

            vars.costs[id_a][id_b] = cost;
            vars.edges[id_a][id_b] = model.NewBoolVar();

            // Keep track of which edges would make id_a a parent
            parent_edge_sums[id_a] += vars.edges[id_a][id_b];

            // Keep track of which edges would make id_b a child
            child_edge_sums[id_b] += vars.edges[id_a][id_b];

            // Cost of selecting an edge
            vars.cost_d += cost*vars.edges[id_a][id_b];
            n_edges++;
        }
    }

    for (int64_t id=0; id<coords.size(); id++){
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
    }

    cerr << "Number of nodes: " << n_nodes << '\n';
    cerr << "Number of edges: " << n_edges << '\n';
}


void jointly_optimize(CpModelBuilder& model, Variables& vars, CpSolverResponse& response){
    cerr << "Finding first pareto bound..." << '\n';

    // First find one extreme of the pareto set
    model.Minimize(vars.cost_d);

    const CpSolverResponse response_d = Solve(model.Build());

    int64_t n_max = SolutionIntegerValue(response_d, vars.cost_n);
    int64_t d_min = SolutionIntegerValue(response_d, vars.cost_d);

    cerr << "Finding second pareto bound..." << '\n';

    // Then find the other extreme of the pareto set
    model.Minimize(vars.cost_n);

    const CpSolverResponse response_n = Solve(model.Build());

    int64_t n_min = SolutionIntegerValue(response_n, vars.cost_n);
    int64_t d_max = SolutionIntegerValue(response_n, vars.cost_d);

    cerr << "Finding joint optimum..." << '\n';

    // Use pareto extremes to normalize the ranges of each objective and then jointly minimize distance from (n_min,d_min)
    auto n_range = n_max - n_min;
    n_range *= n_range;
    auto d_range = d_max - d_min;
    d_range *= d_range;

    auto d_norm = model.NewIntVar({0, d_max*d_max*n_range*2});
    auto n_norm = model.NewIntVar({0, d_max*d_max*n_range*2});

    model.AddEquality(d_norm, (vars.cost_d-d_min));
    model.AddEquality(n_norm, (vars.cost_n-n_min));

    auto d_square = model.NewIntVar({0, d_max*d_max*n_range*2});
    model.AddMultiplicationEquality(d_square,{d_norm,d_norm});

    auto n_square = model.NewIntVar({0, d_max*d_max*n_range*2});
    model.AddMultiplicationEquality(n_square,{n_norm,n_norm});

    model.Minimize(n_square*d_range + d_square*n_range);

//    SatParameters parameters;
//    parameters.set_log_search_progress(true);
//    const CpSolverResponse response_n_d = SolveWithParameters(model.Build(), parameters);

    response = Solve(model.Build());
}


void plot(
        size_t precision,
        int64_t v_min,
        int64_t v_max,
        const vector <pair <double,double> >& coords,
        const vector<color_t>& colors,
        string filename
        ){
    ofstream file(filename);

    if (not (file.good() and file.is_open())){
        throw runtime_error("ERROR: could not write to file: " + filename);
    }

    unordered_map <int64_t, unordered_map<int64_t, int64_t> > coord_map;

    for (int64_t i=0; i<coords.size(); i++){
        auto [x,y] = coords[i];

        auto x_pixel = int64_t(round(x*pow(10,precision)));
        auto y_pixel = int64_t(round(y*pow(10,precision)));

        coord_map[x_pixel][y_pixel] = i;
    }

    auto size = (v_max - v_min)*pow(10,precision);

    color_t white = {255,255,255};

    file << "P3" << '\n';
    file << size << ' ' << size << '\n';
    file << "255" << '\n';
    for (size_t i=0; i<size; i++){
        for (size_t j=0; j<size; j++){
            auto result = coord_map.find(i);
            color_t color;

            if (result != coord_map.end()){
                auto result2 = result->second.find(j);
                if (result2 != result->second.end()){
                    auto c = result2->second;
                    color = colors[c];
                }
                else{
                    color = white;
                }
            }
            else{
                color = white;
            }

            file << int(color[0]) << ' ' << int(color[1]) << ' ' << int(color[2]) << "    ";
        }
        file << '\n';
    }
}


void parse_solution(
        int64_t v_min,
        int64_t v_max,
        vector <pair <double,double> >& coords,
        vector <int>& labels,
        Variables& vars,
        CpSolverResponse& response
        ){

    vector<color_t> color_key = {
            {255,140,0},
            {154,205,50},
            {72,61,139},
            {0,139,139},
            {0,0,139},
            {0,100,0},
            {139,0,139},
            {255,255,0},
            {124,252,0},
            {0,255,127},
            {220,20,60},
            {0,191,255},
            {255,0,255},
            {30,144,255},
            {219,112,147},
            {240,230,140},
            {176,224,230},
            {255,20,147},
            {255,160,122},
            {238,130,238},
            {127,255,212}
    };

    vector<color_t> colors;

    for (auto& l: labels){
        colors.push_back(color_key[l]);
    }

    plot(1, v_min, v_max, coords, colors, "input.ppm");

    unordered_map<int64_t,int64_t> id_to_parent;

    ofstream nodes_file("nodes.csv");
    ofstream edges_file("edges.csv");

    nodes_file << "id,x,y,is_parent,true_label" << '\n';
    edges_file << "id_a,x_a,y_a,id_b,x_b,y_b" << '\n';

    if (response.status() == CpSolverStatus::OPTIMAL || response.status() == CpSolverStatus::FEASIBLE) {
        cerr << "Maximum of objective function: " << response.objective_value() << '\n';

        // Iterate the parent assignment variables and print them
        for (const auto& [id,var]: vars.is_parent){
            bool is_parent = (SolutionIntegerValue(response, var) == 1);
            nodes_file << id << ',' << coords[id].first << ',' << coords[id].second << ',' << is_parent << ',' << labels[id] << '\n';

            if (is_parent){
                auto p = id_to_parent.size();
                id_to_parent[id] = p;
            }
        }

        // Iterate the edges and plot/print the solution
        for (const auto& [a,item]: vars.edges){
            auto color = color_key[id_to_parent[a]];

            for (const auto& [b,var]: item){
                bool is_assigned = (SolutionIntegerValue(response, var) == 1);
                if (is_assigned) {
                    cerr << a << "->" << b << '\n';
                    edges_file << a << ',' << coords[a].first << ',' << coords[a].second << b << ',' << coords[b].first << ',' << coords[b].second << '\n';
                    colors[a] = color;
                    colors[b] = color;
                }
            }
        }
    }
    else {
        cerr << "No solution found." << '\n';
    }

    cerr << '\n';
    cerr << "Statistics" << '\n';
    cerr << CpSolverResponseStats(response) << '\n';

    plot(1, v_min, v_max, coords, colors, "result.ppm");
}


int main(){
    vector <pair <double,double> > coords;
    vector <pair <double,double> > offset_coords;
    vector <int> labels;
    Variables vars;
    CpModelBuilder model;
    CpSolverResponse response;

    int64_t v_min = 0;
    int64_t v_max = 10;

    cerr << "Generating data..." << '\n';
    generate_data(v_min, v_max, coords, offset_coords, labels);

    cerr << "Constructing model..." << '\n';
    construct_model(coords, model, vars);

    cerr << "Optimizing..." << '\n';
    jointly_optimize(model, vars, response);

    cerr << "Generating output..." << '\n';
    parse_solution(v_min, v_max, coords, labels, vars, response);

    cerr << '\n';
}

