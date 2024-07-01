#include <iostream>
#include <stdexcept>

using std::runtime_error;
using std::cerr;

#include "path_optimizer_mathopt.hpp"
#include "TransitiveMap.hpp"
#include "interval_tree.hpp"
#include "VariantGraph.hpp"
#include "VcfReader.hpp"
#include "fasta.hpp"
#include "Timer.hpp"
#include "CLI11.hpp"
#include "fetch.hpp"
#include "misc.hpp"
#include "gaf.hpp"
#include "bed.hpp"

using sv_merge::get_uuid;
using lib_interval_tree::interval_tree_t;

#include "bdsg/hash_graph.hpp"

using bdsg::HandleGraph;

#include <filesystem>

using std::filesystem::path;
using std::filesystem::create_directories;


#include <unordered_map>
#include <exception>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <thread>
#include <limits>
#include <cmath>
#include <tuple>

using std::this_thread::sleep_for;
using std::numeric_limits;
using std::unordered_map;
using std::runtime_error;
using std::exception;
using std::ifstream;
using std::ofstream;
using std::thread;
using std::atomic;
using std::mutex;
using std::log10;
using std::tuple;
using std::cerr;
using std::max;
using std::min;
using std::cref;
using std::ref;


using namespace sv_merge;


void for_each_row_in_csv(path csv_path, const function<void(const vector<string>& items)>& f){
    if (not (csv_path.extension() == ".csv")){
        throw runtime_error("ERROR: file does not have compatible csv extension: " + csv_path.string());
    }

    ifstream file(csv_path);

    if (not (file.is_open() and file.good())){
        throw runtime_error("ERROR: could not read file: " + csv_path.string());
    }

    char c;
    vector<string> items = {""};

    int64_t n_char_in_line = 0;
    char delimiter = ',';

    while (file.get(c)){
        if (c == delimiter){
            items.emplace_back();
            continue;
        }
        if (c == '\r'){
            throw runtime_error("ERROR: carriage return not supported: " + csv_path.string());
        }

        if (c == '\n'){
            if (n_char_in_line == 0){
                continue;
            }

            f(items);

            items.resize(1);
            items[0].clear();
            n_char_in_line = 0;
            continue;
        }

        items.back() += c;

        n_char_in_line++;
    }
}


void optimize(TransMap& transmap, const SolverType& solver_type){
    // Make tmp dir
    path output_dir = "/tmp/" + get_uuid();

    if (not exists(output_dir)){
        create_directories(output_dir);
    }
    else{
        throw runtime_error("Directory already exists: " + output_dir.string());
    }

    optimize_reads_with_d_and_n(transmap, 1, 1, 1, output_dir, solver_type);
}


void solve_from_csv(path csv, const SolverType& solver_type, size_t max_reads_per_sample){
    TransMap transmap;

    ifstream csv_file(csv);

    if (not csv_file.is_open() or csv_file.bad()){
        throw runtime_error("Could not open CSV file " + csv.string());
    }

    for_each_row_in_csv(csv, [&](const vector<string>& items){
        if (items.size() != 4){
            throw runtime_error("ERROR: CSV row does not have 3 items: " + csv.string());
        }

        const string& sample = items[0];
        const string& read = items[1];
        const string& path = items[2];
        float weight = stof(items[3]);

        if (not transmap.has_node(sample)){
//            cerr << "adding: " << sample << '\n';
            transmap.add_sample(sample);
        }

        if (not transmap.has_node(read)){
//            cerr << "adding: " << read << '\n';
            transmap.add_read(read);
        }

        if (not transmap.has_node(path)){
//            cerr << "adding: " << path << '\n';
            transmap.add_path(path);
        }

//        cerr << "adding edge: " << read << ',' << path << ',' << weight << '\n';
        transmap.add_edge(read, path, weight);
        transmap.add_edge(sample, read);
    });

    if (max_reads_per_sample != numeric_limits<size_t>::max()){
        transmap.for_each_sample([&](const string& sample_name, int64_t sample_id){
            size_t counter = 0;
            transmap.for_each_read_of_sample(sample_id, [&](int64_t read_id){
                if (counter >= max_reads_per_sample){
                    transmap.remove_edge(sample_id, read_id);
                }
                counter++;
            });
        });
    }

    optimize(transmap, solver_type);
}


void solve_from_csv_samplewise(path csv, const SolverType& solver_type){
    TransMap transmap;

    ifstream csv_file(csv);

    if (not csv_file.is_open() or csv_file.bad()){
        throw runtime_error("Could not open CSV file " + csv.string());
    }

    string prev_sample;
    vector <vector <string> > sample_data;

    for_each_row_in_csv(csv, [&](const vector<string>& items){
        if (items.size() != 4){
            throw runtime_error("ERROR: CSV row does not have 3 items: " + csv.string());
        }

        const string& sample = items[0];
        const string& read = items[1];
        const string& hap = items[2];
        float weight = stof(items[3]);

        if (sample != prev_sample and not sample_data.empty()){
            optimize(transmap, solver_type);
            if (transmap.empty()){
                cerr << "WARNING: no result for sample: " << prev_sample << '\n';
                for (const auto& item: sample_data){
                    for (const string& x: item){
                        cerr << x << ',';
                    }
                    cerr << '\n';
                }
            }

            transmap = {};
            sample_data.clear();
        }

        sample_data.emplace_back(items);

        if (not transmap.has_node(sample)){
//            cerr << "adding: " << sample << '\n';
            transmap.add_sample(sample);
        }

        if (not transmap.has_node(read)){
//            cerr << "adding: " << read << '\n';
            transmap.add_read(read);
        }

        if (not transmap.has_node(hap)){
//            cerr << "adding: " << hap << '\n';
            transmap.add_path(hap);
        }

//        cerr << "adding edge: " << read << ',' << hap << ',' << weight << '\n';
        transmap.add_edge(read, hap, weight);
        transmap.add_edge(sample, read);

        prev_sample = sample;
    });

    optimize(transmap, solver_type);
}


int main(int argc, char** argv){
    CLI::App app{"Solve from CSV"};

    path input_csv;
    string solver;
    SolverType solver_type;
    bool samplewise = false;
    size_t max_reads_per_sample = numeric_limits<size_t>::max();

    app.add_option("-i,--input", input_csv, "Input CSV file with sample-read-path data for optimizer")->required();
    app.add_option("--solver", solver, "Solver to use, must be one of: scip, glop, pdlp")->required();
    app.add_flag("-s,--samplewise", samplewise, "Optimize each sample separately (default: optimize all samples together)");
    app.add_option("-m,--max-reads", max_reads_per_sample, "Maximum number of reads to optimize per sample (default: all reads). Does NOT appy to samplewise optimization.");

    CLI11_PARSE(app, argc, argv);

    if (solver == "scip"){
        solver_type = SolverType::kGscip;
    }
    else if (solver == "glop"){
        solver_type = SolverType::kGlop;
    }
    else if (solver == "pdlp"){
        solver_type = SolverType::kPdlp;
    }
    else{
        throw runtime_error("ERROR: unknown solver: " + solver);
    }

    if (samplewise){
        solve_from_csv_samplewise(input_csv, solver_type);
    }
    else{
        solve_from_csv(input_csv, solver_type, max_reads_per_sample);
    }
}
