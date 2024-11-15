#include <iostream>
#include <stdexcept>

using std::runtime_error;
using std::cerr;

#include "path_optimizer_mathopt.hpp"
#include "TransitiveMap.hpp"
#include "interval_tree.hpp"
#include "VcfReader.hpp"
#include "fasta.hpp"
#include "Timer.hpp"
#include "CLI11.hpp"
#include "fetch.hpp"
#include "misc.hpp"
#include "gaf.hpp"

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


void solve_from_csv(path csv, const OptimizerConfig& config, size_t max_reads_per_sample){
    TransMap transmap;

    ifstream csv_file(csv);

    if (not csv_file.is_open() or csv_file.bad()){
        throw runtime_error("Could not open CSV file " + csv.string());
    }

    bool skip_header = true;

    for_each_row_in_csv(csv, [&](const vector<string>& items){
        if (skip_header) {
            skip_header = false;
            return;
        }

        if (items.size() != 6){
            throw runtime_error("ERROR: CSV row does not have 6 items: " + csv.string());
        }

        // sample,read,read_length,path,path_length,weight
        const string& sample = items[0];
        const string& read = items[1];
        const string& path = items[3];
        float weight = stof(items[5]);

        if (not transmap.has_node(sample)){
            transmap.add_sample(sample);
        }

        if (not transmap.has_node(read)){
            transmap.add_read(read);
        }

        if (not transmap.has_node(path)){
            transmap.add_path(path);
        }

        transmap.add_edge(read, path, weight);
        transmap.add_edge(sample, read);
    });

    vector <pair <int64_t, int64_t> > edges_to_remove;

    if (max_reads_per_sample != numeric_limits<size_t>::max()){
        transmap.for_each_sample([&](const string& sample_name, int64_t sample_id){
            size_t counter = 0;
            transmap.for_each_read_of_sample(sample_id, [&](int64_t read_id){
                if (counter >= max_reads_per_sample){
                    edges_to_remove.emplace_back(sample_id, read_id);
                }
                counter++;
            });
        });
    }

    for (const auto& edge: edges_to_remove){
        transmap.remove_edge(edge.first, edge.second);
    }

    // Make tmp dir
    path output_dir = "/tmp/" + get_uuid();

    if (not exists(output_dir)){
        create_directories(output_dir);
        cerr << output_dir << '\n';
    }
    else{
        throw runtime_error("Directory already exists: " + output_dir.string());
    }

    optimize(transmap, config, output_dir);
}


int main(int argc, char** argv){
    CLI::App app{"Solve from CSV"};

    path input_csv;
    string solver;
    size_t max_reads_per_sample = numeric_limits<size_t>::max();
    size_t n_threads = 1;

    path windows_bed;
    path tandem_bed;
    string bam_csv;
    path ref_fasta;
    path vcf;

    OptimizerConfig optimizer_config;

    app.add_option("--d_weight", optimizer_config.d_weight, "Scalar coefficient to apply to d_norm^2 in the optimization step, used to priotize or deprioritize the distance term");

    app.add_flag("--rescale_weights", optimizer_config.rescale_weights, "Invoke this to use quadratic difference-from-best match rescaling for read-to-path edges.");

    app.add_flag("--quadratic_objective", optimizer_config.use_quadratic_objective, "Invoke this to use quadratic objective which minimizes the square distance from the 'utopia point'. May incur large run time cost.");

    app.add_flag("--samplewise", optimizer_config.samplewise, "Use samplewise solver instead of global solver");

    app.add_flag("--prune_with_d_min", optimizer_config.prune_with_d_min, "Use d_min solution to remove all edges not used");

    app.add_option("-i,--input", input_csv, "Input CSV file with sample-read-path data for optimizer")->required();

    app.add_option("--solver", solver, "Solver to use, must be one of: scip, glop, pdlp")->required();

    app.add_flag("!--no_ploidy", optimizer_config.use_ploidy_constraint, "If invoked, do not enforce a ploidy <= 2 constraint per sample (w.r.t. paths).");

    app.add_flag("-g,--use_golden_search", optimizer_config.use_golden_search, "If invoked, use explicit n-centric optimization with golden search to find joint minimum instead of encoding joint objective in the solver.");

    app.add_option("-m,--max-reads", max_reads_per_sample, "Maximum number of reads to optimize per sample (default: all reads). Does NOT appy to samplewise optimization.");

    app.add_option("-t,--n_threads", n_threads, "Maximum number of threads to use for solver (default: 1)");

    CLI11_PARSE(app, argc, argv);

    if (solver == "scip"){
        optimizer_config.solver_type = SolverType::kGscip;
    }
    else if (solver == "glop"){
        optimizer_config.solver_type = SolverType::kGlop;
    }
    else if (solver == "pdlp"){
        optimizer_config.solver_type = SolverType::kPdlp;
    }
    else if (solver == "gurobi"){
        optimizer_config.solver_type = SolverType::kGurobi;
    }
    else{
        throw runtime_error("ERROR: unknown solver: " + solver);
    }

    if (optimizer_config.use_ploidy_constraint){
        cerr << "Using ploidy constraint...\n";
    }
    else{
        cerr << "NOT using ploidy constraint...\n";
    }

    solve_from_csv(input_csv, optimizer_config, max_reads_per_sample);
}
