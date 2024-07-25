#include <string>
#include <iostream>
#include <stdexcept>
#include <unordered_map>

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
#include <fstream>

using std::filesystem::path;
using std::filesystem::create_directories;
using std::filesystem::recursive_directory_iterator;


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


void optimize(TransMap& transmap, const SolverType& solver_type, size_t n_threads, bool use_ploidy_constraint, path output_dir){
    if (not exists(output_dir)){
        create_directories(output_dir);
    }
    else{
        throw runtime_error("Directory already exists: " + output_dir.string());
    }

    optimize_reads_with_d_and_n(transmap, 1, 1, n_threads, output_dir, solver_type, use_ploidy_constraint);
}


size_t solve_from_csv(
        path csv,
        const SolverType& solver_type,
        size_t max_reads_per_sample,
        size_t n_threads,
        bool use_ploidy_constraint
        ){

    TransMap transmap;

    ifstream csv_file(csv);

    if (not csv_file.is_open() or csv_file.bad()){
        throw runtime_error("Could not open CSV file " + csv.string());
    }

    size_t i = 0;
    for_each_row_in_csv(csv, [&](const vector<string>& items){
        if (i == 0){
            i++;
            return;
        }

        if (items.size() != 6){
            throw runtime_error("ERROR: CSV row does not have 6 items: " + csv.string());
        }

        //sample,read,read_length,path,path_length,weight
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
        i++;
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

    // Construct a random seed that is based on the csv file
    auto h = std::hash<std::string>{}(csv.string());

    // Make tmp dir
    path output_dir = "/tmp/" + to_string(h);

    optimize(transmap, solver_type, 1, use_ploidy_constraint, output_dir);

    return i-1;
}


void thread_fn(
        atomic<size_t>& job_index,
        const vector<path>& jobs,
        mutex& io_mutex,
        const SolverType& solver_type,
        size_t max_reads_per_sample,
        size_t n_threads,
        bool use_ploidy_constraint,
        path output_dir
        ){

    size_t i = job_index.fetch_add(1);
    path log_path = output_dir / "log.csv";

    while (i < jobs.size()){
        Timer t;
        size_t n_vars = 0;
        bool success = true;

        try {
            n_vars = solve_from_csv(jobs[i], solver_type, max_reads_per_sample, n_threads, use_ploidy_constraint);
        }
        catch (const exception& e) {
            cerr << e.what() << '\n';
            cerr << "ERROR caught at " << jobs[i].string() << '\n';
            success = false;
        }

        io_mutex.lock();
        ofstream log(log_path, std::ios_base::app);

        // get directory name of job (excluding non-leaf dirs)
        string job_dir = jobs[i].parent_path().filename();

        log << job_dir << "," << n_vars << ',' << success << ',' << t.to_csv() << "\n";
        log.close();

        io_mutex.unlock();

        i = job_index.fetch_add(1);
    }
}


// A function which opens a directory and iterates all subdirectories, calling solve_from_csv on each reads_to_paths.csv
void solve_from_directory(
        path directory,
        const SolverType& solver_type,
        size_t max_reads_per_sample,
        size_t n_threads,
        bool use_ploidy_constraint,
        path output_dir
        ){

    if (not exists(output_dir)){
        create_directories(output_dir);
    }
    else{
        throw runtime_error("ERROR: output directory already exists: " + output_dir.string());
    }

    if (not exists(directory)){
        throw runtime_error("ERROR: input directory does not exist: " + directory.string());
    }

    // Find all reads_to_paths.csv files in the subdirectories of the input directory and add them to a vector of jobs
    vector<path> jobs;
    for (const auto& entry : recursive_directory_iterator(directory)){
        if (entry.path().filename() == "reads_to_paths.csv"){
            jobs.push_back(entry.path());
        }
    }

    // Launch n_threads using the given args and the jobs vector
    vector<thread> threads;
    atomic<size_t> job_index(0);
    mutex io_mutex;

    // Write header to log csv
    ofstream log(output_dir / "log.csv");
    log << "name,n_path_to_read_vars,success,h,m,s,ms\n";
    log.close();

    for (size_t i = 0; i < n_threads; i++){
        threads.emplace_back(thread_fn,
                ref(job_index),
                cref(jobs),
                ref(io_mutex),
                solver_type,
                max_reads_per_sample,
                n_threads,
                use_ploidy_constraint,
                output_dir
            );
    }

    for (auto& t : threads){
        t.join();
    }


}


int main(int argc, char** argv){
    CLI::App app{"Solve from CSV"};

    path input_csv;
    path output_dir;
    string solver;
    SolverType solver_type;
    bool use_ploidy_constraint = true;
    size_t max_reads_per_sample = numeric_limits<size_t>::max();
    size_t n_threads = 1;

    app.add_option("-i,--input", input_csv, "Input CSV file with sample-read-path data for optimizer")->required();
    app.add_option("-o,--output_dir", output_dir, "Output directory (must not exist)")->required();
    app.add_option("--solver", solver, "Solver to use, must be one of: scip, glop, pdlp")->required();
    app.add_flag("!--no_ploidy", use_ploidy_constraint, "If invoked, do not enforce a ploidy <= 2 constraint per sample (w.r.t. paths).");
    app.add_option("-m,--max-reads", max_reads_per_sample, "Maximum number of reads to optimize per sample (default: all reads). Does NOT appy to samplewise optimization.");
    app.add_option("-t,--n_threads", n_threads, "Maximum number of threads to use for solver (default: 1)");

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

    if (use_ploidy_constraint){
        cerr << "Using ploidy constraint...\n";
    }
    else{
        cerr << "NOT using ploidy constraint...\n";
    }

    solve_from_directory(input_csv, solver_type, max_reads_per_sample, n_threads, use_ploidy_constraint, output_dir);
}
