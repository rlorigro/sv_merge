#include <iostream>
#include <stdexcept>
#include <filesystem>
#include <string>

using std::ofstream;
using std::string;
using std::runtime_error;
using std::filesystem::path;
using std::filesystem::create_directories;
using std::filesystem::exists;
using std::exception;
using std::cerr;

#include "TransitiveMap.hpp"
#include "path_optimizer_mathopt.hpp"
#include "misc.hpp"

using sv_merge::get_uuid;
using sv_merge::TransMap;
using sv_merge::HeteroNode;
using sv_merge::optimize_reads_with_d_and_n;


void test_optimization(){
    TransMap transmap;

    transmap.add_sample("HG001");
    transmap.add_sample("HG002");

    transmap.add_read("read_01");
    transmap.add_read("read_02");
    transmap.add_read("read_03");
    transmap.add_read("read_04");
    transmap.add_read("read_05");
    transmap.add_read("read_06");

    transmap.add_path("a");
    transmap.add_path("b");
    transmap.add_path("c");

    transmap.add_edge("read_01", "HG001");
    transmap.add_edge("read_03", "HG001");
    transmap.add_edge("read_05", "HG001");

    transmap.add_edge("read_02", "HG002");
    transmap.add_edge("read_04", "HG002");
    transmap.add_edge("read_06", "HG002");

    // HG001 (odd)
    transmap.add_edge("read_01", "a", 1);
    transmap.add_edge("read_01", "b", 2);
    transmap.add_edge("read_01", "c", 3);

    transmap.add_edge("read_03", "a", 1);
    transmap.add_edge("read_03", "b", 2);
    transmap.add_edge("read_03", "c", 3);

    transmap.add_edge("read_05", "a", 1);
    transmap.add_edge("read_05", "b", 2);
    transmap.add_edge("read_05", "c", 3);

    // HG002 (even)
    transmap.add_edge("read_02", "a", 3);
    transmap.add_edge("read_02", "b", 2);
    transmap.add_edge("read_02", "c", 1);

    transmap.add_edge("read_04", "a", 2);
    transmap.add_edge("read_04", "b", 1);
    transmap.add_edge("read_04", "c", 3);

    transmap.add_edge("read_06", "a", 1);
    transmap.add_edge("read_06", "b", 2);
    transmap.add_edge("read_06", "c", 3);

    // Make tmp dir
    path output_dir = "/tmp/" + get_uuid();

    if (not exists(output_dir)){
        create_directories(output_dir / "conventional");
        create_directories(output_dir / "golden_search");
    }
    else{
        throw runtime_error("Directory already exists: " + output_dir.string());
    }

    vector <pair <SolverType,string> > solver_type = {
            {SolverType::kGscip,"kGscip"},
//            {SolverType::kGurobi,"kGurobi"},
//            {SolverType::kGlop,"kGlop"},
            {SolverType::kCpSat,"kCpSat"},
            {SolverType::kPdlp,"kPdlp"},
//            {SolverType::kGlpk,"kGlpk"},
//            {SolverType::kEcos,"kEcos"},
//            {SolverType::kScs,"kScs"},
//            {SolverType::kHighs,"kHighs"},
//            {SolverType::kSantorini,"kSantorin"}
    };

    for (const auto& [t,name]: solver_type) {
        cerr << "-------- solver: " << name << " --------" << '\n';

        try {
            TransMap transmap_copy = transmap;
            optimize_reads_with_d_and_n(transmap_copy, 1, 1, 1, output_dir/"conventional", t);
            cerr << '\n';

            transmap_copy = transmap;
            optimize_reads_with_d_and_n(transmap_copy, 1, 10, 1, output_dir/"conventional", t);
            cerr << '\n';

            transmap_copy = transmap;
            optimize_reads_with_d_and_n(transmap_copy, 10, 1, 1, output_dir/"conventional", t);
            cerr << '\n';
        }
        catch (const exception& e){
            cerr << "ERROR: " << e.what() << '\n';
        }

        cerr << '\n';
    }


    cerr << "-- golden search\n";
    try {
        TransMap transmap_copy = transmap;
        optimize_reads_with_d_and_n_using_golden_search(transmap_copy, 1, 1, 1, output_dir/"golden_search", SolverType::kGscip);
    }
    catch (const exception& e){
        cerr << "ERROR: " << e.what() << '\n';
    }
    cerr << '\n';

    try {
        TransMap transmap_copy = transmap;
        optimize_reads_with_d_and_n_using_golden_search(transmap_copy, 1, 10, 1, output_dir/"golden_search", SolverType::kGscip);
    }
    catch (const exception& e){
        cerr << "ERROR: " << e.what() << '\n';
    }
    cerr << '\n';

    try {
        TransMap transmap_copy = transmap;
        optimize_reads_with_d_and_n_using_golden_search(transmap_copy, 10, 1, 1, output_dir/"golden_search", SolverType::kGscip);
    }
    catch (const exception& e){
        cerr << "ERROR: " << e.what() << '\n';
    }
    cerr << '\n';

}


void test_infeasible_by_ploidy(){
    TransMap transmap;

    transmap.add_sample("HG001");
    transmap.add_sample("HG002");

    transmap.add_read("read_01");
    transmap.add_read("read_02");
    transmap.add_read("read_03");
    transmap.add_read("read_04");
    transmap.add_read("read_05");
    transmap.add_read("read_06");

    transmap.add_path("a");
    transmap.add_path("b");
    transmap.add_path("c");

    transmap.add_edge("read_01", "HG001");
    transmap.add_edge("read_03", "HG001");
    transmap.add_edge("read_05", "HG001");

    transmap.add_edge("read_02", "HG002");
    transmap.add_edge("read_04", "HG002");
    transmap.add_edge("read_06", "HG002");

    // HG001 (odd)
    transmap.add_edge("read_01", "a", 1);

    transmap.add_edge("read_03", "b", 2);

    transmap.add_edge("read_05", "c", 3);

    // Make tmp dir
    path output_dir = "/tmp/" + get_uuid();

    if (not exists(output_dir)){
        create_directories(output_dir);
    }
    else{
        throw runtime_error("Directory already exists: " + output_dir.string());
    }

    vector <pair <SolverType,string> > solver_type = {
            {SolverType::kGscip,"kGscip"},
            {SolverType::kGurobi,"kGurobi"},
            {SolverType::kGlop,"kGlop"},
            {SolverType::kCpSat,"kCpSat"},
            {SolverType::kPdlp,"kPdlp"},
            {SolverType::kGlpk,"kGlpk"},
            {SolverType::kEcos,"kEcos"},
            {SolverType::kScs,"kScs"},
            {SolverType::kHighs,"kHighs"},
            {SolverType::kSantorini,"kSantorin"}
    };

    for (const auto& [t,name]: solver_type) {
        cerr << "solver: " << name << '\n';

        try {
            TransMap transmap_copy = transmap;
            optimize_reads_with_d_and_n(transmap_copy, 1, 1, 1, output_dir, t);
        }
        catch (const exception& e){
            cerr << "ERROR: " << e.what() << '\n';
        }
    }

}


int main(){
    test_optimization();
//    test_infeasible_by_ploidy();
}
