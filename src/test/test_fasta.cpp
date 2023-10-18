#include "Filesystem.hpp"
#include "Sequence.hpp"
#include "Region.hpp"
#include "fasta.hpp"

using sv_merge::for_sequence_in_fasta_file;
using ghc::filesystem::path;
using sv_merge::Sequence;
using sv_merge::Region;

#include <stdexcept>
#include <iostream>
#include <fstream>
#include <vector>

using std::runtime_error;
using std::to_string;
using std::ofstream;
using std::string;
using std::vector;
using std::cerr;


int main(){
    path project_directory = path(__FILE__).parent_path().parent_path().parent_path();
    path data_directory = project_directory / "data";

    cerr << data_directory << '\n';

    path fasta_path = data_directory / "test.fasta";

    vector<Sequence> expected_results = {
            Sequence("a", "AAAA"),
            Sequence("c", "CCCCCCCC"),
            Sequence("g", "GGGGGGGGGGGG")
    };

    size_t i=0;
    for_sequence_in_fasta_file(fasta_path,[&](const Sequence& r){
        cerr << r.name << '\t' << r.sequence << '\n';

        const auto& e = expected_results[i];

        if (r.name != e.name or r.sequence != e.sequence){
            throw runtime_error("FAIL: result not expected on line: " + to_string(i));
        }

        i++;
    });

    if (i != expected_results.size()){
        throw runtime_error("FAIL: not enough results iterated");
    }

    return 0;
}

