#include <filesystem>
using std::filesystem::path;

#include "Region.hpp"
#include "bed.hpp"

using sv_merge::for_region_in_bed_file;
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

    path bed_path = data_directory / "test.bed";

    vector<Region> expected_results = {
            Region("chr1",1,100),
            Region("chr1",1,101),
            Region("chr2",2,200),
            Region("chr2",2,202)
    };

    size_t i=0;
    for_region_in_bed_file(bed_path,[&](const Region& r){
        cerr << r.name << '\t' << r.start << '\t' << r.stop << '\n';

        const auto& e = expected_results[i];

        if (r.name != e.name or r.start != e.start or r.stop != e.stop){
            throw runtime_error("FAIL: result not expected on line: " + to_string(i));
        }

        i++;
    });

    if (i != expected_results.size()){
        throw runtime_error("FAIL: not enough results iterated");
    }

    return 0;
}

