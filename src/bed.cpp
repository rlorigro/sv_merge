#include "bed.hpp"

#include <iostream>

using std::ifstream;
using std::cerr;


namespace sv_merge{


/**
 * Iterate over the regions of a bed file, ignoring annotation fields. Must end in new line if last line is to be
 * considered. Must use actual tab characters (not spaces)
 * @param bed_path path to bed file on disk
 * @param f lambda function to operate on each interval
 */
void for_region_in_bed_file(path bed_path, const function<void(const Region& r)>& f){
    ifstream file(bed_path);

    char c;
    string region_name;
    string start_token;
    string stop_token;

    int64_t n_delimiters = 0;
    char delimiter = '\t';

    while (file.get(c)){
//        cerr << ((c == delimiter) ? string(1,c) : "[tab]") << ' ' << region_name << ' ' << start_token << ' ' << stop_token << '\n';

        if (c == delimiter){
            n_delimiters++;
            continue;
        }
        if (c == '\n'){
            int64_t start = stoll(start_token);
            int64_t stop = stoll(stop_token);
            f(Region(region_name, start, stop));

            region_name.clear();
            start_token.clear();
            stop_token.clear();
            n_delimiters = 0;
            continue;
        }

        if (n_delimiters == 0){
            region_name += c;
        }
        else if (n_delimiters == 1){
            start_token += c;
        }
        else if (n_delimiters == 2){
            stop_token += c;
        }
    }
}


}