#include "fasta.hpp"

#include <iostream>

using std::ifstream;
using std::cerr;


namespace sv_merge {

/**
 * Iterate over the regions of a bed file, ignoring annotation fields. Must end in new line if last line is to be
 * considered. Must use actual tab characters (not spaces)
 * @param bed_path path to bed file on disk
 * @param f lambda function to operate on each interval
 */
void for_sequence_in_fasta_file(path fasta_path, const function<void(const Sequence& s)>& f){
    ifstream file(fasta_path);

    char c;
    Sequence s;

    char header_char = '>';
    bool is_sequence = false;


    while (file.get(c)){
//        cerr << c << ' ' << region_name << ' ' << start_token << ' ' << stop_token << '\n';

        if (c == header_char){
            if (not s.name.empty()){
                f(s);
                s.name.clear();
                s.sequence.clear();
                is_sequence = false;
            }
            continue;
        }

        if (c == '\n'){
            is_sequence = true;
            continue;
        }

        if (is_sequence){
            s.sequence += c;
        }
        else{
            s.name += c;
        }
    }
}

} // sv_merge