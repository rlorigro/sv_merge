#include "fasta.hpp"

#include <stdexcept>
#include <iostream>
#include <fstream>

using std::runtime_error;
using std::ifstream;
using std::cerr;


namespace sv_merge {

/**
 * Iterate over the sequences of a Fasta file, ignoring header annotation fields. Must end in new line if last line is
 * to be considered.
 * @param fasta_path path to fasta file on disk
 * @param f lambda function to operate on each sequence
 */
void for_sequence_in_fasta_file(path fasta_path, const function<void(const Sequence& s)>& f){
    if (not (fasta_path.extension() == ".fasta"
            or fasta_path.extension() != ".fa"
            or fasta_path.extension() != ".fna")){
        throw runtime_error("ERROR: file does not have compatible fasta extension: " + fasta_path.string());
    }

    ifstream file(fasta_path);

    if (not (file.is_open() and file.good())){
        throw runtime_error("ERROR: could not read file: " + fasta_path.string());
    }

    char c;
    Sequence s;

    char header_char = '>';
    bool is_sequence = false;


    while (file.get(c)){
//        cerr << c << " - " << s.name << ' ' << s.sequence << '\n';

        if (c == header_char){
            if (is_sequence){
                f(s);
                s.name.clear();
                s.sequence.clear();
                is_sequence = false;
            }
            continue;
        }

        // Skip whatever tags come after the name in the fasta header
        if (c != '\n' and isspace(c)){
            while (file.peek() != '\n' or file.eof()){
                file.get(c);
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

    if (is_sequence){
        f(s);
    }
}

} // sv_merge