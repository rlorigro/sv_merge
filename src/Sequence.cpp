#include "Sequence.hpp"

#include <stdexcept>

using std::runtime_error;


namespace sv_merge {


Sequence::Sequence(const string& name, const string& sequence):
        name(name),
        sequence(sequence)
{}


char get_complement(char c){
    if (isupper(c)) {
        switch (c){
            case 'A': return 'T';
            case 'C': return 'G';
            case 'G': return 'C';
            case 'T': return 'A';
            case 'N': return 'N';
            default:
                throw runtime_error("ERROR: cannot complement base: " + string(1,c));
        }
    }

    else {
        switch (c){
            case 'a': return 't';
            case 'c': return 'g';
            case 'g': return 'c';
            case 't': return 'a';
            case 'n': return 'n';
            default:
                throw runtime_error("ERROR: cannot complement base: " + string(1,c));
        }
    }
}


void reverse_complement(string& seq){
    string s;

    for (auto iter = seq.rbegin(); iter != seq.rend(); iter++){
        s += get_complement(*iter);
    }

    seq = std::move(s);
}

} // sv_merge
