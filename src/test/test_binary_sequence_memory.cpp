#include "BinarySequence.hpp"
#include "Sequence.hpp"
#include <unordered_set>
#include <filesystem>
#include <fstream>
#include <list>
#include <map>

using sv_merge::BinarySequence;
using sv_merge::reverse_complement;
using std::filesystem::path;
using std::unordered_set;
using std::ifstream;
using std::list;
using std::map;
using std::runtime_error;



int main(){
    {
        path absolute_fasta_path = "/home/ryan/data/human/reference/chm13v2.0_chr1.fa";

        ifstream file(absolute_fasta_path);

        string line;
        string seq;

        while (getline(file, line)) {
            if (line[0] != '>'){
                seq += line;
            }
        }

        cerr << seq.size() << '\n';

        BinarySequence<uint64_t> bs;

        cerr << "fixing ambiguous bases" << '\n';
        for (size_t i=0; i<seq.size(); i++){
            if (seq[i] != 'A' or seq[i] != 'T' or seq[i] != 'C' or seq[i] != 'G'){
                seq[i] = 'A';
            }
        }

        cerr << "encoding as 2bit" << '\n';

        bs.encode(seq);

        cerr << bs.sequence.size() << '\n';
        cerr << bs.get_byte_length() << '\n';
        cerr << bs.size() << '\n';
        cerr << "Done" << '\n';

        bs = {};

        seq = {};

        cerr << '\n';

    }

    return 0;
}

