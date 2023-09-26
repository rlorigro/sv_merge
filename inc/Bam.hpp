#include <functional>
#include <cstdlib>
#include <string>
#include <array>

using std::function;
using std::string;
using std::array;

#include "Filesystem.hpp"

using ghc::filesystem::path;


namespace sv_merge {


class Sequence{
public:
    string name;
    string sequence;
};

class Region {
public:
    string name;
    int64_t start;
    int64_t stop;

    Region(string& region_string);
};


const array<char,16> bases = {
        '=', // 0
        'A', // 1
        'C', // 2
        'M', // 3
        'G', // 4
        'R', // 5
        'S', // 6
        'V', // 7
        'T', // 8
        'W', // 9
        'Y', // 10
        'H', // 11
        'K', // 12
        'D', // 13
        'B', // 14
        'N'  // 15
};


// Each nt in the sequence is stored within a uint8 with 8 bits, XXXXYYYY, where XXXX is nt1 and YYYY is nt2
static const uint8_t bam_sequence_shift = 4;
static const uint8_t bam_sequence_mask = 15;       // 0b1111
static const uint8_t bam_cigar_shift = 4;
static const uint8_t bam_cigar_mask = 15;          // 0b1111

void decompress_bam_sequence(uint8_t* compressed_sequence, int64_t length, string& sequence);

void for_read_in_bam_region(path bam_path, string region, const function<void(Sequence& sequence)>& f);


}
