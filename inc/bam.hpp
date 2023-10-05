#include <functional>
#include <cstdlib>
#include <utility>
#include <string>
#include <array>

using std::function;
using std::string;
using std::array;
using std::pair;
using std::pair;

#include "htslib/include/htslib/sam.h"
#include "Filesystem.hpp"
#include "Sequence.hpp"

using ghc::filesystem::path;


namespace sv_merge {


const array <string, 2> bases = {"=ACMGRSVTWYHKDBN", "=TGKCYSBAWRDKHVN"};


// Each nt in the sequence is stored within a uint8 with 8 bits, XXXXYYYY, where XXXX is nt1 and YYYY is nt2
static const uint8_t bam_sequence_shift = 4;
static const uint8_t bam_sequence_mask = 15;       // 0b1111
static const uint8_t bam_cigar_shift = 4;
static const uint8_t bam_cigar_mask = 15;          // 0b1111


class CigarTuple{
public:
    int8_t code = -1;
    int64_t length = -1;
};

void decompress_bam_sequence(uint8_t* compressed_sequence, int64_t length, bool is_reverse, string& sequence);

void decompress_cigar_bytes(uint32_t bytes, CigarTuple& cigar_tuple);

void iterate_cigar_tuples(bam1_t *alignment, function<void(CigarTuple& cigar_tuple)>& f);

void for_read_in_bam_region(path bam_path, string region, const function<void(Sequence& sequence)>& f);


}
