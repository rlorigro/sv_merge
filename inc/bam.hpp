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


static const array <string, 2> bases = {"=ACMGRSVTWYHKDBN", "=TGKCYSBAWRDKHVN"};

/**
 * ****************************
 * *** CIGAR related macros ***
 * ****************************
 *
 * #define BAM_CMATCH      0
 * #define BAM_CINS        1
 * #define BAM_CDEL        2
 * #define BAM_CREF_SKIP   3
 * #define BAM_CSOFT_CLIP  4
 * #define BAM_CHARD_CLIP  5
 * #define BAM_CPAD        6
 * #define BAM_CEQUAL      7
 * #define BAM_CDIFF       8
 * #define BAM_CBACK       9
 *
**/


static const array <char,9> cigar_code_to_char = {
        'M', // 0
        'I', // 1
        'D', // 2
        'N', // 3
        'S', // 4
        'H', // 5
        'P', // 6
        '=', // 7
        'X', // 8
};


static const array <char,9> cigar_code_to_format_char = {
        '|', // M
        ' ', // I
        ' ', // D
        ' ', // N
        ' ', // S
        ' ', // H
        ' ', // P
        '|', // =
        '*', // X
};


static const array <bool,9> is_query_move = {
        1, // M
        1, // I
        0, // D
        0, // N
        1, // S
        0, // H
        0, // P
        1, // =
        1  // X
};


static const array <bool,9> is_ref_move = {
        1, // M
        0, // I
        1, // D
        1, // N
        0, // S
        0, // H
        0, // P
        1, // =
        1  // X
};



// Each nt in the sequence is stored within a uint8 with 8 bits, XXXXYYYY, where XXXX is nt1 and YYYY is nt2
static const uint8_t bam_sequence_shift = 4;
static const uint8_t bam_sequence_mask = 15;       // 0b1111
static const uint8_t bam_cigar_shift = 4;
static const uint8_t bam_cigar_mask = 15;          // 0b1111


class CigarTuple{
public:
    int64_t length = -1;
    int8_t code = -1;
};


class CigarInterval{
public:
    int64_t length = -1;
    int64_t ref_start = -1;
    int64_t ref_stop = -1;
    int64_t query_start = -1;
    int64_t query_stop = -1;
    int8_t code = -1;
    bool is_reverse;

    pair<int64_t,int64_t> get_forward_ref_interval() const;
    pair<int64_t,int64_t> get_forward_query_interval() const;
    void set_ref_interval_forward();
    void set_query_interval_forward();
    void set_ref_interval_reverse();
    void set_query_interval_reverse();
    bool is_softclip() const;
    bool is_hardclip() const;
    bool is_clip() const;
};


char get_complement(char c);

void reverse_complement(string& seq);

void decompress_bam_sequence(const bam1_t* alignment, string& sequence);

void decompress_cigar_bytes(uint32_t bytes, CigarTuple& cigar);

void for_cigar_tuple_in_alignment(const bam1_t *alignment, const function<void(CigarTuple& cigar)>& f);

void for_cigar_interval_in_alignment(const bam1_t *alignment, const function<void(CigarInterval& cigar)>& f);

void for_cigar_interval_in_alignment(
        const bam1_t *alignment,
        vector<interval_t>& ref_intervals,
        vector<interval_t>& query_intervals,
        const function<void(const CigarInterval& intersection, const CigarInterval& original, const interval_t& interval)>& f_ref,
        const function<void(const CigarInterval& intersection, const CigarInterval& original, const interval_t& interval)>& f_query
);

void for_alignment_in_bam_region(path bam_path, string region, const function<void(const bam1_t *alignment)>& f);

void for_read_in_bam_region(path bam_path, string region, const function<void(Sequence& sequence)>& f);



}
