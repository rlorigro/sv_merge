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
#include "Alignment.hpp"
#include "Sequence.hpp"

using ghc::filesystem::path;


namespace sv_merge {

/// Htslib encoding for bytes in a sequence
static const array <string, 2> bases = {"=ACMGRSVTWYHKDBN", "=TGKCYSBAWRDKHVN"};

// Each nt in the sequence is stored within a uint8 with 8 bits, XXXXYYYY, where XXXX is nt1 and YYYY is nt2
static const uint8_t bam_sequence_shift = 4;
static const uint8_t bam_sequence_mask = 15;       // 0b1111
static const uint8_t bam_cigar_shift = 4;
static const uint8_t bam_cigar_mask = 15;          // 0b1111


class HtsAlignment: public Alignment{
private:
    string query_sequence;
    bam1_t* hts_alignment;
    bool is_decompressed = false;
    bool reverse;

public:
    HtsAlignment(bam1_t* a);

    /// Iterating
    void for_each_cigar_interval(const function<void(const CigarInterval&)>& f) override;
    void for_each_cigar_tuple(const function<void(const CigarTuple&)>& f) override;

    /// Accessing
    // For the htslib bam implementation of this, we can lazily evaluate, only decompress and store the sequence as
    // needed. Unfortunately this makes the getter non-const...
    void get_query_sequence(string& result) override;
    void get_query_name(string& result) const override;
    int64_t get_query_length() const override;
    int64_t get_ref_start() const override;
    int64_t get_query_start() const override;
    bool is_unmapped() const override;
    bool is_reverse() const override;
};


char get_complement(char c);

void reverse_complement(string& seq);

void decompress_bam_sequence(const bam1_t* alignment, string& sequence);

void decompress_cigar_bytes(uint32_t bytes, CigarTuple& cigar);

void for_alignment_in_bam_region(path bam_path, string region, const function<void(Alignment& alignment)>& f);

void for_read_in_bam_region(path bam_path, string region, const function<void(Sequence& sequence)>& f);



}
