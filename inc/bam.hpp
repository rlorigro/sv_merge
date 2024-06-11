#pragma once

#include <filesystem>
#include <functional>
#include <cstdlib>
#include <utility>
#include <string>
#include <array>
#include <span>

using std::filesystem::path;
using std::function;
using std::string;
using std::array;
using std::pair;
using std::span;

#include "htslib/include/htslib/sam.h"
#include "Alignment.hpp"
#include "Sequence.hpp"
#include "Region.hpp"



namespace sv_merge {

/// Htslib encoding for bytes in a sequence
static const array <string, 2> bases = {"=ACMGRSVTWYHKDBN", "=TGKCYSBAWRDKHVN"};

// Each nt in the sequence is stored within a uint8 with 8 bits, XXXXYYYY, where XXXX is nt1 and YYYY is nt2
static const uint8_t bam_sequence_shift = 4;
static const uint8_t bam_sequence_mask = 15;       // 0b1111
static const uint8_t bam_cigar_shift = 4;
static const uint8_t bam_cigar_mask = 15;          // 0b1111


/**
 * Wrapper for htslib alignment struct. This is derived from Alignment which is intended to be used interchangeably with
 * other implementations of an "alignment" such as GAF or PAF
 */
class HtsAlignment: public Alignment{
private:
    string query_sequence;
    vector<uint8_t> qualities;
    bam1_t* hts_alignment;
    bool is_decompressed = false;
    bool reverse;

public:
    HtsAlignment(bam1_t* a);

    /// Iterating
    void for_each_cigar_interval(bool unclip_coords, const function<void(const CigarInterval&)>& f) override;
    void for_each_cigar_tuple(const function<void(const CigarTuple&)>& f) override;

    /// Accessing

    // For the htslib bam implementation of this, we can lazily evaluate, only decompress and store the sequence as
    // needed. Unfortunately this makes the getter non-const...
    void get_query_sequence(string& result, int32_t start, int32_t stop) override;
    void get_query_sequence(string& result) override;
    void get_qualities(vector<uint8_t>& result) override;
    void get_query_name(string& result) const override;
    void get_tag_as_string(const string& name, string& result, bool allow_missing=false) const override;
    [[nodiscard]] int32_t get_query_length() const override;
    [[nodiscard]] int32_t get_ref_start() const override;
    [[nodiscard]] int32_t get_ref_stop() const override;
    [[nodiscard]] int32_t get_query_start() const override;
    [[nodiscard]] bool is_unmapped() const override;
    [[nodiscard]] bool is_reverse() const override;
    [[nodiscard]] bool is_primary() const override;
    [[nodiscard]] bool is_supplementary() const override;
    [[nodiscard]] bool is_not_primary() const;
};


/**
 * Inefficient complementation, good luck to the optimizer TODO: switch to array based or char ordinal offset math
 * @param c - the character to be complemented, must be ACGTN or acgtn
 * @return The complement of c
 */
char get_complement(char c);

void reverse_complement(string& seq);

/**
 * Convert the compressed representation of an aligned sequence into a string.
 * WARNING: reverse complements the sequence to be F orientation
 * @param alignment - pointer to htslib struct
 * @param sequence - string to be filled with resulting sequence
 */
void decompress_bam_sequence(const bam1_t* alignment, string& sequence);


/**
 * Convert the compressed representation of an aligned sequence into a string.
 * WARNING: reverse complements the sequence to be F orientation
 * @param alignment - pointer to htslib struct
 * @param sequence - string to be filled with resulting sequence
 * @param start - start index, in the F orientation of the sequence
 * @param stop - stop index, in the F orientation of the sequence
 */
void decompress_bam_sequence(const bam1_t* alignment, string& sequence, int32_t start, int32_t stop);

void decompress_cigar_bytes(uint32_t bytes, CigarTuple& cigar);

void for_alignment_in_bam(path bam_path, const function<void(Alignment& alignment)>& f);

void for_alignment_in_bam_region(path bam_path, string region, const function<void(Alignment& alignment)>& f);

void for_read_in_bam(path bam_path, const function<void(Sequence& sequence)>& f);

void for_read_in_bam_region(path bam_path, string region, const function<void(Sequence& sequence)>& f);

void for_alignment_in_bam_subregions(
        path bam_path,
        string region,
        const span<const Region>& subregions,
        const function<void(Alignment& alignment, span<const Region>& overlapping_regions)>& f
);


}
