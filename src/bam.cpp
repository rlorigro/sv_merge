#include "bam.hpp"

#include <ranges>
#include <stdexcept>
#include <iostream>

using std::runtime_error;
using std::to_string;
using std::cerr;
using std::cout;


#include "htslib/include/htslib/hts.h"


namespace sv_merge {

/**
 * Inefficient complementation, good luck to the optimizer
 * @param c - the character to be complemented, must be ACGTN or acgtn
 * @return The complement of c
 */
char get_complement(char c){
    if (isupper(c)) {
        if (c == 'A') {
            return 'T';
        }
        else if (c == 'C') {
            return 'G';
        }
        else if (c == 'G') {
            return 'C';
        }
        else if (c == 'T') {
            return 'A';
        }
        else if (c == 'N') {
            return 'N';
        }
        else {
            throw runtime_error("ERROR: cannot complement base: " + string(1,c));
        }
    }

    else {
        if (c == 'a') {
            return 't';
        }
        else if (c == 'c') {
            return 'g';
        }
        else if (c == 'g') {
            return 'c';
        }
        else if (c == 't') {
            return 'a';
        }
        else if (c == 'n') {
            return 'n';
        }
        else {
            throw runtime_error("ERROR: cannot complement base: " + string(1,c));
        }
    }
}


void reverse_complement(string& seq){
    string s;

    for (char & iter : std::ranges::reverse_view(seq)){
        s += get_complement(iter);
    }

    seq = s;
}


pair<int64_t,int64_t> CigarInterval::get_forward_ref_interval() const{
    if (ref_start > ref_stop){
        return {ref_stop,ref_start};
    }
    else {
        return {ref_start,ref_stop};
    }
}


pair<int64_t,int64_t> CigarInterval::get_forward_query_interval() const{
    if (query_start > query_stop){
        return {query_stop,query_start};
    }
    else {
        return {query_start,query_stop};
    }
}


void decompress_cigar_bytes(uint32_t bytes, CigarTuple& cigar){
    cigar.code = int8_t(bytes & bam_cigar_mask);
    cigar.length = int64_t(bytes >> bam_cigar_shift);
}


void decompress_cigar_bytes(uint32_t bytes, CigarInterval& cigar){
    cigar.code = int8_t(bytes & bam_cigar_mask);
    cigar.length = int64_t(bytes >> bam_cigar_shift);
}


void for_cigar_tuple_in_alignment(const bam1_t *alignment, const function<void(CigarTuple& cigar)>& f){
    auto cigar_bytes = bam_get_cigar(alignment);
    for (uint32_t i=0; i<alignment->core.n_cigar; i++){
        CigarTuple c;
        decompress_cigar_bytes(cigar_bytes[i], c);
        f(c);
    }
}


void for_cigar_interval_in_alignment(const bam1_t *alignment, const function<void(CigarInterval& cigar)>& f){
    auto cigar_bytes = bam_get_cigar(alignment);

    CigarInterval c;
    c.query_start = 0;
    c.ref_start = alignment->core.pos;
    c.is_reverse = bam_is_rev(alignment);

    if (c.is_reverse){
        c.query_start = alignment->core.l_qseq;
    }

    for (uint32_t i=0; i<alignment->core.n_cigar; i++){
        decompress_cigar_bytes(cigar_bytes[i], c);

        if (c.is_reverse) {
            c.query_stop = c.query_start - is_query_move[c.code]*c.length;
        }
        else{
            c.query_stop = c.query_start + is_query_move[c.code]*c.length;
        }

        c.ref_stop = c.ref_start + is_ref_move[c.code]*c.length;

        f(c);

        c.ref_start = c.ref_stop;
        c.query_start = c.query_stop;
    }
}


void decompress_bam_sequence(const bam1_t* alignment, string& sequence){
    ///
    /// Convert the compressed representation of an aligned sequence into a string.
    /// Does NOT reverse complement the sequence
    ///
    auto length = alignment->core.l_qseq;
    auto compressed_sequence = bam_get_seq(alignment);
    bool is_reverse = bam_is_rev(alignment);

    sequence.clear();

    uint8_t base_code;
    string base;

    // Fetch 4 bit base code from the correct 8 bit integer and convert to a char
    int64_t start = 0;
    int64_t stop = length;
    int64_t increment = 1;

    if (is_reverse){
        start = length - 1;
        stop = -1;
        increment = -1;
    }

    for (int64_t i=start; i!=stop; i+=increment){
        uint64_t index = i/2;

        if (i%2 == 0){
            // Perform bit SHIFT and decode using the standard or complemented base map
            base_code = compressed_sequence[index] >> bam_sequence_shift;
            sequence += bases[is_reverse][base_code];
        }
        else {
            // Perform bit MASK and decode using the standard or complemented base map
            base_code = compressed_sequence[index] & bam_sequence_mask;
            sequence += bases[is_reverse][base_code];
        }
    }
}


void for_read_in_bam_region(path bam_path, string region, const function<void(Sequence& sequence)>& f) {
    for_alignment_in_bam_region(bam_path, region, [&](const bam1_t *alignment){
        if (alignment->core.tid < 0) {
            return;
        }

        Sequence s;
        s.name = bam_get_qname(alignment);

        decompress_bam_sequence(alignment, s.sequence);

        f(s);
    });
}


void for_alignment_in_bam_region(path bam_path, string region, const function<void(const bam1_t *alignment)>& f) {
    samFile *bam_file;
    bam_hdr_t *bam_header;
    hts_idx_t *bam_index;
    bam1_t *alignment;

    bam_file = nullptr;
    bam_index = nullptr;
    alignment = bam_init1();

    if ((bam_file = hts_open(bam_path.string().c_str(), "r")) == nullptr) {
        throw runtime_error("ERROR: Cannot open bam file: " + bam_path.string());
    }

    // bam index
    if ((bam_index = sam_index_load(bam_file, bam_path.string().c_str())) == nullptr) {
        throw runtime_error("ERROR: Cannot open index for bam file: " + bam_path.string() + "\n");
    }

    // bam header
    if ((bam_header = sam_hdr_read(bam_file)) == nullptr) {
        throw runtime_error("ERROR: Cannot open header for bam file: " + bam_path.string() + "\n");
    }

    hts_itr_t *itr = sam_itr_querys(bam_index, bam_header, region.c_str());

    while (sam_itr_next(bam_file, itr, alignment) >= 0) {
        if (alignment->core.tid < 0) {
            continue;
        }

        f(alignment);
    }

    bam_destroy1(alignment);
    hts_itr_destroy(itr);
    bam_hdr_destroy(bam_header);
    hts_close(bam_file);
    hts_idx_destroy(bam_index);
}


}
