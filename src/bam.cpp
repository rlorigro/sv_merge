#include "bam.hpp"

#include <stdexcept>
#include <iostream>

using std::runtime_error;
using std::to_string;
using std::cerr;
using std::cout;


#include "htslib/include/htslib/hts.h"


namespace sv_merge {

void decompress_cigar_bytes(uint32_t bytes, CigarTuple& cigar_tuple){
    cigar_tuple.code = int8_t(bytes & bam_cigar_mask);
    cigar_tuple.length = int64_t(bytes >> bam_cigar_shift);
}


void iterate_cigar_tuples(bam1_t *alignment, function<void(CigarTuple& cigar_tuple)>& f){
    auto cigar_bytes = bam_get_cigar(alignment);
    for (int32_t i=0; i<alignment->core.n_cigar; i++){
        CigarTuple c;
        decompress_cigar_bytes(cigar_bytes[i], c);
        f(c);
    }
}


void decompress_bam_sequence(uint8_t* compressed_sequence, int64_t length, bool is_reverse, string& sequence){
    ///
    /// Convert the compressed representation of an aligned sequence into a string.
    /// Does NOT reverse complement the sequence
    ///

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

        auto length = alignment->core.l_qseq;
        auto compressed_sequence = bam_get_seq(alignment);

        Sequence s;
        s.name = bam_get_qname(alignment);
        bool is_reverse = bam_is_rev(alignment);

        decompress_bam_sequence(compressed_sequence, length, is_reverse, s.sequence);

        f(s);
    }

    bam_destroy1(alignment);
    hts_itr_destroy(itr);
    bam_hdr_destroy(bam_header);
    hts_close(bam_file);
    hts_idx_destroy(bam_index);
}


}
