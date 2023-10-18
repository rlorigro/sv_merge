#include "IntervalGraph.hpp"
#include "bam.hpp"

#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <ranges>
#include <deque>
#include <cmath>

using std::runtime_error;
using std::to_string;
using std::deque;
using std::swap;
using std::cerr;
using std::cout;
using std::min;
using std::max;


#include "htslib/include/htslib/hts.h"


namespace sv_merge {

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

    for (char & iter : std::ranges::reverse_view(seq)){
        s += get_complement(iter);
    }

    seq = s;
}



void decompress_cigar_bytes(uint32_t bytes, CigarTuple& cigar){
    cigar.code = int8_t(bytes & bam_cigar_mask);
    cigar.length = int64_t(bytes >> bam_cigar_shift);
}


void decompress_cigar_bytes(uint32_t bytes, CigarInterval& cigar){
    cigar.code = int8_t(bytes & bam_cigar_mask);
    cigar.length = int64_t(bytes >> bam_cigar_shift);
}


HtsAlignment::HtsAlignment(bam1_t* a):
    query_sequence(""),
    hts_alignment(a),
    is_decompressed(false),
    reverse(bam_is_rev(a))
{}


int64_t HtsAlignment::get_query_length() const {
    return int64_t(hts_alignment->core.l_qseq);
}


int64_t HtsAlignment::get_ref_start() const {
    return int64_t(hts_alignment->core.pos);
}


int64_t HtsAlignment::get_query_start() const {
    if (is_reverse()){
        return int64_t(get_query_length());
    }
    else{
        return 0;
    }
}


void HtsAlignment::get_query_name(string& result) const {
    result = bam_get_qname(hts_alignment);
}


void HtsAlignment::get_query_sequence(string& result){
    if (not is_decompressed){
        decompress_bam_sequence(hts_alignment, query_sequence);
        is_decompressed = true;
    }
    result = query_sequence;
}


void HtsAlignment::for_each_cigar_interval(const function<void(const CigarInterval&)>& f) {
    auto cigar_bytes = bam_get_cigar(hts_alignment);
    CigarInterval c;

    // Initialize the cigar interval
    c.query_start = 0;
    c.ref_start = get_ref_start();
    c.is_reverse = is_reverse();

    if (c.is_reverse){
        c.query_start = get_query_length();
    }

    for (uint32_t i=0; i < hts_alignment->core.n_cigar; i++) {
        decompress_cigar_bytes(cigar_bytes[i], c);

        // Update interval bounds for this cigar interval
        if (c.is_reverse) {
            c.query_stop = c.query_start - is_query_move[c.code]*c.length;
        }
        else{
            c.query_stop = c.query_start + is_query_move[c.code]*c.length;
        }

        c.ref_stop = c.ref_start + is_ref_move[c.code]*c.length;

        // Temporarily flip the start/stop so that it is conventionally interpretable
        c.set_query_interval_forward();

        f(c);

        // Revert to backwards intervals for iteration/update
        if (c.is_reverse){
            c.set_query_interval_reverse();
        }

        c.ref_start = c.ref_stop;
        c.query_start = c.query_stop;
    }
}


void HtsAlignment::for_each_cigar_tuple(const function<void(const CigarTuple&)>& f) {
    auto cigar_bytes = bam_get_cigar(hts_alignment);
    CigarTuple c;

    for (uint32_t i=0; i < hts_alignment->core.n_cigar; i++) {
        decompress_cigar_bytes(cigar_bytes[i], c);

        f(c);
    }
}


bool HtsAlignment::is_reverse() const {
    return reverse;
}


bool HtsAlignment::is_unmapped() const {
    return (hts_alignment->core.tid < 0);
}


void decompress_bam_sequence(const bam1_t* alignment, string& sequence){
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
    for_alignment_in_bam_region(bam_path, region, [&](Alignment& alignment){
        if (alignment.is_unmapped()) {
            return;
        }

        Sequence s;

        alignment.get_query_name(s.name);
        alignment.get_query_sequence(s.sequence);

        f(s);
    });
}


void for_alignment_in_bam_region(path bam_path, string region, const function<void(Alignment& alignment)>& f) {
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

    HtsAlignment a(alignment);

    while (sam_itr_next(bam_file, itr, alignment) >= 0) {
        if (alignment->core.tid < 0) {
            continue;
        }

        a = HtsAlignment(alignment);

        f(a);
    }

    bam_destroy1(alignment);
    hts_itr_destroy(itr);
    bam_hdr_destroy(bam_header);
    hts_close(bam_file);
    hts_idx_destroy(bam_index);
}


}
