#include "IntervalGraph.hpp"
#include "Sequence.hpp"
#include "Region.hpp"
#include "bam.hpp"

#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <ranges>
#include <deque>
#include <cmath>
#include <span>

using std::priority_queue;
using std::runtime_error;
using std::to_string;
using std::ofstream;
using std::deque;
using std::swap;
using std::cerr;
using std::cout;
using std::span;
using std::min;
using std::max;


#include "htslib/include/htslib/hts.h"
#include "htslib/include/htslib/sam.h"
#include "htslib/include/htslib/hfile.h"


namespace sv_merge {

string strip_bam_extension(const path& bam){
    string name_prefix = bam.string();

    string suffix_a = ".bam.bai";
    string suffix_b = ".bam";

    if (name_prefix.ends_with(suffix_a)){
        name_prefix = name_prefix.substr(0,name_prefix.size() - suffix_a.size());
    }
    if (name_prefix.ends_with(suffix_b)){
        name_prefix = name_prefix.substr(0,name_prefix.size() - suffix_b.size());
    }

    std::replace(name_prefix.begin(), name_prefix.end(), '.', '_');

    return name_prefix;
}


path get_cache_dir(const path& index_path){
    path cache_dir = index_path;

    // Replace all non alphanumeric chars with '_'
    string name = strip_bam_extension(index_path);
    std::replace_if(name.begin(), name.end(), [&](char c){return !isalnum(c);}, '_');

    cache_dir = std::filesystem::temp_directory_path() / "hapestry" / "bai" / name;

    return cache_dir;
}


hts_idx_t* get_index(const path& bam_path, samFile* bam_file){
    path index_path = bam_path.string() + ".bai";
    auto cache_dir = get_cache_dir(index_path);
    auto cache_path = cache_dir / index_path.filename();
    hts_idx_t *bam_index;

    // first check if the BAI is local
    if (std::filesystem::exists(index_path)){
        cerr << "loading local bai: " << index_path << '\n';

        // Just load the index if it is local
        bam_index = sam_index_load3(bam_file, bam_path.string().c_str(), index_path.c_str(), 0);

        if (bam_index == nullptr) {
            throw runtime_error("ERROR: Cannot open index for bam file: " + bam_path.string() + "\n");
        }
    }
    // Check if the bai is cached already
    else if (std::filesystem::exists(cache_path)) {
        cerr << "cache found: " << cache_path << '\n';

        bam_index = sam_index_load3(bam_file, bam_path.string().c_str(), cache_path.c_str(), 0);

        if (bam_index == nullptr) {
            throw runtime_error("ERROR: Cannot open index for bam file: " + bam_path.string() + "\n");
        }
    }
    else{
        cerr << "cache not found, loading remote bai: " << index_path << '\n';

        std::filesystem::create_directories(cache_dir);

        // Load WITHOUT saving
        bam_index = sam_index_load3(bam_file, bam_path.string().c_str(), index_path.c_str(), 0);

        if (bam_index == nullptr) {
            throw runtime_error("ERROR: Cannot open index for bam file: " + bam_path.string() + "\n");
        }

        /// Save an index to a specific file
        /** @param idx    Index to be written
            @param fn     Input BAM/BCF/etc filename
            @param fnidx  Output filename, or NULL to add .bai/.csi/etc to @a fn
            @param fmt    One of the HTS_FMT_* index formats
            @return  0 if successful, or negative if an error occurred.
        */
        int result = hts_idx_save_as(bam_index, bam_path.c_str(), cache_path.c_str(), HTS_FMT_BAI);

        if (result != 0){
            throw runtime_error("ERROR: failed to write index to path: " + cache_path.string());
        }
    }

    return bam_index;
}


void decompress_cigar_bytes(uint32_t bytes, CigarTuple& cigar){
    cigar.code = int8_t(bytes & bam_cigar_mask);
    cigar.length = int32_t(bytes >> bam_cigar_shift);
}


void decompress_cigar_bytes(uint32_t bytes, CigarInterval& cigar){
    cigar.code = int8_t(bytes & bam_cigar_mask);
    cigar.length = int32_t(bytes >> bam_cigar_shift);
}


HtsAlignment::HtsAlignment(bam1_t* a):
    query_sequence(),
    qualities(),
    hts_alignment(a),
    is_decompressed(false),
    reverse(bam_is_rev(a))
{}


int32_t HtsAlignment::get_query_length() const {
    return int32_t(hts_alignment->core.l_qseq);
}


int32_t HtsAlignment::get_ref_start() const {
    return int32_t(hts_alignment->core.pos);
}


int32_t HtsAlignment::get_ref_stop() const {
    return int32_t(bam_endpos(hts_alignment));
}


int32_t HtsAlignment::get_query_start() const {
    if (is_reverse()){
        return int32_t(get_query_length());
    }
    else{
        return 0;
    }
}


void HtsAlignment::get_query_name(string& result) const {
    result = bam_get_qname(hts_alignment);
}


void HtsAlignment::get_query_sequence(string& result, int32_t start, int32_t stop){
    if (not is_decompressed){
        decompress_bam_sequence(hts_alignment, result, start, stop);
    }
    else{
        result = query_sequence.substr(start, stop-start);
    }
}


void HtsAlignment::get_query_sequence(string& result){
    if (not is_decompressed){
        decompress_bam_sequence(hts_alignment, query_sequence);
        is_decompressed = true;
    }
    result = query_sequence;
}


void HtsAlignment::get_qualities(vector<uint8_t>& result){
    uint8_t* q = bam_get_qual(hts_alignment);

    // Inefficient copy of data
    if (is_reverse()){
        // TODO: any way to initialize this directly using reverse iterators, instead of two steps?
        result = vector<uint8_t>(q, q + hts_alignment->core.l_qseq);
        std::reverse(result.begin(), result.end());
    }
    else{
        result = vector<uint8_t>(q, q + hts_alignment->core.l_qseq);
    }
}


void HtsAlignment::get_tag_as_string(const string& tag_name, string& result, bool allow_missing) const{
    result.clear();

    kstring_t s = {0,0, nullptr};
    int error_code = bam_aux_get_str(hts_alignment, tag_name.c_str(), &s);

    if (error_code != 1){
        if (allow_missing){
            // Return early if the tag is missing and the user allows it (text field will be empty)
            return;
        }

        string query_name;
        get_query_name(query_name);
        throw runtime_error("ERROR: could not fetch tag " + tag_name + " from alignment " + query_name + ", error code: " + to_string(error_code));
    }

    result.assign(s.s, s.l);
}


void HtsAlignment::for_each_cigar_interval(bool unclip_coords, const function<void(const CigarInterval&)>& f) {
    auto cigar_bytes = bam_get_cigar(hts_alignment);
    CigarInterval c;

    // Initialize the cigar interval
    c.query_start = 0;
    c.ref_start = get_ref_start();
    c.is_reverse = is_reverse();

    if (c.is_reverse){
        if (unclip_coords) {
            for_each_cigar_tuple([&](const CigarTuple& c2){
                if (c2.code == 5){
                    c.query_start += c2.length;
                }

                c.query_start += is_query_move[c2.code] * c2.length;
            });
        }
        else {
            c.query_start = get_query_length();
        }
    }

    for (uint32_t i=0; i < hts_alignment->core.n_cigar; i++) {
        decompress_cigar_bytes(cigar_bytes[i], c);

        // Optionally advance the query coord for hardclip (H/5) operations if native/unclipped coords are desired
        if (unclip_coords and c.is_hardclip()){
            c.code = 4;
        }

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


bool HtsAlignment::is_not_primary() const {
    return (hts_alignment->core.flag >> 8) & uint16_t(1);
}


bool HtsAlignment::is_primary() const {
    return (not is_not_primary());
}


bool HtsAlignment::is_supplementary() const {
    return ((hts_alignment->core.flag&BAM_FSUPPLEMENTARY) != 0);
}


bool HtsAlignment::is_reverse() const {
    return reverse;
}


bool HtsAlignment::is_unmapped() const {
    return (hts_alignment->core.tid < 0);
}


void decompress_bam_sequence(const bam1_t* alignment, string& sequence){
    auto length = alignment->core.l_qseq;

    // Fetch 4 bit base code from the correct 8 bit integer and convert to a char
    int32_t start = 0;
    int32_t stop = length;

    decompress_bam_sequence(alignment, sequence, start, stop);
}


void decompress_bam_sequence(const bam1_t* alignment, string& sequence, int32_t start, int32_t stop){
    auto compressed_sequence = bam_get_seq(alignment);
    bool is_reverse = bam_is_rev(alignment);

    sequence.clear();
    sequence.reserve(stop-start);

    uint8_t base_code;

    // Fetch 4 bit base code from the correct 8 bit integer and convert to a char
    int32_t increment = 1;

    if (is_reverse){
        std::swap(start,stop);
        start -= 1;
        stop -= 1;
        increment = -1;
    }

    for (int32_t i=start; i!=stop; i+=increment){
        uint32_t index = i/2;

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


void for_read_in_bam(path bam_path, const function<void(Sequence& sequence)>& f) {
    for_alignment_in_bam(bam_path, [&](Alignment& alignment){
        // The goal is to collect all query sequences, so skip any that may be hardclipped (not primary)
        if ((not alignment.is_primary()) or alignment.is_supplementary()) {
            return;
        }

        Sequence s;

        alignment.get_query_name(s.name);
        alignment.get_query_sequence(s.sequence);

        f(s);
    });
}


void for_alignment_in_bam(path bam_path, const function<void(Alignment& alignment)>& f) {
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

    bam_index = get_index(bam_path, bam_file);

    // bam header
    if ((bam_header = sam_hdr_read(bam_file)) == nullptr) {
        throw runtime_error("ERROR: Cannot open header for bam file: " + bam_path.string() + "\n");
    }

    HtsAlignment a(alignment);

    int result;

    while (true) {
        result = sam_read1(bam_file, bam_header, alignment);

        if (result == -1){
            break;
        }
        else if (result < -1){
            throw runtime_error("ERROR: sam read failed with error code: " + to_string(result));
        }

        if (alignment->core.tid < 0) {
            continue;
        }

        a = HtsAlignment(alignment);

        f(a);
    }

    bam_destroy1(alignment);
    bam_hdr_destroy(bam_header);
    hts_close(bam_file);
    hts_idx_destroy(bam_index);
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

    bam_index = get_index(bam_path, bam_file);

    // bam header
    if ((bam_header = sam_hdr_read(bam_file)) == nullptr) {
        throw runtime_error("ERROR: Cannot open header for bam file: " + bam_path.string() + "\n");
    }

    hts_itr_t *itr = sam_itr_querys(bam_index, bam_header, region.c_str());

    HtsAlignment a(alignment);

    int result;
    while (true) {
        result = sam_itr_next(bam_file, itr, alignment);

        if (result == -1){
            break;
        }
        else if (result < -1){
            throw runtime_error("ERROR: sam read failed with error code: " + to_string(result));
        }

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


/**
 *
 * @param bam_path : path to local (on filesystem) or remote (GS URI) BAM to be iterated
 * @param region : the region which envelopes all the subregions
 * @param subregions : subregions that the user is interested in fetching
 * @param f : returns the SAM alignment, and also a vector of Regions which the alignment overlaps, corresponding to
 *              whichever regions the user provided, as a read may span multiple
 */
void for_alignment_in_bam_subregions(
        path bam_path,
        string region,
        const span<const Region>& subregions,
        const function<void(Alignment& alignment, span<const Region>& overlapping_regions)>& f
        ){

    // Don't want to rely on user compliance with documentation to guarantee that these are sorted
    // TODO: fix the requirement for contiguous region!!
    Region r_prev = subregions[0];
    for (const auto& r: subregions){
        if (r_prev.start > r.start){
            throw runtime_error("ERROR: subregions must be sorted by start position");
        }
        if (r_prev.name != r.name){
            throw runtime_error("ERROR: non-contiguous region, names do not match: " + r.name + " != " + r_prev.name);
        }
    }

    samFile *bam_file;
    bam_hdr_t *bam_header;
    hts_idx_t *bam_index;
    bam1_t *alignment;

    bam_file = nullptr;
    bam_index = nullptr;
    alignment = bam_init1();

    if ((bam_file = hts_open(bam_path.c_str(), "r")) == nullptr) {
        throw runtime_error("ERROR: Cannot open bam file: " + bam_path.string());
    }

    bam_index = get_index(bam_path, bam_file);

    // bam header
    if ((bam_header = sam_hdr_read(bam_file)) == nullptr) {
        throw runtime_error("ERROR: Cannot open header for bam file: " + bam_path.string() + "\n");
    }

    hts_itr_t *itr = sam_itr_querys(bam_index, bam_header, region.c_str());

    HtsAlignment a(alignment);

    size_t i_start = 0;

    // Iterate a SORTED bam
    int result;
    while (true) {
        result = sam_itr_next(bam_file, itr, alignment);

        if (result == -1){
            break;
        }
        else if (result < -1){
            throw runtime_error("ERROR: sam read failed with error code: " + to_string(result));
        }

        if (alignment->core.tid < 0) {
            continue;
        }

        if (i_start >= subregions.size()){
            continue;
        }

        a = HtsAlignment(alignment);

        auto alignment_start = a.get_ref_start();
        auto alignment_stop = a.get_ref_stop();

        // Set the iterator beyond regions that have been passed by the alignments already
        while (alignment_start > subregions[i_start].stop){
            i_start++;

            if (i_start >= subregions.size()){
                break;
            }
        }

        auto i_stop = i_start;

        // Iterate all the subregions that intersect the alignment
        while (i_stop < subregions.size() and min(subregions[i_stop].stop, alignment_stop) - max(subregions[i_stop].start, alignment_start) > 0){
            i_stop++;
        }

        auto s = span(subregions).subspan(i_start, i_stop - i_start);

        f(a, s);
    }

    bam_destroy1(alignment);
    hts_itr_destroy(itr);
    bam_hdr_destroy(bam_header);
    hts_close(bam_file);
    hts_idx_destroy(bam_index);
}


}
