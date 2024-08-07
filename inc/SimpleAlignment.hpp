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

#include "Alignment.hpp"
#include "Sequence.hpp"
#include "Region.hpp"



namespace sv_merge {


/**
 * Wrapper for a generic alignment which is any two sequences in F orientation and their corresponding cigar string.
 * This is derived from Alignment which is intended to be used interchangeably with other implementations of an
 * "alignment" such as GAF or PAF.
 *
 * This class cannot represent "reverse" alignments. Strings must be reversed before providing them.
 *
 * WARNING: this class does NOT take ownership of the sequences or cigar string. They must outlive this object to
 * guarantee that the pointers remain valid.
 */
class SimpleAlignment: public Alignment{
private:
    const string& ref_sequence;
    const string& query_sequence;
    const string& cigar;

public:
    SimpleAlignment(const string& ref_sequence, const string& query_sequence, const string& cigar);

    /// Iterating
    void for_each_cigar_interval(bool unclip_coords, const function<void(const CigarInterval&)>& f) override;
    void for_each_cigar_tuple(const function<void(const CigarTuple&)>& f) override;

    /// Accessing
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
};


}
