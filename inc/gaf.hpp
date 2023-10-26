#pragma once

#include "Alignment.hpp"

#include <functional>
#include <stdexcept>
#include <cstdlib>
#include <utility>
#include <fstream>
#include <string>
#include <array>

using std::runtime_error;
using std::function;
using std::ifstream;
using std::string;
using std::array;
using std::pair;
using std::pair;

namespace sv_merge{


class GafAlignment: public Alignment{
    vector <pair <string,bool> > path;
    string query_name;
    string ref_name;
    int64_t query_length;
    int64_t query_start;
    int64_t query_stop;
    int64_t path_length;
    int64_t path_start;
    int64_t path_stop;
    int64_t n_match;
    int64_t alignment_length;
    int64_t map_quality;
    bool reversal;

public:
    void set_path(const string& path);
    string get_path_string() const;
    void set_query_name(const string& query_name);
    string get_query_name() const;
    void set_ref_name(const string& ref_name);
    string get_ref_name() const;
    void set_query_length(int64_t query_length);
    int64_t get_query_length() const;
    void set_query_start(int64_t query_start);
    int64_t get_query_start() const;
    void set_query_stop(int64_t query_stop);
    int64_t get_query_stop() const;
    void set_path_length(int64_t path_length);
    int64_t get_path_length() const;
    void set_path_start(int64_t path_start);
    int64_t get_path_start() const;
    void set_path_stop(int64_t path_stop);
    int64_t get_path_stop() const;
    void set_reversal(bool reversal);
    bool get_reversal() const;
    void set_n_match(int64_t n);
    int64_t get_n_match() const;
    void set_alignment_length(int64_t length);
    int64_t get_alignment_length() const;
    void set_map_quality(int64_t q);
    int64_t get_map_quality() const;

    bool parse_path_reversal_token(char c);

    void for_each_cigar_interval(const function<void(const CigarInterval& cigar)>& f);
    void for_each_cigar_tuple(const function<void(const CigarTuple& cigar)>& f);
    void get_query_sequence(string& result);
    void get_query_name(string& result) const;
    int64_t get_ref_start() const;
    bool is_unmapped() const;
    bool is_reverse() const;
};


bool parse_reversal_token(const string& token);

void for_alignment_in_gaf(const path& gaf_path, const function<void(GafAlignment& alignment)>& f);

void for_alignment_in_gaf(const path& gaf_path, const function<void(Alignment& alignment)>& f);


}
