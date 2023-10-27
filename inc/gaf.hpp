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
    vector <CigarTuple> cigar;
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
    bool primary;

public:
    /// Modifying
    void set_path(const vector<pair<string,bool> >& p);
    void set_path(const string& path);
    void set_query_name(const string& query_name);
    void set_ref_name(const string& ref_name);
    void set_query_length(int64_t query_length);
    void set_query_start(int64_t query_start);
    void set_query_stop(int64_t query_stop);
    void set_path_length(int64_t path_length);
    void set_path_start(int64_t path_start);
    void set_path_stop(int64_t path_stop);
    void set_reversal(bool reversal);
    void set_n_match(int64_t n);
    void set_alignment_length(int64_t length);
    void set_map_quality(int64_t q);
    void set_is_primary(bool p);
    void add_tag(const string& tag);

    /// Accessing
    string get_path_string() const;
    string get_query_name() const;
    string get_ref_name() const;
    int64_t get_query_length() const override;
    int64_t get_query_start() const override;
    int64_t get_query_stop() const;
    int64_t get_path_length() const;
    int64_t get_path_start() const;
    int64_t get_path_stop() const;
    int64_t get_n_match() const;
    int64_t get_alignment_length() const;
    int64_t get_map_quality() const;
    bool is_primary() const;
    bool is_reverse() const override;

    /// Helper
    bool parse_path_reversal_token(char c) const;
    void load_cigar(const string& cigar_string);

    /// Iterators and Alignment API fulfillment
    void for_each_cigar_interval(const function<void(const CigarInterval& cigar)>& f) override;
    void for_each_cigar_tuple(const function<void(const CigarTuple& cigar)>& f) override;
    void get_query_sequence(string& result) override;
    void get_query_name(string& result) const override;
    int64_t get_ref_start() const override;
    bool is_unmapped() const override;
};


bool parse_reversal_token(const string& token);

void for_alignment_in_gaf(const path& gaf_path, const function<void(GafAlignment& alignment)>& f);

void for_alignment_in_gaf(const path& gaf_path, const function<void(Alignment& alignment)>& f);


}
