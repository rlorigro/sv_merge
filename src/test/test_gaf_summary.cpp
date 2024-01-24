#include "gaf.hpp"

using namespace sv_merge;

#include <iostream>
#include <stdexcept>

using std::runtime_error;
using std::cerr;


void test_simple(){
    GafSummary g;
    CigarInterval c;

    c.ref_start = 0;
    c.ref_stop = 10;
    c.code = cigar_char_to_code.at('=');
    c.length = c.ref_stop - c.ref_start;

    g.update_ref("a", 20, c, true);

    c.ref_start = 5;
    c.ref_stop = 15;
    c.code = cigar_char_to_code.at('X');
    c.length = c.ref_stop - c.ref_start;
    g.update_ref("a", 20, c, true);

    for (const auto& [name, summaries]: g.ref_summaries){
        for (const auto& s: summaries) {
            cerr << name << '\n';
            cerr << '\t' << "start:\t" << s.start << '\n';
            cerr << '\t' << "stop:\t" << s.stop << '\n';
            cerr << '\t' << "n_match:\t" << s.n_match << '\n';
            cerr << '\t' << "n_mismatch:\t" << s.n_mismatch << '\n';
            cerr << '\t' << "n_delete:\t" << s.n_delete << '\n';
            cerr << '\t' << "n_insert:\t" << s.n_insert << '\n';
            cerr << '\n';
        }
    }

    g.resolve_all_overlaps();

    for (const auto& [name, summaries]: g.ref_summaries){
        for (const auto& s: summaries) {
            cerr << name << '\n';
            cerr << '\t' << "start:\t" << s.start << '\n';
            cerr << '\t' << "stop:\t" << s.stop << '\n';
            cerr << '\t' << "n_match:\t" << s.n_match << '\n';
            cerr << '\t' << "n_mismatch:\t" << s.n_mismatch << '\n';
            cerr << '\t' << "n_delete:\t" << s.n_delete << '\n';
            cerr << '\t' << "n_insert:\t" << s.n_insert << '\n';
            cerr << '\n';
        }
    }
}


void test_reverse(){
    GafSummary g;
    CigarInterval c;

    c.ref_start = 0;
    c.ref_stop = 10;
    c.code = cigar_char_to_code.at('=');
    c.length = c.ref_stop - c.ref_start;

    g.update_ref("a", 20, c, true);

    c.ref_start = 5;
    c.ref_stop = 15;
    c.code = cigar_char_to_code.at('X');
    c.length = c.ref_stop - c.ref_start;
    g.update_ref("a", 20, c, true);

    for (const auto& [name, summaries]: g.ref_summaries){
        for (const auto& s: summaries) {
            cerr << name << '\n';
            cerr << '\t' << "start:\t" << s.start << '\n';
            cerr << '\t' << "stop:\t" << s.stop << '\n';
            cerr << '\t' << "n_match:\t" << s.n_match << '\n';
            cerr << '\t' << "n_mismatch:\t" << s.n_mismatch << '\n';
            cerr << '\t' << "n_delete:\t" << s.n_delete << '\n';
            cerr << '\t' << "n_insert:\t" << s.n_insert << '\n';
            cerr << '\n';
        }
    }

    g.resolve_all_overlaps();

    for (const auto& [name, summaries]: g.ref_summaries){
        for (const auto& s: summaries) {
            cerr << name << '\n';
            cerr << '\t' << "start:\t" << s.start << '\n';
            cerr << '\t' << "stop:\t" << s.stop << '\n';
            cerr << '\t' << "n_match:\t" << s.n_match << '\n';
            cerr << '\t' << "n_mismatch:\t" << s.n_mismatch << '\n';
            cerr << '\t' << "n_delete:\t" << s.n_delete << '\n';
            cerr << '\t' << "n_insert:\t" << s.n_insert << '\n';
            cerr << '\n';
        }
    }
}


void test_coincidental(){
    GafSummary g;
    CigarInterval c;

    c.ref_start = 0;
    c.ref_stop = 10;
    c.code = cigar_char_to_code.at('=');
    c.length = c.ref_stop - c.ref_start;

    g.update_ref("a", 20, c, true);

    c.ref_start = 5;
    c.ref_stop = 15;
    c.code = cigar_char_to_code.at('X');
    c.length = c.ref_stop - c.ref_start;
    g.update_ref("a", 20, c, true);

    c.ref_start = 10;
    c.ref_stop = 20;
    c.code = cigar_char_to_code.at('D');
    c.length = c.ref_stop - c.ref_start;
    g.update_ref("a", 20, c, true);

    for (const auto& [name, summaries]: g.ref_summaries){
        for (const auto& s: summaries) {
            cerr << name << '\n';
            cerr << '\t' << "start:\t" << s.start << '\n';
            cerr << '\t' << "stop:\t" << s.stop << '\n';
            cerr << '\t' << "n_match:\t" << s.n_match << '\n';
            cerr << '\t' << "n_mismatch:\t" << s.n_mismatch << '\n';
            cerr << '\t' << "n_delete:\t" << s.n_delete << '\n';
            cerr << '\t' << "n_insert:\t" << s.n_insert << '\n';
            cerr << '\n';
        }
    }

    g.resolve_all_overlaps();

    for (const auto& [name, summaries]: g.ref_summaries){
        for (const auto& s: summaries) {
            cerr << name << '\n';
            cerr << '\t' << "start:\t" << s.start << '\n';
            cerr << '\t' << "stop:\t" << s.stop << '\n';
            cerr << '\t' << "n_match:\t" << s.n_match << '\n';
            cerr << '\t' << "n_mismatch:\t" << s.n_mismatch << '\n';
            cerr << '\t' << "n_delete:\t" << s.n_delete << '\n';
            cerr << '\t' << "n_insert:\t" << s.n_insert << '\n';
            cerr << '\n';
        }
    }
}


int main(){
    cerr << "test_simple" << '\n';
    test_simple();

    cerr << "test_coincidental" << '\n';
    test_coincidental();

    return 0;
}
