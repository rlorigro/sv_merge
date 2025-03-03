#include "gaf.hpp"

using namespace sv_merge;

#include <iostream>
#include <stdexcept>

using std::runtime_error;
using std::cerr;


void test_simple(){
    TransMap t;
    GafSummary g(t);
    CigarInterval c;

    c.ref_start = 0;
    c.ref_stop = 10;
    c.code = cigar_char_to_code.at('=');
    c.length = c.ref_stop - c.ref_start;

    g.update_node("a", c, true);

    c.ref_start = 5;
    c.ref_stop = 15;
    c.code = cigar_char_to_code.at('X');
    c.length = c.ref_stop - c.ref_start;
    g.update_node("a", c, true);

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

    cerr << "---- resolved overlaps ----\n";

    unordered_map<string,float> a1 = {
        {"start",      0},
        {"stop",       5},
        {"n_match",    5},
        {"n_mismatch", 0},
        {"n_delete",   0},
        {"n_insert",   0}
    };

    unordered_map<string,float> a2 = {
        {"start",      5},
        {"stop",       10},
        {"n_match",    2.5},
        {"n_mismatch", 2.5},
        {"n_delete",   0},
        {"n_insert",   0}
    };

    unordered_map<string,float> a3 = {
        {"start",      10},
        {"stop",       15},
        {"n_match",    0},
        {"n_mismatch", 5},
        {"n_delete",   0},
        {"n_insert",   0}
    };

    vector <unordered_map <string,float> > expected_results = {a1,a2,a3};

    for (const auto& [name, summaries]: g.ref_summaries){
        int i = 0;
        for (const auto& s: summaries) {
            cerr << name << ' ' << i << '\n';
            auto r = expected_results[i];

            cerr << '\t' << "start:\t" << s.start << '\n';
            if (r["start"] != s.start){
                cerr << "expected: " << r["start"] << '\n';
                throw runtime_error("FAIL: start mismatch");
            }
            cerr << '\t' << "stop:\t" << s.stop << '\n';
            if (r["stop"] != s.stop){
                cerr << "expected: " << r["stop"] << '\n';
                throw runtime_error("FAIL: stop mismatch");
            }
            cerr << '\t' << "n_match:\t" << s.n_match << '\n';
            if (r["n_match"] != s.n_match){
                cerr << "expected: " << r["n_match"] << '\n';
                throw runtime_error("FAIL: n_match mismatch");
            }
            cerr << '\t' << "n_mismatch:\t" << s.n_mismatch << '\n';
            if (r["n_mismatch"] != s.n_mismatch){
                cerr << "expected: " << r["n_mismatch"] << '\n';
                throw runtime_error("FAIL: n_mismatch mismatch");
            }
            cerr << '\t' << "n_delete:\t" << s.n_delete << '\n';
            if (r["n_delete"] != s.n_delete){
                cerr << "expected: " << r["n_delete"] << '\n';
                throw runtime_error("FAIL: n_delete mismatch");
            }
            cerr << '\t' << "n_insert:\t" << s.n_insert << '\n';
            if (r["n_insert"] != s.n_insert){
                cerr << "expected: " << r["n_insert"] << '\n';
                throw runtime_error("FAIL: n_insert mismatch");
            }
            i++;
            cerr << '\n';
        }
    }
}

void test_coincidental(){
    TransMap t;
    GafSummary g(t);
    CigarInterval c;

    c.ref_start = 0;
    c.ref_stop = 10;
    c.code = cigar_char_to_code.at('=');
    c.length = c.ref_stop - c.ref_start;

    g.update_node("a", c, true);

    c.ref_start = 5;
    c.ref_stop = 15;
    c.code = cigar_char_to_code.at('X');
    c.length = c.ref_stop - c.ref_start;
    g.update_node("a", c, true);

    c.ref_start = 10;
    c.ref_stop = 20;
    c.code = cigar_char_to_code.at('D');
    c.length = c.ref_stop - c.ref_start;
    g.update_node("a", c, true);

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



    cerr << "---- resolved overlaps ----\n";

    unordered_map<string,float> a1 = {
            {"start",      0},
            {"stop",       5},
            {"n_match",    5},
            {"n_mismatch", 0},
            {"n_delete",   0},
            {"n_insert",   0}
    };

    unordered_map<string,float> a2 = {
            {"start",      5},
            {"stop",       10},
            {"n_match",    2.5},
            {"n_mismatch", 2.5},
            {"n_delete",   0},
            {"n_insert",   0}
    };

    unordered_map<string,float> a3 = {
            {"start",      10},
            {"stop",       15},
            {"n_match",    0},
            {"n_mismatch", 2.5},
            {"n_delete",   2.5},
            {"n_insert",   0}
    };

    unordered_map<string,float> a4 = {
            {"start",      15},
            {"stop",       20},
            {"n_match",    0},
            {"n_mismatch", 0},
            {"n_delete",   5},
            {"n_insert",   0}
    };

    vector <unordered_map <string,float> > expected_results = {a1,a2,a3,a4};

    for (const auto& [name, summaries]: g.ref_summaries){
        int i = 0;
        for (const auto& s: summaries) {
            cerr << name << ' ' << i << '\n';
            auto r = expected_results[i];

            cerr << '\t' << "start:\t" << s.start << '\n';
            if (r["start"] != s.start){
                cerr << "expected: " << r["start"] << '\n';
                throw runtime_error("FAIL: start mismatch");
            }
            cerr << '\t' << "stop:\t" << s.stop << '\n';
            if (r["stop"] != s.stop){
                cerr << "expected: " << r["stop"] << '\n';
                throw runtime_error("FAIL: stop mismatch");
            }
            cerr << '\t' << "n_match:\t" << s.n_match << '\n';
            if (r["n_match"] != s.n_match){
                cerr << "expected: " << r["n_match"] << '\n';
                throw runtime_error("FAIL: n_match mismatch");
            }
            cerr << '\t' << "n_mismatch:\t" << s.n_mismatch << '\n';
            if (r["n_mismatch"] != s.n_mismatch){
                cerr << "expected: " << r["n_mismatch"] << '\n';
                throw runtime_error("FAIL: n_mismatch mismatch");
            }
            cerr << '\t' << "n_delete:\t" << s.n_delete << '\n';
            if (r["n_delete"] != s.n_delete){
                cerr << "expected: " << r["n_delete"] << '\n';
                throw runtime_error("FAIL: n_delete mismatch");
            }
            cerr << '\t' << "n_insert:\t" << s.n_insert << '\n';
            if (r["n_insert"] != s.n_insert){
                cerr << "expected: " << r["n_insert"] << '\n';
                throw runtime_error("FAIL: n_insert mismatch");
            }
            i++;
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
