#pragma once

#include <string>
#include <vector>
#include <BinarySequence.hpp>

using std::string;
using std::vector;


namespace sv_merge {

class Sequence{
public:
    string name;
    string sequence;

    Sequence(const string& name, const string& sequence);
    Sequence()=default;
};


class StrandedSequence{
public:
    string name;
    string sequence;
    bool is_reverse;

    StrandedSequence(const string& name, const string& sequence, bool is_reverse);
    StrandedSequence(const string& name, const string& sequence);
    StrandedSequence();
};


class StrandedQSequence{
public:
    BinarySequence<uint64_t> sequence;
    string name;
    vector<uint8_t> qualities;
    bool is_reverse;
    string tags;

    StrandedQSequence(const string& name, const string& sequence, const vector<uint8_t>& qualities, bool is_reverse);
    StrandedQSequence(const string& name, const string& sequence, const vector<uint8_t>& qualities);
    StrandedQSequence(const string& name, const string& sequence, bool is_reverse);
    StrandedQSequence(const string& name, const string& sequence);
    StrandedQSequence(const string& name, const BinarySequence<uint64_t>& sequence);
    StrandedQSequence();
};


char get_complement(char c);

void reverse_complement(string& seq);

}
