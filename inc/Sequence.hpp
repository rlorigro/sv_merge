#pragma once

#include <string>

using std::string;


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


char get_complement(char c);

void reverse_complement(string& seq);

}
