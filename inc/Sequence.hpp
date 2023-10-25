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


}
