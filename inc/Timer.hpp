#ifndef GFASE_DURATION_HPP
#define GFASE_DURATION_HPP

#include <iostream>
#include <chrono>

using std::chrono::duration_cast;
using std::chrono::steady_clock;
using std::chrono::duration;

using std::chrono::hours;
using std::chrono::minutes;
using std::chrono::seconds;
using std::chrono::milliseconds;

using std::ostream;
using std::string;


namespace sv_merge{

class Timer{
    steady_clock::time_point start;

public:
    Timer();
    string elapsed() const;
    string to_csv() const;
    void reset();
};

}

ostream& operator<<(ostream& o, const sv_merge::Timer& t);


#endif //GFASE_DURATION_HPP
