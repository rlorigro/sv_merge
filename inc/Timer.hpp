#pragma once

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
using std::to_string;


namespace sv_merge{

class Timer{
    steady_clock::time_point start;

public:
    Timer();
    string elapsed() const;
    milliseconds elapsed_milliseconds() const;
    string to_csv() const;

    void reset();
};


template <class T> void duration_to_csv(T& duration, string& result) {

    const auto h = duration_cast<hours>(duration);
    const auto m = duration_cast<minutes>(duration - h);
    const auto s = duration_cast<seconds>(duration - h - m);
    const auto ms = duration_cast<milliseconds>(duration - h - m - s);

    result.append(to_string(h.count()));
    result.append(",");
    result.append(to_string(m.count()));
    result.append(",");
    result.append(to_string(s.count()));
    result.append(",");
    result.append(to_string(ms.count()));
}

}

ostream& operator<<(ostream& o, const sv_merge::Timer& t);
