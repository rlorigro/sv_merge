#include "Timer.hpp"

using std::to_string;

namespace sv_merge{

Timer::Timer():
        start(std::chrono::steady_clock::now())
{}


milliseconds Timer::elapsed_milliseconds() const {
    milliseconds d = duration_cast<milliseconds>(std::chrono::steady_clock::now() - start);
    return d;
}


string Timer::elapsed() const {
    auto d = std::chrono::steady_clock::now() - start;

    const auto h = duration_cast<hours>(d);
    const auto m = duration_cast<minutes>(d - h);
    const auto s = duration_cast<seconds>(d - h - m);
    const auto ms = duration_cast<milliseconds>(d - h - m - s);

    string ds;
    ds.append("[");
    ds.append(to_string(h.count()));
    ds.append("h ");
    ds.append(to_string(m.count()));
    ds.append("m ");
    ds.append(to_string(s.count()));
    ds.append("s ");
    ds.append(to_string(ms.count()));
    ds.append("ms");
    ds.append("] ");

    return ds;
}


string Timer::to_csv() const {
    auto d = std::chrono::steady_clock::now() - start;

    string ds;
    duration_to_csv(d, ds);

    return ds;
}


void Timer::reset() {
    start = std::chrono::steady_clock::now();
}


}

ostream& operator<<(ostream& o, const sv_merge::Timer& t) {
    o << t.elapsed();

    return o;
}


