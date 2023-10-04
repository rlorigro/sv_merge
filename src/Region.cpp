#include "Region.hpp"

#include <iostream>

using std::cerr;


namespace sv_merge {


Region::Region(string& name, int64_t start, int64_t stop):
    name(name),
    start(start),
    stop(stop)
{}


Region::Region(string& region_string) {
    bool found_colon = false;
    bool found_hyphen = false;

    string start_token;
    string stop_token;

    for (char &c: region_string) {
//        cerr << c << ' ' << found_colon << ' ' << found_hyphen << ' ' << start_token << ' ' << stop_token << '\n';

        if (c == ':') {
            found_colon = true;
            continue;
        }
        else if (c == '-') {
            found_hyphen = true;
            continue;
        }

        if (not found_colon) {
            name += c;
        }
        else {
            if (not found_hyphen) {
                start_token += c;
            }
            else {
                stop_token += c;
            }
        }
    }

    cerr << start_token.size() << '\n';

    start = stoi(start_token);
    stop = stoi(stop_token);
}

}
