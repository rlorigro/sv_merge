#include <iostream>
#include <stdexcept>
#include <thread>

using std::this_thread::sleep_for;
using std::runtime_error;
using std::cerr;

#ifdef __linux__

#include "bdsg/include/bdsg/internal/hash_map.hpp"
#include "TransitiveMap.hpp"
#include "misc.hpp"

using sv_merge::get_peak_memory_usage;
using sv_merge::TransMap;
using sv_merge::HeteroNode;

void test_full_mapping(size_t n_samples){
    TransMap transmap;

    size_t coverage = 8;
    size_t n_reads = n_samples*coverage;
    size_t n_paths = 100;

    for (size_t r=0; r<n_reads; r++) {
        auto r_name = "r" + std::to_string(r);
        transmap.add_read(r_name);

        size_t s = r / coverage;

        if (r % coverage == 0) {
            auto s_name = "s" + std::to_string(s);
            transmap.add_sample(s_name);
            transmap.add_edge(r_name, s_name);
        }
    }

    for (size_t p=0; p<n_paths; p++) {
        auto p_name = "p" + std::to_string(p);
        transmap.add_path(p_name);
    }

    for (size_t r=0; r<n_reads; r++) {
        auto r_name = "r" + std::to_string(r);

        for (size_t p=0; p<n_paths; p++) {
            auto p_name = "p" + std::to_string(p);

            transmap.add_edge(r_name, p_name);
        }
    }

    sleep_for(std::chrono::seconds(2));
}


void test_sample_read_mapping(size_t n_samples){
    TransMap transmap;

    size_t coverage = 8;
    size_t n_reads = n_samples*coverage;
    size_t n_paths = 100;

    for (size_t r=0; r<n_reads; r++) {
        auto r_name = "r" + std::to_string(r);
        transmap.add_read(r_name);

        size_t s = r / coverage;

        if (r % coverage == 0) {
            auto s_name = "s" + std::to_string(s);
            transmap.add_sample(s_name);
            transmap.add_edge(r_name, s_name);
        }
    }

    sleep_for(std::chrono::seconds(2));
}


int main(){
    size_t n_samples = 1;

    cerr << "n_samples: " << n_samples << '\n';

    test_full_mapping(n_samples);

    cerr << get_peak_memory_usage() << " used" << '\n';

    n_samples *= 10;

    cerr << "n_samples: " << n_samples << '\n';

    test_full_mapping(n_samples);

    cerr << get_peak_memory_usage() << " used" << '\n';

    n_samples *= 10;

    cerr << "n_samples: " << n_samples << '\n';

    test_full_mapping(n_samples);

    cerr << get_peak_memory_usage() << " used" << '\n';

    n_samples *= 10;

    cerr << "n_samples: " << n_samples << '\n';

    test_full_mapping(n_samples);

    cerr << get_peak_memory_usage() << " used" << '\n';

    n_samples *= 10;

    cerr << "n_samples: " << n_samples << '\n';

    test_full_mapping(n_samples);

    cerr << get_peak_memory_usage() << " used" << '\n';

    n_samples *= 10;
}

#else

int main(){
    cerr << "This exe is only supported on Linux" << '\n';
    return 0;
}

#endif