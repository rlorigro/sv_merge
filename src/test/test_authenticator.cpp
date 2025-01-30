#include "Authenticator.hpp"
#include "Timer.hpp"
#include "bam.hpp"

using namespace sv_merge;

#include <thread>
#include <atomic>
#include <exception>

using std::function;
using std::thread;
using std::atomic;
using std::exception;


void thread_fn(int64_t thread_id, mutex& m, Timer& t, Authenticator& authenticator){
    string region_string = "chr1:10000000-10005000";
    path bam_path = "gs:/fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/231ea80e-bef3-4e96-90a4-07efc80eb523/minimap2/20583855-420b-4027-97a0-b521fa94934e/call-alignAndSortBAM/HG00673.bam";

    int64_t duration = 0;

    while (duration < 60*64) {
        authenticator.try_with_authentication(2, [&]() {
            for_read_in_bam_region(bam_path, region_string, [&](Sequence &s) {
                m.lock();
                cerr << t << ' ' << "thread id: " << thread_id << " result: " << s.name << '\n';
                m.unlock();

                // Stop after the first read
                return;
            });
        });

        sleep(120);
        duration += 120;
    }
}


int main(){
    string command = "gcloud auth print-access-token";
    string token;

    // Will throw error if fails, otherwise returns token
    run_command(command, token);

    Authenticator authenticator;

    Timer t;
    mutex m;

    // Thread-related variables
    size_t n_threads = 32;
    atomic<size_t> job_index = 0;
    vector<thread> threads;

    threads.reserve(n_threads);

    // Launch threads
    for (size_t n=0; n<n_threads; n++) {
        try {
            cerr << "launching: " << n << '\n';
            threads.emplace_back(thread_fn,
                                 n,
                                 std::ref(m),
                                 std::ref(t),
                                 std::ref(authenticator)
            );
        } catch (const exception &e) {
            throw e;
        }
    }

    // Wait for threads to finish
    for (auto &n: threads) {
        n.join();
    }

}
