//
// Created by sylwester on 10/8/20.
//
#include <mutex>
#include <vector>
#include <future>
#include <cmath>
#include <iostream>
#include <zconf.h>
#include "Utils/WorkloadManager.h"

void WorkloadManager::parallelBlockExecution(unsigned n, unsigned N, unsigned blocks, unsigned threads,
                                             function<void(unsigned, unsigned, unsigned)> job) {
    unsigned N0 = N - n;

    if (blocks < 1) blocks = 1;
    if (blocks > N0) blocks = N0;

    unsigned W = (unsigned) ceil((double) N0 / blocks);

    vector<mutex> locks(blocks);

    auto worker = [&](int thread_id) {
        for (unsigned i = 0; i < blocks; i++) {
            auto &l = locks[i];
            if (l.try_lock()) {
                unsigned a = i * W;
                unsigned b = min((i + 1) * W - 1, N0);
                a += n;
                b += n;
                job(a, b, thread_id);
//                l.unlock(); // do not unlock! otherwise the same block may be processed many times!
            }
        }
    };

    vector<future<void> > futures(threads - 1);
    for (unsigned i = 1; i < threads; i++) {
        futures[i - 1] = std::async(std::launch::async, worker, i);
    }
    worker(0);
    for (auto &f : futures) f.get();
}


void WorkloadManager::test() {

    int N = 100;
    int T = 5;
    int B = 50;

    auto job = [](unsigned a, unsigned b, unsigned id) {
        cerr << "Executing job in thread id = " << id << "  a = " << a << "  b = " << b << endl;
        sleep(5);
    };

    parallelBlockExecution(0, N, B, T, job);


}