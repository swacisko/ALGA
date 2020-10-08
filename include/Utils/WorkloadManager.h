//
// Created by sylwester on 10/8/20.
//

#ifndef ALGA_WORKLOADMANAGER_H
#define ALGA_WORKLOADMANAGER_H

#include <functional>

using namespace std;

class WorkloadManager {
public:

    /**
     * This is used to parallelly execute jobs on given interval [n,N].
     * The interval [n,N] is divided into [blocks] roughly equal blocks, then each block is processed using given function [job] by one of [threads] threads.
     * @param n
     * @param N
     * @param blocks
     * @param threads
     */
    static void parallelBlockExecution(unsigned n, unsigned N, unsigned blocks, unsigned threads,
                                       function<void(unsigned, unsigned, unsigned)> job);

    static void test();

private:

};

#endif //ALGA_WORKLOADMANAGER_H
