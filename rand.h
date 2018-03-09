/*
 * rand.h
 *
 *  Created on: Mar 27, 2015
 *      Author: jahnka
 */

#ifndef RAND_H
#define RAND_H

#include <vector>
#include <random>
//#include <string>
#include "matrices.h"

bool changeBeta(double prob);
int sampleRandomMove(std::vector<double> prob);
int* sampleTwoElementsWithoutReplacement(int n);
int pickRandomNumber(int n);
double sample_0_1();
int* getRandTreeCode(int n);
bool samplingByProb(double prob);

class RandomGenerator
{
    public:
        template<class T>
        auto seed(T const &seed)
        {
            sEngine.seed(seed);
        }

        auto& operator()()
        {
            return sEngine;
        }

    private:
        static std::default_random_engine sEngine;
};

template<class T>
inline auto getRandomNumber(T const low, T const high)
{
    auto dist = [&]() {
        if constexpr(std::is_integral<T>())
        {
            return std::uniform_int_distribution<T>(low, high);
        }
        else if constexpr(std::is_floating_point<T>())
        {
            return std::uniform_real_distribution<T>(low, high);
        }
    }();

    return dist(RandomGenerator()());
}

template<class T>
inline auto getRandomNumber(T const size)
{
    return getRandomNumber<T>(0, size-1);
}

// This creates the parent vector of a random binary tree. Entries 0...m-1 are for the leafs.
// Entries m....2m-3 are for the inner nodes except the root, the root has index 2m-2 which has no parent
// and therefore has no entry in the parent vector
template<class T = int, class U = std::size_t>
inline auto getRandomBinaryTree(U count)
{
    auto leafsAndInnerNodesParents = DynamicArray<T>((2 * count) - 2);

    DynamicArray<U> queue(count);
    for(int i = 0; i < queue.size(); i++)
    {
        queue[i] = i;
    }

    while(queue.size() > 1)
    {
        auto pos = getRandomNumber<U>(queue.size());
        auto child1 = queue[pos];
        queue[pos] = queue.back();
        queue.pop_back();

        pos = getRandomNumber<U>(queue.size());
        auto child2 = queue[pos];
        queue[pos] = count;

        leafsAndInnerNodesParents[child1] = count;
        leafsAndInnerNodesParents[child2] = count;
        ++count;
    }

    return leafsAndInnerNodesParents;
}

#endif
