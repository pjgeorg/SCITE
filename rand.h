/*
 * rand.h
 *
 *  Created on: Mar 27, 2015
 *      Author: jahnka
 *  Modified on: Mar, 2018
 *      Author: pjgeorg
 */

#ifndef RAND_H
#define RAND_H

#include <random>

#include "matrices.h"

class RandomGenerator
{
public:
    template<class T>
    auto seed(T const &seed)
    {
        sEngine.seed(seed);
    }

    auto &operator()()
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
    return getRandomNumber<T>(0, size - 1);
}

template<class T>
inline auto sampleTwoElementsWithoutReplacement(T const low, T const high)
{
    StaticArray<T, 2> values;
    values[0] = getRandomNumber(low, high);
    do
    {
        values[1] = getRandomNumber(low, high);
    } while(values[0] == values[1]);

    return values;
}

template<class T>
inline auto sampleTwoElementsWithoutReplacement(T const size)
{
    return sampleTwoElementsWithoutReplacement<T>(0, size - 1);
}

template<class T>
inline auto changeBeta(T const prob)
{
    if(getRandomNumber<T>(0, 1) <= prob)
    {
        return true;
    }
    else
    {
        return false;
    }
}

/* This function gets a number of nodes n, and creates a random pruefer code for
 * a rooted tree with n+1 nodes (root is always node n+1) */
template<class T = int, class U = std::size_t>
inline auto getRandTreeCode(U const n)
{
    // as usual n is the number of mutations
    // #nodes = n mutations plus root (wildtype) = n+1
    // #codeLength = #nodes-2 = n-1
    DynamicArray<T> code(n - 1);
    for(int i = 0; i < code.size(); ++i)
    {
        code[i] = getRandomNumber<T>(n + 1);
    }
    return code;
}

// This creates the parent vector of a random binary tree. Entries 0...m-1 are
// for the leafs. Entries m....2m-3 are for the inner nodes except the root, the
// root has index 2m-2 which has no parent and therefore has no entry in the
// parent vector
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

// picks randomly one of the tree moves based on the move probabilities
template<class T>
inline auto sampleRandomMove(DynamicArray<T> const &prob)
{
    auto percent = getRandomNumber<T>(0, 1);
    auto probSum = prob[1];

    // start at index 1; the probability at prob[0] is for changing the error
    // rate (which is treated separately)
    for(std::size_t i = 1; i < prob.size() - 1; ++i)
    {
        if(percent <= probSum)
        {
            return i;
        }
        probSum += prob[i + 1];
    }
    return prob.size() - 1;
}

#endif // RAND_H
