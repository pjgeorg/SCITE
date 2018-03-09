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

bool changeBeta(double prob);
int sampleRandomMove(std::vector<double> prob);
int* sampleTwoElementsWithoutReplacement(int n);
int pickRandomNumber(int n);
double sample_0_1();
int* getRandTreeCode(int n);
bool samplingByProb(double prob);
int* getRandomBinaryTree(int m);

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
#endif
