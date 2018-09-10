/*
 * rand.hpp
 *
 *  Created on: Mar 27, 2015
 *      Author: jahnka
 *  Modified on: Sep, 2018
 *      Author: pjgeorg
 */

#ifndef RAND_H
#define RAND_H

#include <array>
#include <random>

using RNG = std::default_random_engine;

template<class T>
static inline auto getRandomNumber(RNG &rng, T const low, T const high)
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

    return dist(rng);
}

template<class T>
static inline auto getRandomNumber(RNG &rng, std::vector<T> const &container)
{
    return getRandomNumber<std::size_t>(rng, 0, container.size() - 1);
}

template<class T>
static inline auto getRandomNumber(RNG &rng, T const size)
{
    return getRandomNumber<T>(rng, T{0}, size - 1);
}

template<class T>
static inline auto sampleNormal(RNG &rng, T const mean, T const stddev)
{
    return std::normal_distribution<T>(mean, stddev)(rng);
}

template<class T>
static inline auto sampleTwoElementsWithoutReplacement(
    RNG &rng, T const low, T const high)
{
    std::array<T, 2> values;
    values[0] = getRandomNumber(rng, low, high);
    do
    {
        values[1] = getRandomNumber(rng, low, high);
    } while(values[0] == values[1]);

    return values;
}

template<class T>
static inline auto sampleTwoElementsWithoutReplacement(RNG &rng, T const size)
{
    return sampleTwoElementsWithoutReplacement<T>(rng, T{0}, size - 1);
}
#endif // RAND_H
