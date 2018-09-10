/*
 * stopwatch.hpp
 *
 *  Created on: Sep, 2018
 *      Author: pjgeorg
 */

#ifndef STOPWATCH_H
#define STOPWATCH_H

extern "C" {
#include <time.h>
}

#include <iostream>

template<class T = double>
static inline auto getTimeInSeconds()
{
    timespec timeSpec;
    if(clock_gettime(CLOCK_MONOTONIC, &timeSpec))
    {
        std::cerr << "Unable to get current time." << std::endl;
        std::abort();
    }
    return T(timeSpec.tv_nsec) / T(1e9) + T(timeSpec.tv_sec);
}

template<class T = double>
class Stopwatch
{
public:
    Stopwatch(std::size_t const flops = 0, std::size_t const bytes = 0)
        : mFlops(flops), mBytes(bytes)
    {
    }
    Stopwatch(Stopwatch const &) = delete;
    Stopwatch(Stopwatch &&) = default;
    Stopwatch &operator=(Stopwatch const &) = delete;
    Stopwatch &operator=(Stopwatch &&) = default;
    ~Stopwatch() = default;

    auto start()
    {
        mSeconds -= getTimeInSeconds();
    }

    auto stop()
    {
        mSeconds += getTimeInSeconds();
    }

    auto addFlops(std::size_t const flops)
    {
        mFlops += flops;
    }

    auto addBytes(std::size_t const bytes)
    {
        mBytes += bytes;
    }

    auto seconds() const
    {
        return mSeconds;
    }

    auto flops() const
    {
        return mFlops;
    }

    auto flopsPerSecond() const
    {
        return T(flops()) / seconds() / T(1e9);
    }

    auto bytes() const
    {
        return mBytes;
    }

    auto bytesPerSecond() const
    {
        return T(bytes()) / seconds() / T(1e9);
    }

private:
    T mSeconds = 0.0;
    std::size_t mFlops;
    std::size_t mBytes;
};
#endif // STOPWATCH_H
