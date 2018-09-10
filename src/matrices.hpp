/*
 * matrices.hpp
 *
 *  Created on: Mar 27, 2015
 *      Author: jahnka
 *  Modified on: Sep, 2018
 *      Author: pjgeorg
 */

#ifndef MATRICES_H
#define MATRICES_H

#include <cstdint>

template<class T, std::size_t N>
class Vector
{
public:
    auto &operator[](std::size_t const i)
    {
        return mData[i];
    }

    auto &operator[](std::size_t const i) const
    {
        return mData[i];
    }

    auto &fill(T const &value)
    {
        for(std::size_t i = 0; i < N; ++i)
        {
            mData[i] = value;
        }
        return *this;
    }

private:
    T mData[N];
};

template<class T, std::size_t M, std::size_t N>
class Matrix
{
public:
    auto &operator[](std::size_t const i)
    {
        return mData[i];
    }

    auto &operator[](std::size_t const i) const
    {
        return mData[i];
    }

    auto &fill(T const &value)
    {
        for(std::size_t i = 0; i < M; ++i)
        {
            mData[i].fill(value);
        }
        return *this;
    }

    auto &operator+=(Matrix const &rhs)
    {
        for(std::size_t i = 0; i < M; ++i)
        {
            for(std::size_t j = 0; j < N; ++j)
            {
                mData[i][j] += rhs[i][j];
            }
        }

        return *this;
    }

private:
    Vector<T, N> mData[M];
};

template<class T, std::size_t N>
static inline auto max(Vector<T, N> const &arg)
{
    auto max = arg[0];
    for(std::size_t i = 1; i < N; ++i)
    {
        max = std::max(max, arg[i]);
    }
    return max;
}

template<class T, std::size_t N>
static inline auto operator==(Vector<T, N> const &lhs, Vector<T, N> const &rhs)
{
    for(std::size_t i = 0; i < N; ++i)
    {
        if(lhs[i] != rhs[i])
        {
            return false;
        }
    }
    return true;
}

template<class T, std::size_t N>
static inline auto operator!=(Vector<T, N> const &lhs, Vector<T, N> const &rhs)
{
    return !(lhs == rhs);
}
#endif // MATRICES_H
