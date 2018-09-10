/*
 * parameter.hpp
 *
 *  Created on: Sep, 2018
 *      Author: pjgeorg
 */

#ifndef PARAMETER_H
#define PARAMETER_H

#include <algorithm>
#include <array>
#include <sstream>
#include <string>

static inline auto parameterExists(
    char *const *begin, char *const *end, std::string const &parameter)
{
    return std::find(begin, end, parameter) != end;
}

template<class T>
static inline auto parameterOption(
    char **begin, char **end, std::string const &parameter, T &value)
{
    char **input = std::find(begin, end, parameter);
    if(input != end)
    {
        if(++input != end)
        {
            std::istringstream stream(*input);
            stream >> value;
        }
        else
        {
            std::cerr << "Invalid argument: " << parameter << std::endl;
            std::abort();
        }
    }
    else
    {
        std::cerr << "Parameter not found: " << parameter << std::endl;
        std::abort();
    }
}

template<class T, std::size_t N>
static inline auto parameterOption(char **begin, char **end,
    std::string const &parameter, std::array<T, N> &option)
{
    char **input = std::find(begin, end, parameter);
    if(input != end)
    {
        for(auto &value : option)
        {
            if(++input != end)
            {
                std::istringstream stream(*input);
                stream >> value;
            }
            else
            {
                std::cerr << "Invalid argument: " << parameter << std::endl;
                std::abort();
            }
        }
    }
}

template<class T, std::size_t N = 1>
static inline auto parameterOption(
    char **begin, char **end, std::string const &parameter)
{
    std::conditional_t<N == 1, T, std::array<T, N>> value;
    parameterOption(begin, end, parameter, value);
    return value;
}
#endif // PARAMETER_H
