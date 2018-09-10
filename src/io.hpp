/*
 * io.hpp
 *
 *  Created on: Oct 12, 2015
 *      Author: jahnka
 *  Modified on: Sep, 2018
 *      Author: pjgeorg
 */

#ifndef IO_H
#define IO_H

#include <algorithm>
#include <array>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "matrices.hpp"

template<class T, std::size_t M, std::size_t S>
static inline auto getInputData(std::string const &fileName)
{
    Matrix<T, S, M> data;
    std::ifstream input(fileName);

    if(!input)
    {
        std::cerr << "Cannot open file " << fileName << std::endl;
        std::abort();
    }

    for(std::size_t m = 0; m < M; ++m)
    {
        for(std::size_t s = 0; s < S; ++s)
        {
            if(!(input >> data[s][m]))
            {
                std::cerr << "Unable to read data from file " << fileName
                          << std::endl;
                std::abort();
            }
        }
    }

    if(!input.eof())
    {
        std::string tmp;
        input >> tmp;
        if(!input.eof())
        {
            std::cerr << "There seems to be more data in input file "
                      << fileName << ". Wrong parameters?" << std::endl;
            std::abort();
        }
    }

    return data;
}

template<class T, std::size_t N, std::size_t G>
static inline auto getTrueTree(
    std::string const &fileName, std::array<std::string, G> const &geneNames)
{
    std::vector<std::string> lines;
    std::ifstream input(fileName);
    if(!input)
    {
        std::cerr << "Cannot open true tree file " << fileName << std::endl;
        std::abort();
    }

    std::string line;
    while(std::getline(input, line))
    {
        if(!line.empty())
        {
            lines.push_back(line);
        }
    }

    Vector<T, N> trueTree;
    std::array<bool, N> check;
    check.fill(false);
    for(auto const &line : lines)
    {
        if(auto found = line.find(" -> "); found != std::string::npos)
        {
            auto parent = line.substr(0, found);
            auto child = line.substr(found + 4, line.size() - found - 5);

            auto childPos =
                std::find(geneNames.begin(), geneNames.end(), child);
            auto parentPos =
                std::find(geneNames.begin(), geneNames.end(), parent);

            if(childPos == geneNames.end() || parentPos == geneNames.end())
            {
                std::cerr << "Unable to find gene name in gene names."
                          << std::endl;
                std::abort();
            }

            auto childIndex = std::distance(geneNames.begin(), childPos);
            auto parentIndex = std::distance(geneNames.begin(), parentPos);

            trueTree[childIndex] = parentIndex;
            check[childIndex] = true;
        }
    }
    if(std::find(check.begin(), check.end(), false) != check.end())
    {
        std::cerr << "Corrupted true tree file " << fileName << "."
                  << std::endl;
        std::abort();
    }

    return trueTree;
}

template<std::size_t M>
static inline auto getGeneNames(std::string const &fileName)
{
    std::array<std::string, M + 1> names;
    std::ifstream input(fileName);

    if(!input)
    {
        std::cerr << "Cannot open gene names file " << fileName << std::endl;
        std::abort();
    }

    for(std::size_t i = 0; i < M; ++i)
    {
        input >> names[i];
    }
    names[M] = "Root";

    if(!input.eof())
    {
        std::cerr << "There seems to be more data in gene names file "
                  << fileName << " than expected. Corrupted file?" << std::endl;
        std::abort();
    }

    return names;
}

template<std::size_t M>
static inline auto getGeneNames()
{
    std::array<std::string, M + 1> names;

    for(std::size_t i = 0; i <= M; ++i)
    {
        std::stringstream id;
        id << i + 1;
        names[i] = id.str();
    }

    return names;
}

template<char ScoreType>
static inline auto getOutputFileName(
    std::size_t i, std::string const &prefix, std::string const &extension)
{
    std::stringstream fileName;
    fileName << prefix;

    if constexpr(ScoreType == 'm')
    {
        fileName << "_ml";
    }
    else if constexpr(ScoreType == 's')
    {
        fileName << "_map";
    }
    else
    {
        static_assert(ScoreType != ScoreType, "Unknown score type.");
    }

    fileName << i << extension;
    return fileName.str();
}

static inline auto writeToFile(
    std::string const &content, std::string const &fileName)
{
    if(!content.empty())
    {
        std::ofstream output(fileName);

        if(!output)
        {
            std::cerr << "Cannot open output file." << fileName << std::endl;
            std::abort();
        }

        output << content;
    }
}

// Write sample taken in MCMC to file
static inline auto writeSampleOutput(
    std::vector<std::string> &samples, std::string const &fileName)
{
    if(!samples.empty())
    {
        std::ofstream output(fileName);

        if(!output)
        {
            std::cerr << "Cannot open samples output file." << fileName
                      << std::endl;
            std::abort();
        }

        for(auto const &sample : samples)
        {
            output << sample;
        }

        samples.clear();
    }
}

// Converts a tree given as lists of children to the Newick tree format
// Note: Only works if the recursion is started with the root = M
template<class T>
static inline std::string getNewickCode(
    std::vector<std::vector<T>> const &children, std::size_t const root)
{
    std::stringstream newick;
    auto const &rootChildren = children[root];

    if(!rootChildren.empty())
    {
        newick << "(";
        bool first = true;

        for(auto const &child : rootChildren)
        {
            if(!first)
            {
                newick << ",";
            }
            first = false;

            newick << getNewickCode(children, child);
        }
        newick << ")";
    }
    newick << root + 1;

    return newick.str();
}

template<class T>
static inline auto getNewickCode(std::vector<std::vector<T>> const &children)
{
    return getNewickCode(children, children.size() - 1) + "\n";
}

// Creates the content for the GraphViz file from a tree parent vector using the
// gene names as node labels
template<class T, class F, std::size_t M, std::size_t S, std::size_t N>
static inline auto getGraphViz(Vector<T, N> const &tree,
    std::array<std::string, M + 1> const &names, bool const attachSamples,
    Matrix<bool, N, N> const &ancestors, Matrix<F, 4, 2> const &logScores,
    Matrix<T, S, M> const &data)
{
    std::stringstream content;
    content << "digraph G {\n";
    content << "node [color=deeppink4, style=filled, fontcolor=white];\n";

    for(std::size_t i = 0; i < M; ++i)
    {
        content << names[tree[i]] << " -> " << names[i] << ";\n";
    }

    if(attachSamples)
    {
        content << "node [color=lightgrey, style=filled, fontcolor=black];\n";
        auto attachment =
            getBestAttachmentString(ancestors, logScores, data, names);
        content << attachment;
    }

    content << "}\n";

    return content.str();
}

// Creates the attachment string for the samples, the optimal attachment points
// are recomputed from scratch based on error log scores
template<class T, class F, std::size_t M, std::size_t S, std::size_t N>
static inline auto getBestAttachmentString(Matrix<bool, N, N> const &ancestors,
    Matrix<F, 4, 2> const &logScores, Matrix<T, S, M> const &data,
    std::array<std::string, M + 1> const &names)
{
    auto points = getBestAttachmentPoints(ancestors, logScores, data);

    std::stringstream attachment;

    for(std::size_t m = 0; m <= M; ++m)
    {
        for(std::size_t s = 0; s < S; ++s)
        {
            if(points[m][s] == true)
            {
                attachment << names[m] << " -> s" << s << ";\n";
            }
        }
    }

    return attachment.str();
}

// Recomputes the best attachment points of the samples to a tree.
template<class T, class F, std::size_t M, std::size_t S, std::size_t N>
static inline auto getBestAttachmentPoints(Matrix<bool, N, N> const &ancestors,
    Matrix<F, 4, 2> const &logScores, Matrix<T, S, M> const &data)
{
    Matrix<bool, M + 1, S> points;
    points.fill(false);

    for(std::size_t s = 0; s < S; ++s)
    {
        F bestScore = 0.0;

        // Attach node to root (no genes mutated)
        for(std::size_t m = 0; m < M; ++m)
        {
            bestScore += logScores[data[s][m]][0];
        }

        // Loop over all attachment points (genes)
        for(std::size_t parent = 0; parent < M; ++parent)
        {
            F score = 0.0;

            for(std::size_t m = 0; m < M; ++m)
            {
                score += logScores[data[s][m]][ancestors[m][parent]];
            }

            bestScore = std::max(score, bestScore);
        }

        bool rootAttachment = true;
        for(std::size_t parent = 0; parent < M; ++parent)
        {
            if(points[parent][s])
            {
                rootAttachment = false;
                break;
            }
        }

        if(rootAttachment)
        {
            points[M][s] = true;
        }
    }

    return points;
}
#endif // IO_H
