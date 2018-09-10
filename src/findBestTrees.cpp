/*
 * findBestTrees.cpp
 *
 *  Created on: Mar 27, 2015
 *      Author: jahnka
 *  Modified on: Aug, 2018
 *      Author: pjgeorg
 */

extern "C" {
#include <time.h>
}

#include <array>
#include <iostream>

#include "io.hpp"
#include "matrices.hpp"
#include "mcmc.hpp"
#include "parameter.hpp"
#include "stopwatch.hpp"

#include "std.hpp"

int main(int argc, char *argv[])
{
    ///////////////////////////////////////////////////////////////////////////

    using Integer = int;
    using Float = double;

#ifdef MUTATIONS
    constexpr auto cMutations = MUTATIONS;
#else
#error Number of mutations (MUTATIONS) not defined.
#endif // MUTATIONS

#ifdef SAMPLES
    constexpr auto cSamples = SAMPLES;
#else
#error Number of samples (SAMPLES) not defined.
#endif // SAMPLES

    // 'm': mutation tree, 't': rooted binary leaf-labelled tree
    constexpr auto cTreeType =
#ifdef TREE_TYPE
        TREE_TYPE;
#else
        'm';
#endif // TREE_TYPE

    // 'm': default, 's': Sample attachment points are marginalized out
    constexpr auto cScoreType =
#ifdef SCORE_TYPE
        SCORE_TYPE;
#else
        'm';
#endif // SCORE_TYPE

    // True mutation tree is known
    constexpr auto cTrueTree =
#ifdef TRUE_TREE
        true;
#else
        false;
#endif // TRUE_TREE

    ///////////////////////////////////////////////////////////////////////////

    static_assert(cTreeType == 'm' || cTreeType == 't', "Unknown TreeType.");
    static_assert(cScoreType == 'm' || cScoreType == 's', "Unknown ScoreType.");

    // Size of mutation/binary tree
#if __cpp_constexpr >= 201603
    constexpr auto cTreeSize = []() {
        if constexpr(cTreeType == 'm')
        {
            return cMutations;
        }
        else if constexpr(cTreeType == 't')
        {
            return (2 * cSamples) - 2;
        }
    }();
#else
    constexpr auto cTreeSize =
        cTreeType == 'm' ? cMutations : (2 * cSamples) - 2;
#endif // __cpp_constexpr

    Stopwatch<> stopwatch;
    stopwatch.start();

    //  Get command line parameters
    auto inputFile = parameterOption<std::string>(argv, argv + argc, "-i");
    auto repetitions = parameterOption<std::size_t>(argv, argv + argc, "-r");
    auto chainLength = parameterOption<std::size_t>(argv, argv + argc, "-l");
    auto fd = parameterOption<Float>(argv, argv + argc, "-fd");
    auto ad = parameterOption<Float, 2>(argv, argv + argc, "-ad");

    auto cc = [&argv, &argc]() {
        if(parameterExists(argv, argv + argc, "-cc"))
        {
            return parameterOption<Float>(argv, argv + argc, "-cc");
        }
        else
        {
            return Float(0.0);
        }
    }();

    auto errorRates = std::array<Float, 4>{fd, ad[0], ad[1], cc};

    auto sampleStep = [&argv, &argc]() {
        if(parameterExists(argv, argv + argc, "-p"))
        {
            return parameterOption<std::size_t>(argv, argv + argc, "-p");
        }
        else
        {
            return std::size_t{0};
        }
    }();

    auto errorRateMove = [&argv, &argc]() {
        if(parameterExists(argv, argv + argc, "-e"))
        {
            return parameterOption<Float>(argv, argv + argc, "-e");
        }
        else
        {
            return Float(0.0);
        }
    }();

    auto chi = [&argv, &argc]() {
        if(parameterExists(argv, argv + argc, "-x"))
        {
            return parameterOption<Float>(argv, argv + argc, "-x");
        }
        else
        {
            return Float(10.0);
        }
    }();

    auto priorStd = [&argv, &argc]() {
        if(parameterExists(argv, argv + argc, "-sd"))
        {
            return parameterOption<Float>(argv, argv + argc, "-sd");
        }
        else
        {
            return Float(0.1);
        }
    }();

    auto outputName = [&argv, &argc, &inputFile]() {
        if(parameterExists(argv, argv + argc, "-o"))
        {
            return parameterOption<std::string>(argv, argv + argc, "-o");
        }
        else
        {
            return inputFile.substr(0, inputFile.find_last_of("."));
        }
    }();

    auto sampleFile = [&argv, &argc, &outputName]() {
        std::stringstream fileName;
        fileName << outputName << ".samples";
        return fileName.str();
    }();

    auto geneNames = [&argv, &argc]() {
        if(parameterExists(argv, argv + argc, "-names"))
        {
            auto fileName =
                parameterOption<std::string>(argv, argv + argc, "-names");
            return getGeneNames<cMutations>(fileName);
        }
        else
        {
            return getGeneNames<cMutations>();
        }
    }();

    auto attachSamples = parameterExists(argv, argv + argc, "-a");

    auto maxTreeListSize = [&argv, &argc]() {
        if(parameterExists(argv, argv + argc, "-max_treelist_size"))
        {
            return parameterOption<std::size_t>(
                argv, argv + argc, "-max_treelist_size");
        }
        else
        {
            return std::numeric_limits<std::size_t>::max();
        }
    }();

    auto gamma = [&argv, &argc]() {
        if(parameterExists(argv, argv + argc, "-g"))
        {
            return parameterOption<Float>(argv, argv + argc, "-g");
        }
        else
        {
            return Float(1.0);
        }
    }();

    auto treeMoves = [&argv, &argc]() {
        if(parameterExists(argv, argv + argc, "-move_probs"))
        {
            auto probs = parameterOption < Float,
                 cTreeType == 'm' ? 3 : 2 > (argv, argv + argc, "-move_probs");
            if(auto sum = std::accumulate(probs.begin(), probs.end(), 0.0);
                sum != 1.0)
            {
                for(auto &v : probs)
                {
                    v /= sum;
                }
            }
            return probs;
        }
        else
        {
            if constexpr(cTreeType == 'm')
            {
                return std::array<Float, 3>{0.55, 0.4, 0.05};
            }
            else
            {
                return std::array<Float, 2>{0.4, 0.6};
            }
        }
    }();

    // cTreeType == 'm':
    //     change beta / prune & re-attach / swap node labels / swap subtrees
    // cTreeType == 't':
    //     change beta / prune & re-attach / swap leaf labels
    std::array<Float, treeMoves.size() + 1> moveProbs;
    moveProbs[0] = errorRateMove;
    for(std::size_t i = 1; i < moveProbs.size(); ++i)
    {
        moveProbs[i] = treeMoves[i - 1];
    }

    // Read data matrix
    auto inputData = getInputData<Integer, cMutations, cSamples>(inputFile);

    auto optimalTrees = runMCMCbeta<cTreeType, cScoreType, cTreeSize>(inputData,
        repetitions, chainLength, errorRates, moveProbs, gamma, chi, priorStd,
        sampleStep, sampleFile);

    std::cout << optimalTrees.size() << " opt trees \n";

    if constexpr(cTrueTree && cTreeType == 'm')
    {
        auto fileName = parameterOption<std::string>(argv, argv + argc, "-t");
        auto trueTree = getTrueTree<Integer, cTreeSize>(fileName, geneNames);
        auto minDistToTrueTree = getMinDistToTrueTree(optimalTrees, trueTree);
        std::cout << "Minimum distance to true tree: " << minDistToTrueTree
                  << std::endl;
    }

    // Output optimal trees found in individual files
    auto logScores = getLogScores(errorRates);

    auto num = std::min(optimalTrees.size(), maxTreeListSize);
    for(std::size_t i = 0; i < num; ++i)
    {
        auto const &tree = optimalTrees[i].first;
        auto const &beta = optimalTrees[i].second;

        // newick
        {
            auto children = getChildren(tree);
            auto fileName =
                getOutputFileName<cScoreType>(i, outputName, ".newick");
            auto newick = getNewickCode(children);
            writeToFile(newick, fileName);
        }

        updateLogScores(logScores, beta);

        // GraphViz
        if constexpr(cTreeType == 'm')
        {
            auto ancestors = getAncestors(tree);
            auto fileName = getOutputFileName<cScoreType>(i, outputName, ".gv");

            auto content = getGraphViz(tree, geneNames, attachSamples,
                ancestors, logScores, inputData);
            writeToFile(content, fileName);
        }
    }

    stopwatch.stop();
    std::cout << "Time elapsed: " << stopwatch.seconds() << "s" << std::endl;

    return 0;
}
