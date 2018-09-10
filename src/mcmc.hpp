/*
 * mcmc.hpp
 *
 *  Created on: Mar 27, 2015
 *      Author: jahnka
 *  Modified on: Sep, 2018
 *      Author: pjgeorg
 */

#ifndef MCMC_H
#define MCMC_H

#include "move_tree.hpp"
#include "rand_tree.hpp"
#include "score_tree.hpp"
#include "tree.hpp"
#include "treelist.hpp"

template<class F>
static inline auto logBetaPDF(F const x, std::array<F, 2> const &betaPrior)
{
    // f(x,a,b) = gamma(a+b)//gamma(a)gamma(b)) * x^(a-1) * (1-x)^(b-1)
    return std::log(std::tgamma(betaPrior[0] + betaPrior[1])) +
        (betaPrior[0] - 1) * std::log(x) +
        (betaPrior[1] - 1) * std::log(1 - x) -
        std::log(std::tgamma(betaPrior[0])) -
        std::log(std::tgamma(betaPrior[1]));
}

template<class T>
static inline auto changeBeta(RNG &rng, T const prob)
{
    if(getRandomNumber(rng, T{0}, T{1}) <= prob)
    {
        return true;
    }
    else
    {
        return false;
    }
}

template<class T>
static inline auto proposeNewBeta(RNG &rng, T const current, T const jumpStd)
{
    auto proposed = std::abs(current + sampleNormal(rng, T{0}, jumpStd));

    if(proposed > 1)
    {
        return proposed - 2 * (proposed - 1);
    }

    return proposed;
}

template<class T, class F, std::size_t M>
static inline auto sampleFromPosterior(Vector<T, M> const &tree,
    F const betaProb, F const beta, F const score, F const logScore)
{
    std::stringstream content;
    content << logScore << "\t";
    content << countBranches(tree);
    if(betaProb > 0.0)
    {
        content << "\t" << beta;
        content << "\t" << score;
    }
    content << "\t";
    for(std::size_t i = 0; i < M; ++i)
    {
        content << tree[i] << " ";
    }
    content << "\n";

    return content.str();
}

template<char TreeType, char ScoreType, std::size_t N, class T, class F,
    std::size_t M, std::size_t S, std::size_t P>
static inline auto runMCMCbeta(Matrix<T, S, M> const &data,
    std::size_t const repetitions, std::size_t const chainLength,
    std::array<F, 4> const &errorRates, std::array<F, P> const &moveProbs,
    F const gamma, F const chi, F const priorStd, std::size_t const sampleStep,
    std::string const &sampleFile)
{
    std::size_t const burnIn = 0.25 * chainLength;
    auto const betaPriorMean =
        errorRates[1] + errorRates[2]; // AD1 + AD2 prior mean for AD error rate
    auto const &betaPriorStd = priorStd; // prior std for AD error rate

    // chi: scaling of the known error rate for the MH
    auto const jumpStd = betaPriorStd / chi;

    // Turn the mean and std into parameters of the beta distribution
    std::array<F, 2> betaPrior;
    betaPrior[0] = ((1 - betaPriorMean) * betaPriorMean * betaPriorMean /
                       (betaPriorStd * betaPriorStd)) -
        betaPriorMean;
    betaPrior[1] = betaPrior[0] * ((1 / betaPriorMean) - 1);

    std::size_t optStatesAfterBurnIn = 0;

    auto bestTreeLogScore = std::numeric_limits<F>::lowest();
    auto bestScore = std::numeric_limits<F>::lowest();
    auto bestBeta = betaPriorMean;
    std::vector<std::pair<Vector<T, N>, F>> optTrees;

    RNG rng(std::random_device{}());
    std::vector<std::string> sampleOutput;

    for(std::size_t r = 0; r < repetitions; ++r)
    {
        // Start MCMC with random tree
        auto currTree = getRandomTree<TreeType, T, N>(rng);
        auto currAncestors = getAncestors(currTree);
        auto currLogScores = getLogScores(errorRates);
        auto currBeta = betaPriorMean;
        auto currTreeLogScore =
            getTreeLogScore<TreeType, ScoreType>(data, currTree, currLogScores);

        // Zero if beta is fixed
        auto currBetaLogScore =
            (moveProbs[0] == 0) ? F{0} : logBetaPDF(currBeta, betaPrior);

        // Combined score of current tree and current beta
        auto currScore = currTreeLogScore + currBetaLogScore;

        for(std::size_t it = 0; it < chainLength; ++it)
        {
            if(changeBeta(rng, moveProbs[0])) // Move changes beta, not the tree
            {
                // New beta is proposed, log scores change tree is copy of
                // current tree
                auto propBeta = proposeNewBeta(rng, currBeta, jumpStd);
                auto propLogScores = currLogScores;
                updateLogScores(propLogScores, propBeta);
                auto propBetaLogScore = logBetaPDF(propBeta, betaPrior);

                auto propTreeLogScore = getTreeLogScore<TreeType, ScoreType>(
                    data, currTree, propLogScores, bestTreeLogScore);

                if(getRandomNumber<F>(rng, 0, 1) <
                    std::exp((propTreeLogScore + propBetaLogScore -
                                 currTreeLogScore - currBetaLogScore) *
                        gamma))
                {
                    // Proposed move has been accepted
                    currTreeLogScore = propTreeLogScore;
                    currBeta = propBeta;
                    currBetaLogScore = propBetaLogScore;
                    currScore = currTreeLogScore + currBetaLogScore;
                    currLogScores = propLogScores;
                }
            }
            else // Move changed tree
            {
                F correction = 1.0;
                auto propTree = proposeNewTree<TreeType>(
                    rng, currTree, currAncestors, moveProbs, correction);

                auto propTreeLogScore = getTreeLogScore<TreeType, ScoreType>(
                    data, propTree, currLogScores, bestTreeLogScore);

                if(getRandomNumber<F>(rng, 0, 1) < correction *
                        std::exp((propTreeLogScore - currTreeLogScore) * gamma))
                {
                    // Proposed tree has been accepted
                    currAncestors = getAncestors(propTree);
                    currTree = propTree;
                    currTreeLogScore = propTreeLogScore;
                    currScore = currTreeLogScore + currBetaLogScore;
                }
            }

            // Sample from the posterior if required and past the burn-in phase
            if(sampleStep > 1 && it >= burnIn && it % sampleStep == 0)
            {
                sampleOutput.push_back(sampleFromPosterior(currTree,
                    moveProbs[0], currBeta, currScore, currTreeLogScore));
            }

            // Update best tree in case we have found a new best one
            if(currScore > bestScore)
            {
                resetTreeList(optTrees, currTree, currBeta);
                optStatesAfterBurnIn = 1;
                bestTreeLogScore = currTreeLogScore;
                bestScore = currScore;
                bestBeta = currBeta;
            }
            else if(currScore == bestScore)
            {
                addToTreeList(optTrees, currTree, currBeta);

                // Update the number of MCMC steps we spent in an optimal state
                if(it >= burnIn)
                {
                    ++optStatesAfterBurnIn;
                }
            }
            if((it + 1) % 100000 == 0)
            {
                std::cout << "At mcmc repetition " << r + 1 << "/"
                          << repetitions << ", step " << it + 1
                          << ": best tree score " << bestTreeLogScore
                          << " and best beta " << bestBeta
                          << " and best overall score " << bestScore
                          << std::endl;
                writeSampleOutput(sampleOutput, sampleFile);
            }
        }
    }

    writeSampleOutput(sampleOutput, sampleFile);

    auto noStepsAfterBurnin = repetitions * (chainLength - burnIn);
    std::cout.precision(17);
    std::cout << "best log score for tree:\t" << bestTreeLogScore << std::endl;
    std::cout << "#optimal steps after burn-in:\t" << optStatesAfterBurnIn
              << std::endl;
    std::cout << "total #steps after burn-in:\t" << noStepsAfterBurnin
              << std::endl;
    std::cout << "%optimal steps after burn-in:\t"
              << (1.0 * optStatesAfterBurnIn) / noStepsAfterBurnin << std::endl;
    if(moveProbs[0] != 0.0)
    {
        std::cout << "best value for beta:\t" << bestBeta << std::endl;
        std::cout << "best log score for (T, beta):\t" << bestScore
                  << std::endl;
    }

    return optTrees;
}
#endif // MCMC_H
