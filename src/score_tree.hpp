/*
 * score_tree.hpp
 *
 *  Created on: Aug 16, 2015
 *      Author: jahnka
 *  Modified on: Sep, 2018
 *      Author: pjgeorg
 */

#ifndef SCORE_TREE_H
#define SCORE_TREE_H

#include <array>
#include <cmath>
#include <tuple>

#include "std.hpp"

// Computes a table of the log scores of observing one geotype, given that
// another genotype if the true one; for three observed type (+missing
// observation) and two true states
template<class T>
static inline auto getLogScores(std::array<T, 4> const &errorRates)
{
    auto const &fd = errorRates[0];
    auto const &ad1 = errorRates[1];
    auto const &ad2 = errorRates[2];
    auto const &cc = errorRates[3];

    Matrix<T, 4, 2> logScores;

    logScores[0][0] = std::log(1.0 - cc - fd); // observed 0, true 0
    logScores[0][1] = std::log(ad1); // observed 0, true 1
    logScores[1][0] = std::log(fd); // observed 1, true 0
    logScores[1][1] = std::log(1.0 - (ad1 + ad2)); // observed 1, true 1
    logScores[2][0] = cc == 0.0 ? 0.0 : std::log(cc); // observed 2, true 0
    logScores[2][1] = ad2 == 0.0 ? 0.0 : std::log(ad2); // observed 2, true 1
    logScores[3][0] = std::log(1.0); // value N/A, true 0
    logScores[3][1] = std::log(1.0); // value N/A, true 1

    return logScores;
}

// Updates log scores after a new AD rate was accepted
template<class T>
static inline auto updateLogScores(Matrix<T, 4, 2> &logScores, T const AD)
{
    auto const [AD1, AD2] = [&]() {
        if(logScores[2][1] ==
            0) // The default case: There are not homozygous mutations observed
        {
            return std::tuple<T, T>(AD, 0);
        }
        else // Homozygous mutations were called. For simplicity set both
             // dropout rates to 1/2
        {
            return std::tuple<T, T>(AD / 2, AD / 2);
        }
    }();

    logScores[0][1] = std::log(AD1); // observed 0, true 1
    logScores[1][1] = std::log(T{1} - AD); // observed 1, true 1
    logScores[2][1] = AD == 0 ? 0.0 : std::log(AD2); // observed 2, true 1
    logScores[3][1] = std::log(1.0);
}

// computes an approximate scoring for a tree using the max attachment score per
// sample
template<class T, class F, std::size_t M, std::size_t S>
static inline auto maxScoreTreeFast(Matrix<F, 4, 2> const &logScores,
    Matrix<T, S, M> const &data, Vector<T, M> const &parent,
    Vector<T, M + 1> const &bft)
{
    F treeScore = 0.0;

    // sum over the best attachment scores of all samples is tree score
    for(std::size_t s = 0; s < S; ++s)
    {
        auto scores = getAttachmentScoresFast(parent, logScores, data[s], bft);
        treeScore += max(scores);
    }

    return treeScore;
}

template<class T, class F, std::size_t M, std::size_t S>
static inline auto maxScoreTreeAccurate(Matrix<F, 4, 2> const &logScores,
    Matrix<T, S, M> const &data, Vector<T, M> const &parent,
    Vector<T, M + 1> const &bft)
{
    Matrix<T, 4, 2> treeScores;
    treeScores.fill(0);

    for(std::size_t s = 0; s < S; ++s)
    {
        auto bestAttachment =
            getBestAttachmentScoreAccurate(parent, logScores, data[s], bft);
        treeScores += bestAttachment;
    }

    return getTrueScore(treeScores, logScores);
}

// computes an approximate scoring for a tree summing the score over all
// attachment points per sample
template<class T, class F, std::size_t M, std::size_t S>
static inline auto sumScoreTreeFast(Matrix<F, 4, 2> const &logScores,
    Matrix<T, S, M> const &data, Vector<T, M> const &tree,
    Vector<T, M + 1> const &bft)
{
    auto score = F{0};

    for(std::size_t s = 0; s < S; ++s)
    {
        auto scores = getAttachmentScoresFast(tree, logScores, data[s], bft);
        auto maxScore = max(scores);

        auto sum = F{0};
        for(std::size_t m = 0; m <= M; ++m)
        {
            sum += std::exp(scores[bft[m]] - maxScore);
        }
        score += std::log(sum) + maxScore;
    }

    return score;
}

template<class T, class F, std::size_t M>
static inline auto getSumAttachmentScoreAccurate(Vector<T, M> const &tree,
    Matrix<F, 4, 2> const &logScores, Vector<T, M> const &data,
    Vector<T, M + 1> const &bft)
{
    auto attachmentScores = getAttachmentMatrices(tree, data, bft);

    auto trueScores = getTrueScores(attachmentScores, logScores);
    auto maxScore = max(trueScores);

    auto sum = F{0};
    for(std::size_t n = 0; n < M + 1; ++n)
    {
        sum += std::exp(trueScores[n] - maxScore);
    }

    return std::log(sum) + maxScore;
}

// Computes te log score for the complete tree using the sumScore scheme where
// likelihoods of all attachment points of a sample are added
template<class T, class F, std::size_t M, std::size_t S>
static inline auto sumScoreTreeAccurate(Matrix<F, 4, 2> const &logScores,
    Matrix<T, S, M> const &data, Vector<T, M> const &tree,
    Vector<T, M + 1> const &bft)
{
    auto score = F{0};

    for(std::size_t s = 0; s < S; ++s)
    {
        score += getSumAttachmentScoreAccurate(tree, logScores, data[s], bft);
    }

    return score;
}

// Computes an approximate score for a tree. This is fast, but rounding errors
// can occur
template<char ScoreType, class T, class F, std::size_t M, std::size_t S>
static inline auto scoreTreeFast(Matrix<F, 4, 2> const &logScores,
    Matrix<T, S, M> const &data, Vector<T, M> const &parent)
{
    auto result = std::numeric_limits<F>::lowest();

    auto bft = getBreadthFirstTraversal(parent);

    if constexpr(ScoreType == 'm')
    {
        return maxScoreTreeFast(logScores, data, parent, bft);
    }
    else if constexpr(ScoreType == 's')
    {
        return sumScoreTreeFast(logScores, data, parent, bft);
    }
    else
    {
        static_assert(ScoreType != ScoreType, "Unsupported score Type");
    }
}

// Accurate score computation (minimizes rounding errors)
template<char ScoreType, class T, class F, std::size_t M, std::size_t S>
static inline auto scoreTreeAccurate(Matrix<F, 4, 2> const &logScores,
    Matrix<T, S, M> const &data, Vector<T, M> const &parent)
{
    auto result = std::numeric_limits<F>::lowest();

    auto bft = getBreadthFirstTraversal(parent);

    if constexpr(ScoreType == 'm')
    {
        return maxScoreTreeAccurate(logScores, data, parent, bft);
    }
    else if constexpr(ScoreType == 's')
    {
        return sumScoreTreeAccurate(logScores, data, parent, bft);
    }
    else
    {
        static_assert(ScoreType != ScoreType, "Unsupported score Type");
    }
}

// Compute the log likelihood for a single mutation for all subtrees of the
// binary tree where the expected state of the mutation can either be absent or
// present in the while subtree
template<bool State, class T, class F, std::size_t M, std::size_t S,
    std::size_t C>
static inline auto getBinarySubtreeScore(Matrix<T, S, M> const &data,
    Matrix<F, 4, 2> const &logScores,
    std::vector<std::vector<T>> const &children, Vector<T, C> const &bft,
    std::size_t const mutation)
{
    Vector<F, C> score;
    // Required ?
    score.fill(0.0);

    for(std::size_t n = C - 1; n < C; --n)
    {
        auto node = bft[n];

        if(node < S)
        {
            // For leafs the score is P(Dij|Eij)
            score[node] = logScores[data[node][mutation]][State];
        }
        else
        {
            // Is tree actually binary?
            if(children[node].size() != 2)
            {
                std::cerr << "Error: Node " << node << " has "
                          << children[node].size() << " children." << std::endl;
                std::abort();
            }

            score[node] = score[children[node][0]] + score[children[node][1]];
        }
    }
    return score;
}

// Compute the best log likelihood for placing a single mutation in a given
// sample tree. Iterates through all ndoes as possible placements of the
// mutation to find the best one. All samples below the placement of the
// mutation should have it, mutation can also be placed at a leaf, i.e. uniquely
// to the sample
template<class T, class F, std::size_t M, std::size_t S, std::size_t C>
static inline auto getBinaryTreeMutationScore(Matrix<T, S, M> const &data,
    Matrix<F, 4, 2> const &logScores,
    std::vector<std::vector<T>> const &children, Vector<T, C> const &bft,
    std::size_t const mutation)
{
    auto bestScore = std::numeric_limits<F>::lowest();
    auto absentScore =
        getBinarySubtreeScore<false>(data, logScores, children, bft, mutation);
    auto presentScore =
        getBinarySubtreeScore<true>(data, logScores, children, bft, mutation);

    for(std::size_t n = 0; n < C; ++n)
    {
        auto score = absentScore[C - 1] - absentScore[n] + presentScore[n];
        bestScore = std::max(bestScore, score);
    }

    return bestScore;
}

// Compute the maximum log likelihood of a binary tree for a given mutation
// matrix.
template<class T, class F, std::size_t M, std::size_t S, std::size_t C>
static inline auto getBinaryTreeScore(Matrix<T, S, M> const &data,
    Matrix<F, 4, 2> const &logScores, Vector<T, C> const &tree)
{
    auto children = getChildren(tree);
    auto bft = getBreadthFirstTraversal(tree);

    auto sum = F{0};
    for(std::size_t m = 0; m < M; ++m)
    {
        sum += getBinaryTreeMutationScore(data, logScores, children, bft, m);
    }

    return sum;
}

// Computes the score of a new candidate tree. First a fast approximate score is
// computed, then if the new score is better, or slightly worse than the best
// score so far, a more accurate but more costly score computation is done.
template<char ScoreType, class T, class F, std::size_t M, std::size_t S>
static inline auto scoreTree(Matrix<F, 4, 2> const &logScores,
    Matrix<T, S, M> const &data, Vector<T, M> const &parent,
    F const bestTreeLogScore)
{
    constexpr F epsilon = 0.000000000001;
    auto approx = scoreTreeFast<ScoreType>(logScores, data, parent);

    if(approx > bestTreeLogScore - epsilon)
    {
        return scoreTreeAccurate<ScoreType>(logScores, data, parent);
    }

    return approx;
}

template<char TreeType, char ScoreType, class T, class F, std::size_t M,
    std::size_t S, std::size_t N>
static inline auto getTreeLogScore(Matrix<T, S, M> const &data,
    Vector<T, N> const &tree, Matrix<F, 4, 2> const &logScores)
{
    if constexpr(TreeType == 'm')
    {
        return scoreTreeAccurate<ScoreType>(logScores, data, tree);
    }
    else if constexpr(TreeType == 't')
    {
        return getBinaryTreeScore(data, logScores, tree);
    }
    else
    {
        static_assert(
            TreeType != TreeType, "Unknown TreeType in getTreeLogScore.");
    }
}

template<char TreeType, char ScoreType, class T, class F, std::size_t M,
    std::size_t S, std::size_t N>
static inline auto getTreeLogScore(Matrix<T, S, M> const &data,
    Vector<T, N> const &tree, Matrix<F, 4, 2> const &logScores,
    F const bestTreeLogScore)
{
    if constexpr(TreeType == 'm')
    {
        return scoreTree<ScoreType>(logScores, data, tree, bestTreeLogScore);
    }
    else if constexpr(TreeType == 't')
    {
        return getBinaryTreeScore(data, logScores, tree);
    }
    else
    {
        static_assert(
            TreeType != TreeType, "Unknown TreeType in getTreeLogScore.");
    }
}

// computes the log score for attaching a sample to the root node (this score is
// equal for all trees)
template<class T, class F, std::size_t M>
static inline auto rootAttachmentScore(
    Matrix<F, 4, 2> const &logScores, Vector<T, M> const &data)
{
    F score = 0.0;

    for(std::size_t m = 0; m < M;
        ++m) // sum over log scores for all other nodes in tree
    {
        score += logScores[data[m]][0]; // none of them is ancestor of the
                                        // sample as it is attached to root
    }

    return score;
}

// computes the attachment scores of a sample to all nodes in the tree (except
// root)
template<class T, class F, std::size_t M>
static inline auto getAttachmentScoresFast(Vector<T, M> const &parent,
    Matrix<F, 4, 2> const &logScores, Vector<T, M> const &data,
    Vector<T, M + 1> const &bft)
{
    Vector<F, M + 1> attachmentScores;
    attachmentScores.fill(std::numeric_limits<F>::lowest());

    attachmentScores[M] = rootAttachmentScore(logScores, data);

    for(std::size_t m = 1; m <= M; ++m)
    {
        auto node = bft[m];
        attachmentScores[node] = attachmentScores[parent[node]];
        attachmentScores[node] -= logScores[data[node]][0];
        attachmentScores[node] += logScores[data[node]][1];
    }

    return attachmentScores;
}

template<class T, class F, std::size_t M>
static inline auto getBestAttachmentScoreAccurate(Vector<T, M> const &parent,
    Matrix<F, 4, 2> const &logScores, Vector<T, M> const &data,
    Vector<T, M + 1> const &bft)
{
    auto attachmentScores = getAttachmentMatrices(parent, data, bft);
    auto bestScore = std::numeric_limits<F>::lowest();
    std::size_t bestMatrix = 0;

    for(std::size_t m = 0; m <= M; ++m)
    {
        // now get true attachment scores and find best score among all
        // attachment points
        auto newScore = getTrueScore(attachmentScores[m], logScores);
        if(newScore > bestScore)
        {
            bestMatrix = m;
            bestScore = newScore;
        }
    }

    return attachmentScores[bestMatrix];
}

/* computes the attachment scores of a sample to all nodes in a tree, score is a
 * matrix counting the number of different match/mismatch score types */
template<class T, std::size_t M>
static inline auto getAttachmentMatrices(Vector<T, M> const &parent,
    Vector<T, M> const &data, Vector<T, M + 1> const &bft)
{
    Vector<Matrix<T, 4, 2>, M + 1> attachmentScores;

    // start with attaching node to root (no genes mutated)
    attachmentScores[M].fill(0);

    for(std::size_t m = 0; m < M; ++m)
    {
        ++attachmentScores[M][data[m]][0];
    }

    // now get scores for the other nodes due to bft traversal in an order such
    // that attachment matrix of parent is already filled
    for(std::size_t i = 1; i <= M; ++i)
    {
        auto node = bft[i];
        attachmentScores[node] = attachmentScores[parent[node]];
        --attachmentScores[node][data[node]][0];
        ++attachmentScores[node][data[node]][1];
    }
    return attachmentScores;
}

// computes the attachment score of a sample to a tree from a matrix
// representation of the score (counting the # 0->1, 1->0, 0->0, ...)
template<class T, class F>
static inline auto getTrueScore(
    Matrix<T, 4, 2> const &attachmentScore, Matrix<F, 4, 2> const &logScores)
{
    F score = 0.0;
    for(std::size_t j = 0; j < 4; ++j)
    {
        for(std::size_t k = 0; k < 2; ++k)
        {
            score += attachmentScore[j][k] * logScores[j][k];
        }
    }
    return score;
}

template<class T, class F, std::size_t N>
static inline auto getTrueScores(
    Vector<Matrix<T, 4, 2>, N> const &attachmentScores,
    Matrix<F, 4, 2> const &logScores)
{
    Vector<F, N> scores;
    for(std::size_t n = 0; n < N; ++n)
    {
        scores[n] = getTrueScore(attachmentScores[n], logScores);
    }
    return scores;
}
#endif // SCORE_TREE_H
