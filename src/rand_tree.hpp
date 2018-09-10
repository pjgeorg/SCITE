/*
 * rand_tree.h
 *
 *  Created on: Mar 27, 2015
 *      Author: jahnka
 *  Modified on: Sep, 2018
 *      Author: pjgeorg
 */

#ifndef RAND_TREE_H
#define RAND_TREE_H

#include "rand.hpp"
#include "tree.hpp"

// creates a random parent tree vector for nodes 0, .., M with node M as root
template<class T, std::size_t M>
auto getRandomMutationTree(RNG &rng)
{
    // #nodes = mutations plus root (wildtype) = M + 1
    constexpr auto nodes = M + 1;

    // #codeLength = #nodes-2 = M - 1
    constexpr auto codeLength = M - 1;

    // #nodeCount = #codeLength+1 = M
    constexpr auto nodeCount = M;

    Vector<T, codeLength> code;

    for(std::size_t i = 0; i < codeLength; ++i)
    {
        code[i] = getRandomNumber<T>(rng, nodes);
    }

    Vector<T, nodeCount> parent;

    // node id -> index of last occ in code, -1 if no occurrence or if id=root
    auto lastOcc = getLastOcc(code);

    // queue[node]=true if all children have been attached to this node, or if
    // it is leaf
    auto queue = getInitialQueue(code);

    // this is used for a node that has been passed by the "queue" before all
    // children have been attached
    std::pair<bool, std::size_t> queueCutter(false, 0);
    auto next = getNextInQueue(queue, 0);

    for(std::size_t i = 0; i < codeLength; ++i)
    {
        // add new edge to tree from smallest node with all children attached to
        // its parent
        if(queueCutter.first)
        {
            // this node is queueCutter if the queue has already
            parent[queueCutter.second] = code[i];
            // passed this node
            queueCutter.first = false;
        }
        else
        {
            // use the next smallest node in the queue, otherwise
            parent[next] = code[i];
            // find next smallest element in the queue
            next = getNextInQueue(queue, next + 1);
        }

        if(lastOcc[code[i]] == i)
        {
            // an element is added to the queue, or we have a new queueCutter
            updateQueue(code[i], queue, next);
            queueCutter = updateQueueCutter(code[i], next);
        }
    }
    if(queueCutter.first)
    {
        parent[queueCutter.second] = nodeCount;
    }
    else
    {
        parent[next] = nodeCount;
    }

    return parent;
}

// This creates the parent vector of a random binary tree. Entries 0...S-1 are
// for the leafs. Entries S....2S-3 are for the inner nodes except the root, the
// root has index 2S-2 which has no parent and therefore has no entry in the
// parent vector
template<class T, std::size_t N>
inline auto getRandomBinaryTree(RNG &rng)
{
    constexpr auto samples = (N + 2) / 2;

    Vector<T, N> leafsAndInnerNodesParents;

    auto count = samples;
    std::vector<T> queue(count);
    for(std::size_t i = 0; i < queue.size(); ++i)
    {
        queue[i] = i;
    }

    while(queue.size() > 1)
    {
        auto pos = getRandomNumber(rng, queue);
        auto child1 = queue[pos];
        queue[pos] = queue.back();
        queue.pop_back();

        pos = getRandomNumber(rng, queue);
        auto child2 = queue[pos];
        queue[pos] = count;

        leafsAndInnerNodesParents[child1] = count;
        leafsAndInnerNodesParents[child2] = count;
        ++count;
    }

    return leafsAndInnerNodesParents;
}

template<char TreeType, class T, std::size_t N>
static inline auto getRandomTree(RNG &rng)
{
    if constexpr(TreeType == 'm')
    {
        return getRandomMutationTree<T, N>(rng);
    }
    else if constexpr(TreeType == 't')
    {
        return getRandomBinaryTree<T, N>(rng);
    }
    else
    {
        static_assert(
            TreeType != TreeType, "Unknown TreeType in getRandomTree.");
    }
}

// Picks a parent ramdomly from the set of possible parents (including the root
// (M+1))
template<class T>
static inline auto pickParent(
    RNG &rng, std::vector<T> const &parents, T const root)
{
    if(auto idx = getRandomNumber(rng, parents.size() + 1);
        idx == parents.size())
    {
        return root;
    }
    else
    {
        return parents[idx];
    }
}

// picks randomly one of the tree moves based on the move probabilities
template<class T, std::size_t P>
inline auto sampleRandomMove(RNG &rng, std::array<T, P> const &prob)
{
    auto percent = getRandomNumber<T>(rng, 0, 1);
    auto probSum = prob[1];

    // start at index 1; the probability at prob[0] is for changing the error
    // rate (which is treated separately)
    for(std::size_t i = 1; i < P - 1; ++i)
    {
        if(percent <= probSum)
        {
            return i;
        }
        probSum += prob[i + 1];
    }
    return prob.size() - 1;
}
#endif // RAND_TREE_H
