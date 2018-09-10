/*
 * move_tree.h
 *
 *  Created on: Mar 15, 2016
 *      Author: jahnka
 *  Modified on: Sep, 2018
 *      Author: pjgeorg
 */

#ifndef MOVE_TREE_H
#define MOVE_TREE_H

#include <utility>
#include "rand_tree.hpp"
#include "tree.hpp"

// Creates a new parent tree vector after pruning and reattaching subtree
template<class T, std::size_t M>
static inline auto getNewTree(
    Vector<T, M> tree, T const nodeToMove, T const newParent)
{
    tree[nodeToMove] = newParent;
    return tree;
}

// Creates a new parent tree vector after swapping two nodes
template<class T, std::size_t M>
static inline auto getNewTree(Vector<T, M> tree, std::array<T, 2> const &nodes)
{
    auto &[first, second] = nodes;

    for(std::size_t m = 0; m < M; ++m)
    {
        if(tree[m] == first && m != second)
        {
            tree[m] = second;
        }
        else if(tree[m] == second && m != first)
        {
            tree[m] = first;
        }
    }

    std::swap(tree[first], tree[second]);

    if(tree[first] == first)
    {
        tree[first] = second;
    }
    if(tree[second] == second)
    {
        tree[second] = first;
    }

    return tree;
}

// Re-orders the nodes so that the descendant is first in case the nodes are in
// the same lineage
template<class T, std::size_t M>
static inline auto reorderToStartWithDescendant(
    std::array<T, 2> &nodes, Matrix<bool, M, M> const &ancestor)
{
    if(ancestor[nodes[0]][nodes[1]] == true)
    {
        std::swap(nodes[0], nodes[1]);
    }
}

template<class T, class F, std::size_t M, std::size_t P>
static inline auto proposeNewMutationTree(RNG &rng, Vector<T, M> const &tree,
    Matrix<bool, M, M> const &ancestors, std::array<F, P> const &moveProbs,
    F &correction)
{
    switch(sampleRandomMove(rng, moveProbs))
    {
            // Prune and re-attach
        case 1:
        {
            // Pick a node to move with its subtree
            auto nodeToMove = getRandomNumber<T>(rng, M);
            auto possibleParents = getNonDescendants<T>(ancestors, nodeToMove);
            auto parent = pickParent<T>(rng, possibleParents, M);
            return getNewTree(tree, nodeToMove, parent);
        }
            // Swap two node labels
        case 2:
        {
            auto nodesToSwap = sampleTwoElementsWithoutReplacement<T>(rng, M);
            return getNewTree(tree, nodesToSwap);
        }
            // Swap two subtress
        case 3:
        {
            auto nodesToSwap = sampleTwoElementsWithoutReplacement<T>(rng, M);
            // Make sure we move the descendant first (in case nodes are in same
            // lineage)
            reorderToStartWithDescendant(nodesToSwap, ancestors);

            auto const &[nodeToMove, nextNodeToMove] = nodesToSwap;

            if(ancestors[nextNodeToMove][nodeToMove] ==
                0) // Nodes are in different lineages
            {
                auto propTree = tree;
                std::swap(propTree[nodeToMove], propTree[nextNodeToMove]);

                return propTree;
            }
            else // Nodes are in the same lineage -- need to avoid cycle in the
                 // tree
            {
                auto propTree = tree;
                propTree[nodeToMove] = tree[nextNodeToMove];
                auto descendants = getDescendants<T>(ancestors, nodeToMove);
                auto propAncestors = getAncestors(propTree);
                auto nextDescendants =
                    getDescendants<T>(propAncestors, nextNodeToMove);

                propTree[nextNodeToMove] =
                    descendants[getRandomNumber(rng, descendants.size())];

                // Neighborhood correction needed for MCMC convergence, but not
                // important for simulated annealing
                correction =
                    static_cast<F>(descendants.size()) / nextDescendants.size();

                return propTree;
            }
        }
        default:
            std::cerr << "Unknown move type." << std::endl;
            std::abort();
    }
}

template<class T, std::size_t C>
static inline auto pickNodeToMove(RNG &rng, Vector<T, C> const &tree)
{
    while(true)
    {
        if(auto picked = getRandomNumber<T>(rng, C); tree[picked] != C)
        {
            return picked;
        }
    }
}

template<class T, std::size_t C>
static inline auto getSibling(Vector<T, C> const &tree,
    std::vector<std::vector<T>> const &children, T const nodeToMove)
{
    if(children[tree[nodeToMove]][0] != nodeToMove)
    {
        return children[tree[nodeToMove]][0];
    }
    else
    {
        return children[tree[nodeToMove]][1];
    }
}

template<class T, class F, std::size_t C, std::size_t P>
static inline auto proposeNewBinaryTree(RNG &rng, Vector<T, C> const &tree,
    Matrix<bool, C, C> const &ancestors, std::array<F, P> const &moveProbs)
{
    switch(sampleRandomMove(rng, moveProbs))
    {
            // Prune and re-attach
        case 1:
        {
            auto nodeToMove = pickNodeToMove(rng, tree);
            auto parent = tree[nodeToMove];
            auto children = getChildren(tree);
            auto sibling = getSibling(tree, children, nodeToMove);

            auto propTree = tree;
            propTree[sibling] = tree[parent];

            auto possibleSiblings = getNonDescendants<T>(ancestors, parent);

            if(possibleSiblings.size() == 0)
            {
                std::cerr << "Error: No new siblings found for node "
                          << nodeToMove << " for move type 1 in binary tree."
                          << std::endl;
                std::abort();
            }

            auto newSibling =
                possibleSiblings[getRandomNumber(rng, possibleSiblings)];
            propTree[newSibling] = parent;
            propTree[parent] = tree[newSibling];

            return propTree;
        }
            // Swap two node labels
        case 2:
        {
            constexpr auto mutations = (C + 2) / 2;
            auto first = getRandomNumber(rng, mutations);
            auto second = getRandomNumber(rng, mutations);

            auto propTree = tree;
            std::swap(propTree[first], propTree[second]);
            return propTree;
        }
        default:
            std::cerr << "Unknown move type." << std::endl;
            std::abort();
    }
}

template<char TreeType, class T, class F, std::size_t M, std::size_t P>
static inline auto proposeNewTree(RNG &rng, Vector<T, M> const &tree,
    Matrix<bool, M, M> const &ancestors, std::array<F, P> const &moveProbs,
    F &correction)
{
    if constexpr(TreeType == 'm')
    {
        return proposeNewMutationTree(
            rng, tree, ancestors, moveProbs, correction);
    }
    else if constexpr(TreeType == 't')
    {
        return proposeNewBinaryTree(rng, tree, ancestors, moveProbs);
    }
    else
    {
        static_assert(
            TreeType != TreeType, "Unknown TreeType in proposeNewTree.");
    }
}
#endif // MOVE_TREE_H
