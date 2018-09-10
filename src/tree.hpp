/*
 * tree.hpp
 *
 *  Created on: Oct 12, 2015
 *      Author: jahnka
 *  Modified on: Sep, 2018
 *      Author: pjgeorg
 */

#ifndef TREE_H
#define TREE_H

#include "matrices.hpp"
#include "rand.hpp"

template<class T, std::size_t M>
static inline auto getLastOcc(Vector<T, M> const &code)
{
    Vector<T, M + 2> lastOcc;

    lastOcc.fill(-1);

    constexpr auto root = M + 1;
    for(std::size_t i = 0; i < M; ++i)
    {
        if(code[i] != root)
        {
            lastOcc[code[i]] = i;
        }
    }
    return lastOcc;
}

template<class T, std::size_t M>
static inline auto getInitialQueue(Vector<T, M> const &code)
{
    Vector<bool, M + 2> queue;
    queue.fill(true);

    for(std::size_t i = 0; i < M; ++i)
    {
        queue[code[i]] = false;
    }

    return queue;
}

template<std::size_t M>
static inline auto getNextInQueue(
    Vector<bool, M> const &queue, std::size_t const pos)
{
    for(auto i = pos; i < M - 1; ++i)
    {
        if(queue[i] == true)
        {
            return i;
        }
    }

    // No node left in queue. Possibly a cycle?
    return M - 1;
}

template<std::size_t M>
static inline auto updateQueue(
    std::size_t const node, Vector<bool, M> &queue, std::size_t const next)
{
    if(node >= next)
    {
        // add new node to queue
        queue[node] = true;
    }
}

static inline auto updateQueueCutter(
    std::size_t const node, std::size_t const next)
{
    if(node >= next)
    {
        // new node can be added to the queue
        return std::pair<bool, std::size_t>(false, 0);
    }
    else
    {
        // new node needs to cut the queue, as it has already passed it
        return std::pair<bool, std::size_t>(true, node);
    }
}

// Get ancestors of a parent tree
template<class T, std::size_t M>
static inline auto getAncestors(Vector<T, M> const &tree)
{
    constexpr auto root = M;

    Matrix<bool, M, M> ancestors;

    ancestors.fill(false);

    for(std::size_t i = 0; i < M; ++i)
    {
        auto anc = i;
        std::size_t its = 0;
        while(anc < root)
        {
            // if the ancestor is the root node, it is not represented in the
            // adjacency matrix
            if(tree[anc] < root)
            {
                ancestors[tree[anc]][i] = true;
            }

            anc = tree[anc];
            ++its;
        }
    }

    for(std::size_t i = 0; i < M; ++i)
    {
        ancestors[i][i] = true;
    }

    return ancestors;
}

template<class T, std::size_t N>
static inline auto getChildren(Vector<T, N> const &tree)
{
    std::vector<std::vector<T>> children(N + 1);
    for(std::size_t i = 0; i < N; ++i)
    {
        children[tree[i]].push_back(i);
    }
    return children;
}

// Computes a breadth first traversal of a tree from the parent tree vector
template<class T, std::size_t N>
static inline auto getBreadthFirstTraversal(Vector<T, N> const &tree)
{
    auto childrens = getChildren(tree);

    Vector<T, N + 1> bft;
    bft[0] = N;

    std::size_t k = 0;
    for(std::size_t i = 0; i < N + 1; ++i)
    {
        for(std::size_t j = 0; j < childrens[bft[i]].size(); ++j)
        {
            bft[++k] = childrens[bft[i]][j];
        }
    }

    return bft;
}

// Returns all nodes that are descendants of the given node
// Note: ancestor is 1 at [i][j] if i is an ancestor of j in the tree
template<class T, std::size_t M>
static inline auto getDescendants(
    Matrix<bool, M, M> const &ancestors, std::size_t const node)
{
    std::vector<T> descendants;
    descendants.reserve(M);

    for(std::size_t i = 0; i < M; ++i)
    {
        if(ancestors[node][i] == true)
        {
            descendants.push_back(i);
        }
    }
    return descendants;
}

// Returns all nodes that are not descendants of the given node
// i.e. ancestors and nodes in a different branch of the tree
// note: ancestor is 0 at [i][j] if i is not an ancestor of j in the tree
template<class T, std::size_t M>
static inline auto getNonDescendants(
    Matrix<bool, M, M> const &ancestors, std::size_t const node)
{
    std::vector<T> nonDescendants;
    nonDescendants.reserve(M);

    for(std::size_t i = 0; i < M; ++i)
    {
        if(ancestors[node][i] == false)
        {
            nonDescendants.push_back(i);
        }
    }
    return nonDescendants;
}

// Counts the number of branches in a tree, this is equal to the number of leafs
// in the tree
template<class T, std::size_t M>
static inline auto countBranches(Vector<T, M> const &tree)
{
    auto count = T{0};

    for(auto const &child : getChildren(tree))
    {
        if(child.size() == 0)
        {
            ++count;
        }
    }

    return count;
}

// Returns the distance between two trees where the distance is the number of
// nodes having different parents in the two trees
template<class T, std::size_t M>
static inline auto getSimpleDistance(
    Vector<T, M> const &a, Vector<T, M> const &b)
{
    auto distance = T{0};
    for(std::size_t i = 0; i < M; ++i)
    {
        if(a[i] != b[i])
        {
            ++distance;
        }
    }

    return distance;
}

template<class T, class F, std::size_t M>
static inline auto getMinDistToTrueTree(
    std::vector<std::pair<Vector<T, M>, F>> const &trees,
    Vector<T, M> const &trueTree)
{
    auto distance = std::numeric_limits<T>::max();
    for(auto const &tree : trees)
    {
        distance = std::min(distance, getSimpleDistance(tree.first, trueTree));
    }

    return distance;
}

#endif // TREE_H
