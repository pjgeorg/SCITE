/*
 * treelist.hpp
 *
 *  Created on: Mar 27, 2015
 *      Author: jahnka
 *  Modified on: Sep, 2018
 *      Author: pjgeorg
 */

#ifndef TREELIST_H
#define TREELIST_H

template<class T, class F, std::size_t M>
static inline auto resetTreeList(std::vector<std::pair<Vector<T, M>, F>> &trees,
    Vector<T, M> const &parent, F const beta)
{
    trees.clear();
    trees.emplace_back(parent, beta);
}

template<class T, class F, std::size_t M>
static inline auto addToTreeList(std::vector<std::pair<Vector<T, M>, F>> &trees,
    Vector<T, M> const &parent, F const beta)
{
    if(!treeListContains(trees, parent))
    {
        trees.emplace_back(parent, beta);
    }
}

template<class T, class F, std::size_t M>
static inline auto treeListContains(
    std::vector<std::pair<Vector<T, M>, F>> &trees, Vector<T, M> const &newTree)
{
    for(auto const &tree : trees)
    {
        if(tree.first == newTree)
        {
            return true;
        }
    }
    return false;
}
#endif // TREELIST_H
