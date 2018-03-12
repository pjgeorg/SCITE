/*
 * trees.h
 *
 *  Created on: Oct 12, 2015
 *      Author: jahnka
 */

#ifndef TREES_H_
#define TREES_H_
using namespace std;

std::vector<int> getDescendants(bool** ancMatrix, int node, int n);
std::vector<int> getNonDescendants(bool**& ancMatrix, int node, int n);
string getNewickCode(vector<vector<int> > list, int root);
int* prueferCode2parentVector(int* code, int codeLength);
int* getBreadthFirstTraversal(int* parent, int n);
bool** parentVector2ancMatrix(int* parent, int n);
bool* getInitialQueue(int* code, int codeLength);
int* getLastOcc(int* code, int codeLength);
int getNextInQueue(bool* queue, int pos, int length);
void updateQueue(int node, bool* queue, int next);
int updateQueueCutter(int node, bool* queue, int next);
int* starTreeVec(int n);
bool** starTreeMatrix(int n);

#include "matrices.h"
#include "rand.h"

/* converts a parent vector to the list of children */
template<class T>
auto getChildListFromParentVector(DynamicArray<T> const &parents)
{
    DynamicArray<DynamicArray<T>> childList(parents.size()+1);
    for(std::size_t i = 0; i < parents.size(); ++i)
    {
        childList[parents[i]].push_back(i);
    }
    return childList;
}

/* counts the number of branches in a tree, this is the same as the number of leafs in the tree */
template<class T>
auto countBranches(DynamicArray<T> const & parents)
{
    std::size_t count = 0;
	auto childList = getChildListFromParentVector(parents);
	for(std::size_t i=0; i<childList.size(); ++i)
    {
		if(childList[i].size()==0)
        {
            ++count;
        }
	}
	return count;
}

/* creates a random parent vector for nodes 0, .., n with node n as root*/
template<class T>
auto getRandParentVec(std::size_t n)
{
	auto randCode = getRandTreeCode(n);
    return prueferCode2parentVector(randCode.data(), n-1);
}

// TODO: Temporary functions, remove after refactor
template<class T>
auto getChildListFromParentVector(T *parents, int length)
{
    auto tParents = toDynamicArray(parents, length);
    return getChildListFromParentVector(tParents);
}

template<class T>
auto countBranches(T* parents, int length)
{
    auto tParents = toDynamicArray(parents, length);
    return countBranches(tParents);
}
#endif /* TREES_H_ */
