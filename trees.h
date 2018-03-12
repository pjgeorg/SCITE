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
int getNextInQueue(bool* queue, int pos, int length);
void updateQueue(int node, bool* queue, int next);
int updateQueueCutter(int node, bool* queue, int next);
int* starTreeVec(int n);
bool** starTreeMatrix(int n);

#include <utility>
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

template<class T>
auto getLastOcc(DynamicArray<T> const &code)
{
    DynamicArray<T> lastOcc(code.size() + 2, -1);
    auto root = code.size() + 1;
	for(std::size_t i=0; i<code.size(); ++i)
    {
		if(code[i] != root)
        {
			lastOcc[code[i]] = i;
		}
	}
	return lastOcc;
}

template<class R, class T>
auto getInitialQueue(DynamicArray<T> const &code)
{
    DynamicArray<R> queue(code.size() + 2, true);

	for(std::size_t i=0; i<code.size(); ++i)
    {
		queue[code[i]] = false;
	}
	return queue;
}

template<class T>
auto getNextInQueue(DynamicArray<T> const&queue, std::size_t const pos)
{
    for(auto i = pos; i<queue.size() -1; ++i)
    {
		if(queue[i]==true)
        {
			return i;
		}
	}
	// No node left in queue. Possibly a cycle?";
	return queue.size() - 1;
}

template<class T>
auto updateQueue(std::size_t const node, DynamicArray<T> &queue, std::size_t const next)
{
	if(node>=next)
    {
        // add new node to queue
		queue[node] = true;
	}
}

static inline auto updateQueueCutter(std::size_t const node, std::size_t const next)
{
	if(node>=next)
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

/* given a Pruefer code, compute the corresponding parent vector */
template<class T>
auto prueferCode2parentVector(DynamicArray<T> &code)
{
    auto nodeCount = code.size() + 1;
    DynamicArray<T> parent(nodeCount);

    // node id -> index of last occ in code, -1 if no occurrence or if id=root
	auto lastOcc = getLastOcc(code);    

    // queue[node]=true if all children have been attached to this node, or if it is leaf
    auto queue = getInitialQueue<bool>(code);
  
// this is used for a node that has been passed by the "queue" before all children have been attached
    std::pair<bool, std::size_t> queueCutter(false, 0);    
	auto next = getNextInQueue(queue, 0);

	for(std::size_t i=0; i<code.size(); ++i){               // add new edge to tree from smallest node with all children attached to its parent
		if(queueCutter.first){
			parent[queueCutter.second] = code[i];         // this node is queueCutter if the queue has already passed this node
			queueCutter.first = false;
		}
		else{
			parent[next] = code[i];                               // use the next smallest node in the queue, otherwise
			next = getNextInQueue(queue, next+1);     // find next smallest element in the queue
		}

		if(lastOcc[code[i]]==i){                               // an element is added to the queue, or we have a new queueCutter
			updateQueue(code[i], queue, next);
			queueCutter = updateQueueCutter(code[i], next);
		}
	}
	if(queueCutter.first){
		parent[queueCutter.second] = nodeCount;
	}
	else{
		parent[next] = nodeCount;
	}

	return parent;
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
    return prueferCode2parentVector(randCode);
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
