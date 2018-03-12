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

/* given a Pruefer code, compute the corresponding parent vector */
template<class T>
auto prueferCode2parentVector(DynamicArray<T> &code)
{
    auto nodeCount = code.size() + 1;
    DynamicArray<T> parent(nodeCount);
	int* lastOcc = getLastOcc(code.data(), code.size());    // node id -> index of last occ in code, -1 if no occurrence or if id=root
	bool* queue = getInitialQueue(code.data(), code.size());  // queue[node]=true if all children have been attached to this node, or if it is leaf
	int queueCutter = -1;    // this is used for a node that has been passed by the "queue" before all children have been attached
	int next = getNextInQueue(queue, 0, code.size()+1);

	for(int i=0; i<code.size(); i++){               // add new edge to tree from smallest node with all children attached to its parent
		if(queueCutter >=0){
			parent[queueCutter] = code[i];         // this node is queueCutter if the queue has already passed this node
			//cout << queueCutter << " -> " << code[i] << "\n";
			queueCutter = -1;
		}
		else{
			parent[next] = code[i];                               // use the next smallest node in the queue, otherwise
			//cout << next << " -> " << code[i] << "\n";
			next = getNextInQueue(queue, next+1, code.size()+1);     // find next smallest element in the queue
		}

		if(lastOcc[code[i]]==i){                               // an element is added to the queue, or we have a new queueCutter
			updateQueue(code[i], queue, next);
			queueCutter = updateQueueCutter(code[i], queue, next);
		}
	}
	if(queueCutter>=0){
		parent[queueCutter] = nodeCount;
		//cout << queueCutter << " -> " << nodeCount << "\n";
	}
	else{
		parent[next] = nodeCount;
		//cout << next << " -> " << nodeCount << "\n";
	}

	delete [] lastOcc;
	delete [] queue;
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
