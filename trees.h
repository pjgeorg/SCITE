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
int countBranches(int* parents, int length);
void deleteChildLists(vector<vector<int> > &childLists);
string getNewickCode(vector<vector<int> > list, int root);
int* prueferCode2parentVector(int* code, int codeLength);
int* getBreadthFirstTraversal(int* parent, int n);
bool** parentVector2ancMatrix(int* parent, int n);
int* getRandParentVec(int n);
bool* getInitialQueue(int* code, int codeLength);
int* getLastOcc(int* code, int codeLength);
int getNextInQueue(bool* queue, int pos, int length);
void updateQueue(int node, bool* queue, int next);
int updateQueueCutter(int node, bool* queue, int next);
int* starTreeVec(int n);
bool** starTreeMatrix(int n);
int* reverse(int* array, int length);

#include "matrices.h"

// TODO: Replate input by DynamicArray
/* converts a parent vector to the list of children */
template<class T>
auto getChildListFromParentVector(T *parents, int n)
{
    DynamicArray<DynamicArray<T>> childList(n+1);
    for(std::size_t i = 0; i < n; ++i)
    {
        childList[parents[i]].push_back(i);
    }
    return childList;
}


#endif /* TREES_H_ */
