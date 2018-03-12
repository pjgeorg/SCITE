/*
 * trees.cpp
 *
 *  Created on: Oct 12, 2015
 *      Author: jahnka
 */

#include <stdbool.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <sstream>
#include <math.h>
#include <queue>
#include "matrices.h"
#include "treelist.h"
#include "trees.h"
#include "output.h"

using namespace std;


/* returns all nodes that are descendants of the given node */
/* note: ancMatrix is 1 at [i,j] if i is an ancestor of j in the tree */
std::vector<int> getDescendants(bool** ancMatrix, int node, int n){
  std::vector<int> descendants;
  for(int i=0; i<n; i++){
  	if(ancMatrix[node][i]==true){
			descendants.push_back(i);
		}
	}
	return descendants;
}

/* returns all nodes that are not descendants of the given node */
/* i.e. ancestors and nodes in a different branch of the tree   */
/* note: ancMatrix is 0 at [i,j] if i is not an ancestor of j in the tree */
std::vector<int> getNonDescendants(bool**& ancMatrix, int node, int n){
	std::vector<int> ancestors;
	for(int i=0; i<n; i++){
		if(ancMatrix[node][i]==false){
			ancestors.push_back(i);
		}
	}
	return ancestors;
}

/* converts a tree given as lists of children to the Newick tree format */
/* Note: This works only if the recursion is started with the root node which is n+1 */
string getNewickCode(vector<vector<int> > list, int root){
	stringstream newick;
	vector<int> rootChilds = list.at(root);
	if(!rootChilds.empty()){
		newick << "(";
		bool first = true;
		for(int i=0; i<rootChilds.size(); i++){
			if(!first){
				newick << ",";
			}
			first = false;
			newick << getNewickCode(list, rootChilds.at(i));
		}
		newick << ")";
	}
	newick << root+1;
	return newick.str();
}



/*  computes a breadth first traversal of a tree from the parent vector  */
int* getBreadthFirstTraversal(int* parent, int n){

	vector<vector<int> > childLists = getChildListFromParentVector(parent, n);
	int* bft = new int[n+1];
	bft[0] = n;
	int k = 1;

	for(int i=0; i<n+1; i++){
		for(int j=0; j<childLists[bft[i]].size(); j++){
			bft[k++] = childLists[bft[i]][j];
		}
	}
	for(int i=0; i<childLists.size(); i++){
		childLists[i].clear();
	}
	childLists.clear();
	return bft;
}

/* transforms a parent vector to an ancestor matrix*/
bool** parentVector2ancMatrix(int* parent, int n){
	bool** ancMatrix = init_boolMatrix(n, n, false);
	int root = n;
	for(int i=0; i<n; i++){
		int anc = i;
		int its =0;
		while(anc < root){                              // if the ancestor is the root node, it is not represented in the adjacency matrix
			if(parent[anc]<n){
				ancMatrix[parent[anc]][i] = true;
			}

			anc = parent[anc];
			its++;
		}
	}
	for(int i=0; i<n; i++){
		ancMatrix[i][i] = true;
	}
	return ancMatrix;
}

bool* getInitialQueue(int* code, int codeLength){
	//cout << "code Length: " << codeLength << "\n";
	int queueLength = codeLength+2;
	//cout << "queueLength: " << queueLength << "\n";
	bool* queue = init_boolArray(queueLength, true);

	for(int i=0; i<codeLength; i++){
		queue[code[i]] = false;
	}
	return queue;
}


void updateQueue(int node, bool* queue, int next){

	if(node>=next){                //  add new node to queue
		queue[node] = true;
	}
}

int updateQueueCutter(int node, bool* queue, int next){
	if(node>=next){
		return -1;         // new node can be added to the queue
	}
	else{
		return node;         // new node needs to cut the queue, as it has already passed it
	}
}


int* getLastOcc(int* code, int codeLength){
	int* lastOcc = init_intArray(codeLength+2, -1);
	int root = codeLength+1;
	for(int i=0; i<codeLength; i++){
		if(code[i] != root){
			lastOcc[code[i]] = i;
		}
	}
	return lastOcc;
}

int getNextInQueue(bool* queue, int pos, int length){
	for(int i=pos; i<length; i++){
		if(queue[i]==true){
			return i;
		}
	}
	//cout << "No node left in queue. Possibly a cycle?";
	return length;
}

/* creates the parent vector for a star tree with node n as center and 0,...,n-1 as leafs */
int* starTreeVec(int n){
	int* starTreeVec = new int[n];
	for(int i=0;i<n;i++){
		starTreeVec[i] = n;
	}
	return starTreeVec;
}

/* creates the ancestor matrix for the same tree */
bool** starTreeMatrix(int n){

  bool** starTreeMatrix = init_boolMatrix(n, n, false);
  for(int i=0;i<n;i++){
		starTreeMatrix[i][i] = true;
	}
	return starTreeMatrix;
}
