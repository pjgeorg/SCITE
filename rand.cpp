/*
 * rand_C.cpp
 *
 *  Created on: Mar 27, 2015
 *      Author: jahnka
 */

//#include <string>
//#include <random>
#include "rand.h"
#include <iostream>
#include <random>
#include <stdlib.h>
#include <time.h>
#include "matrices.h"

std::default_random_engine RandomGenerator::sEngine((std::random_device())());

using namespace std;


/* This function gets a number of nodes n, and creates a random pruefer code for a rooted tree with n+1 nodes (root is always node n+1) */
int* getRandTreeCode(int n){                // as usual n is the number of mutations

	int nodes = n+1;                        // #nodes = n mutations plus root (wildtype)
	int codeLength = nodes-2;
	int* code = new int[codeLength];
	for(int i=0; i<codeLength; i++){
		code[i] = rand() % nodes;
	}
	return code;
}

bool changeBeta(double prob){
	 double percent = (rand() % 100)+1;    // between 1 and 100
	 if(percent <= prob*100){
		 return true;
	 }
	 return false;
}

int sampleRandomMove(std::vector<double> prob){ // picks randomly one of the tree moves based on the move probabilities

    double percent = (rand() % 100)+1;    // between 1 and 100
    double probSum = prob[1];
    for(int i=1; i<prob.size()-1; i++){    // start at index 1; the probability at prob[0] is for changing the error rate (which is treated separately)
        if(percent <= probSum*100){
          return i;
        }
        probSum += prob[i+1];
    }
    return prob.size()-1;
}


bool samplingByProb(double prob){
	double percent = rand() % 100;
	if(percent <= prob*100){
		return true;
	}
	return false;
}


int* sampleTwoElementsWithoutReplacement(int n){

    int* result = new int[2];
	  result[0] = rand() % n;
	  result[1] = result[0];
    while(result[0]==result[1]){
      result[1] = rand() % n;
    }
	  return result;
}

int pickRandomNumber(int n){

    return getRandomNumber(0, n-1);
}

double sample_0_1(){

  //return (((double) rand()+0.5) / ((RAND_MAX+1)));
  return ((double) rand() / RAND_MAX);
}

int getElemFromQueue(int index, std::vector<int> queue){
	int elem = queue.at(index);
	if (index != queue.size() - 1)
	{
		queue[index] = std::move(queue.back());
	}

	//cout << queue.size() << " elements in queue in subroutine\n";
	return elem;
}
