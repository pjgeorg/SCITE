/*
 * findBestTrees_noR.cpp
 *
 *  Created on: Mar 27, 2015
 *      Author: jahnka
 *  Modified on: Aug, 2018
 *      Author: pjgeorg
 */

#include <stdbool.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include "matrices.h"
#include "treelist.h"
#include "trees.h"
#include "output.h"
#include "mcmc.h"
#include "rand.h"
#include "parameter.h"
#include "scoreTree.h"

template<std::size_t G, std::size_t S>
static inline auto getDataMatrix(std::string const &fileName)
{
    Matrix<int, S, G> dataMatrix;

    ifstream in(fileName.c_str());

    if (!in) {
        std::cerr << "Cannot open file " << fileName << std::endl;
        std::abort();
    }

    for (std::size_t g = 0; g < G; ++g) {
        for (std::size_t s = 0; s < S; ++s) {
            in >> dataMatrix[s][g];
        }
    }

    in.close();

    return dataMatrix;
}

template<std::size_t G>
static inline auto getTrueTree(std::string const &fileName)
{
    std::array<int, G> trueTree;

	std::vector<std::string> lines;
	std::ifstream input(fileName.c_str());

	std::string line;
	while(std::getline(input, line))
    {
	    if (!line.empty())
        {
	        lines.push_back(line);
        }
	}

	for(auto const &line : lines)
    {
        if(auto found = line.find(" -> "); found != std::string::npos)
        {
            auto parent = std::stoi(line.substr(0, found));
            auto child = std::stoi(line.substr(found + 3));
            trueTree[child - 1] = parent - 1;
        }
	}
	return trueTree;
}

template<char ScoreType>
static inline auto getOutputFileName(std::size_t i, std::string const &prefix, std::string const &extension)
{
	std::stringstream fileName;
    fileName << prefix;

    if constexpr(ScoreType == 'm')
    {
        fileName << "_ml";
	}
    else if constexpr(ScoreType == 's')
    {
        fileName << "_map";
    }
	else
    {
        static_assert(ScoreType != ScoreType, "Unknown score type.");
    }

    fileName << i << extension;
	return fileName.str();
}

int main(int argc, char* argv[])
{
    constexpr auto cMutations = 18;
    constexpr auto cSamples = 58;
    constexpr auto cTreeType = 'm'; // 'm': mutation tree, 't': rooted binary leaf-labelled tree
    constexpr auto cScoreType = 'm'; // 'm': default, 's': Sample attachment points are marginalized out
    constexpr auto cSeed = 1UL;
    constexpr auto cTrueTree = false;

	/****************   begin timing  *********************/
			clock_t begin=clock();
	/****************************************************/

	/**  Get command line parameters  **/
    auto inputFile = parameterOption<std::string>(argv, argv + argc, "-i");
    auto repetitions = parameterOption<std::size_t>(argv, argv + argc, "-r");
    auto chainLength = parameterOption<std::size_t>(argv, argv + argc, "-l");
    auto fd = parameterOption<double>(argv, argv + argc, "-fd");
    auto ad = parameterOption<double, 2>(argv, argv + argc, "-ad");

    auto cc = [&argv, &argc](){
        if(parameterExists(argv, argv + argc, "-cc"))
        {
            return parameterOption<double>(argv, argv + argc, "-cc");
        }
        else
        {
            return 0.0;
        }
    }();

    auto errorRates = std::array{fd, ad[0], ad[1], cc};

    auto sampleStep = [&argv, &argc](){
        if(parameterExists(argv, argv + argc, "-p"))
        {
            return parameterOption<std::size_t>(argv, argv + argc, "-p");
        }
        else
        {
            return std::size_t{0};
        }
    }();

    auto errorRateMove = [&argv, &argc](){
        if(parameterExists(argv, argv + argc, "-e"))
        {
            return parameterOption<double>(argv, argv + argc, "-e");
        }
        else
        {
            return 0.0;
        }
    }();

    auto chi = [&argv, &argc](){
        if(parameterExists(argv, argv+argc, "-x"))
        {
            return parameterOption<double>(argv, argv+argc, "-x");
        }
        else
        {
            return 10.0;
        }
    }();

    auto priorStd = [&argv, &argc](){
        if(parameterExists(argv, argv+argc, "-sd"))
        {
            return parameterOption<double>(argv, argv+argc, "-sd");
        }
        else
        {
            return 0.1;
        }
    }();

    auto outputName = [&argv, &argc, &inputFile](){
        if(parameterExists(argv, argv+argc, "-o"))
        {
            return parameterOption<std::string>(argv, argv+argc, "-o");
        }
        else
        {
            return inputFile.substr(0, inputFile.find_last_of("."));
        }
    }();

    auto geneNames = [&argv, &argc](){
        std::array<std::string, cMutations + 1> names;
        if(parameterExists(argv, argv+argc, "-names"))
        {
            auto fileName = parameterOption<std::string>(argv, argv+argc, "-names");
            std::ifstream input(fileName);
            if(!input)
            {
                std::cerr << "Cannot open gene names file " << fileName << std::endl;
                std::abort();
            }

            for(std::size_t i = 0; i < cMutations; ++i)
            {
                input >> names[i];
            }
            names[cMutations] = "Root";
        }
        else
        {
            for(std::size_t i = 0; i <= cMutations; ++i)
            {
                std::stringstream id;
                id << i+1;
                names[i] = id.str();
            }
        }
        return names;
    }();

    auto attachSamples = parameterExists(argv, argv + argc, "-a");

    auto maxTreeListSize = [&argv, &argc](){
        if(parameterExists(argv, argv+argc, "-max_treelist_size"))
        {
            return parameterOption<std::size_t>(argv, argv + argc, "-max_treelist_size");
        }
        else
        {
            return std::numeric_limits<std::size_t>::max();
        }
    }();

    auto gamma = [&argv, &argc](){
        if(parameterExists(argv, argv+argc, "-g"))
        {
            return parameterOption<double>(argv, argv + argc, "-g");
        }
        else
        {
            return 1.0;
        }
    }();

    auto useTreeList = !parameterExists(argv, argv + argc, "-no_tree_list");

    auto treeMoves = [&argv, &argc](){
        if(parameterExists(argv, argv+argc, "-move_probs"))
        {
            auto probs = parameterOption<double, cTreeType == 'm' ? 3 : 2>(argv, argv+argc, "-move_probs");
            if(auto sum = std::accumulate(probs.begin(), probs.end(), 0.0); sum != 1.0)
            {
                for(auto &v : probs)
                {
                    v /= sum;
                }
            }
            return probs;
        }
        else
        {
            if constexpr(cTreeType=='m')
            {
                return std::array{0.55, 0.4, 0.05};
            }
            else
            {
                return std::array{0.4, 0.6};
            }
        }
    }();

    // moveProbs: cTreeType == 'm': change beta / prune&re-attach / swap node labels / swap subtrees
    // moveProbs: cTreeType == 't': change beta / prune&re-attach / swap leaf labels
	std::array<double, treeMoves.size() + 1> moveProbs;
    moveProbs[0] = errorRateMove;
    for(std::size_t i = 1; i<moveProbs.size(); ++i)
    {
        moveProbs[i] = treeMoves[i-1];
    }

    auto dataMatrix = getDataMatrix<cMutations, cSamples>(inputFile);



	std::vector<struct treeBeta> optimalTrees;            // list of optimal tree/beta combinations found by MCMC
	std::string sampleOutput;                            // the samples taken in the MCMC as a string for outputting

	/* initialize the random number generator, either with a user defined seed, or a random number */
//    RandomGenerator().seed(1);

	/** get the true parent vector from GraphViz file if available (for simulated data only)  **/
	int* trueParentVec = NULL;
//	if(trueTreeComp==true){ trueParentVec = getTrueTree<cMutations>(trueTreeFileName); }

	/**  Find best scoring trees by MCMC  **/
	sampleOutput = runMCMCbeta(optimalTrees, errorRates.data(), repetitions, chainLength, gamma, std::vector(moveProbs.begin(), moveProbs.end()), cMutations, cSamples, toPointer(dataMatrix), cScoreType, trueParentVec, sampleStep, sampleStep != 0, chi, priorStd, useTreeList, cTreeType);


	/***  output results  ***/

	string prefix = outputName;

	/* output the samples taken in the MCMC */
	stringstream sampleOutputFile;
	sampleOutputFile << prefix << ".samples";
	writeToFile(sampleOutput, sampleOutputFile.str());
	cout << "samples from posterior written to: " << sampleOutputFile.str() << "\n";

	/* output the optimal trees found in individual files */

	double** logScores = getLogScores(fd, ad[0], ad[1], cc);
	int parentVectorSize = cMutations;
	if(cTreeType=='t'){parentVectorSize = (2*cSamples)-2;}                     // transposed case: binary tree, m leafs and m-1 inner nodes, root has no parent
	int outputSize = std::min(optimalTrees.size(), maxTreeListSize);
	for(int i=0; i<outputSize; i++){

		int* parentVector = optimalTrees.at(i).tree;
		bool** ancMatrix = parentVector2ancMatrix(parentVector, parentVectorSize);
		vector<vector<int> > childLists = getChildListFromParentVector(parentVector, parentVectorSize);

		stringstream newick;
		string outputFile = getOutputFileName<cScoreType>(i, prefix, ".newick");
		newick << getNewickCode(childLists, parentVectorSize) << "\n";
		writeToFile(newick.str(), outputFile);
		outputFile = getOutputFileName<cScoreType>(i, prefix, ".gv");	                                // print out tree as newick code

		if(errorRateMove != 0.0){
			updateLogScores(logScores, optimalTrees[i].beta);
		}

		if(cTreeType == 'm'){
			string output;
			output = getGraphVizFileContentNames(parentVector, parentVectorSize, std::vector(geneNames.begin(), geneNames.end()), attachSamples, ancMatrix, cSamples, logScores, toPointer(dataMatrix));
			writeToFile(output, outputFile);
		}
		else{
			int* bestPlacement = getHighestOptPlacementVector(toPointer(dataMatrix), cMutations, cSamples, logScores, ancMatrix);
			vector<string> bestBinTreeLabels = getBinTreeNodeLabels((2*cSamples)-1, bestPlacement, cMutations, std::vector(geneNames.begin(), geneNames.end()));
			//cout << getGraphVizBinTree(optimalTrees.at(0).tree, (2*m)-1, m, bestBinTreeLabels);
		}

		free_boolMatrix(ancMatrix);

	}

	delete [] logScores[0];
	delete [] logScores;
	cout << optimalTrees.size() << " opt trees \n";
	emptyVectorFast(optimalTrees, cMutations);


	/****************   end timing  *********************/
  		clock_t end=clock();
  		double diffticks=end-begin;
  		double diffms=(diffticks*1000)/CLOCKS_PER_SEC;
  		cout << "Time elapsed: " << diffms << " ms"<< endl;
  	/****************************************************/
}
