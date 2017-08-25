/****************************************************************************
*
*	cd_spnet.cpp:	Code for finding communities using Genetic Programming. 
						Outputs optimal Qs value.
*
*                       Aditya Karnam
*                       August 2017
*
****************************************************************************/
  

#include "network.h"
#include "cd_spnet.h"

// need to filter which EdgeID exist in original network
// how to access network here?
int GA::removeEdgeByID(int edgeID) {
	gaSparseNetwork->removeEdge(edgeID / networkNumVertices, edgeID % networkNumVertices);
	return 0;
}

GA::GA(Network *sparseNetwork, int populationSize, int numNodes, int numEdges) {
	if ((individuals = new Individual[populationSize]) == NULL) // allocate memory to set of individuals
    	fatal("memory not allocated");
    networkNumVertices = numNodes;
    networkNumEdges = numEdges;
    gaSparseNetwork = sparseNetwork;
}

GA::~GA() {
	delete [] individuals;
}

// generate random number in range (min, max)
int GA::generateRandomNumber(int min, int max) {
	return ((rand() % max) + min); 	
}

// max EdgeID can be (networkNumVertices)^2 for generating random edgeIDs
//num of edges removed can be a max upto (0, networkNumEdges/2)
void GA::generate_GA() {
	std::cout << "Before removing edge : " << gaSparseNetwork->getNumEdges() << std::endl;
	gaSparseNetwork->getOriginalEdgeIDS();
	removeEdgeByID(23);
	std::cout << "After removing edge : " << gaSparseNetwork->getNumEdges() << std::endl;
}