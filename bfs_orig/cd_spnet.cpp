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

// binary search part
// A recursive binary search function. It returns location of x in
// given array arr[l..r] is present, otherwise -1
int binarySearch(int arr[], int l, int r, int x)
{
   if (r >= l)
   {
        int mid = l + (r - l)/2;
 
        // If the element is present at the middle itself
        if (arr[mid] == x)  return mid;
 
        // If element is smaller than mid, then it can only be present
        // in left subarray
        if (arr[mid] > x) return binarySearch(arr, l, mid-1, x);
 
        // Else the element can only be present in right subarray
        return binarySearch(arr, mid+1, r, x);
   }
 
   // We reach here when element is not present in array
   return -1;
}
// ends

// need to filter which EdgeID exist in original network
// how to access network here?
int GA::removeEdgeByID(int edgeID) {
	gaSparseNetwork->removeEdge(edgeID / networkNumVertices, edgeID % networkNumVertices);
	return 1;
}

int GA::removeEdgeByPosition(int v1, int v2) {
	gaSparseNetwork->removeEdge(v1, v2);
	return 1;
}

// get index by using binary search
int GA::getEdgeIDIndex(int v1, int v2) {
	return binarySearch(originalEdgeIDS, 0, networkNumVertices-1, v1*networkNumVertices+v2);
}

GA::GA(Network &sparseNetwork, int popSize, int generations, int numNodes, int numEdges) {
	if ((chromosomes = new Chromosome[popSize]) == NULL) // allocate memory to set of Chromosomes
    	fatal("memory not allocated");
    populationSize = popSize;
    networkNumVertices = numNodes;
    networkNumEdges = numEdges;
    gaSparseNetwork = &sparseNetwork;
    numGenerations = generations;

    if ((simpleGAChromosome = new int*[popSize]) == NULL) // allocate memory to set of Basic GA Chromosomes
    	fatal("memory not allocated");

    // initalize size of basic GA chromosome to be equal to numVertices
    for (int i = 0; i < popSize; ++i) {
    	simpleGAChromosome[i] = new int[networkNumVertices]; // size of array equal to numVertices
    }

    if ((originalEdgeIDS = new int[numEdges]) == NULL)
    	fatal("memory not allocated");

    int edgePos = 0;

    // generate originalEdgeIDS for future use
    for (int i = 0; i < networkNumVertices; ++i)
    {
    	Edge *edgePtr;
    	edgePtr = &gaSparseNetwork->vertices[i].firstEdge;

    	if (edgePtr->next == 0) // no edges for this vertex
    		continue;

    	while(edgePtr->next != 0) {
    		edgePtr = edgePtr->next;
    		
    		// Avoid RETAINSYMMETRIC
    		if(i < edgePtr->target){
    			if(GA_DEBUG)
    				std::cout << "#"<< edgePos << " :-: " << networkNumVertices * (i) << " + " << (edgePtr->target) << std::endl;
    			originalEdgeIDS[edgePos++] = networkNumVertices*(i)+(edgePtr->target);
    		}
    	}
    }

    // sort originalEdgeIDS for binary search use
    std::sort(originalEdgeIDS, originalEdgeIDS + networkNumEdges);

    // generate initial population
    for (int i = 0; i < populationSize; ++i) {
    	for (int j = 0; j < networkNumVertices; ++j) {
    		simpleGAChromosome[i][j] = generateRandomNumber(0, networkNumVertices/4);
    	}
    }

    // initial population using Basic Chromosome (Community ID based)
    if(GA_DEBUG) {
    	for (int i = 0; i < networkNumVertices; ++i) {
	    	for (int j = 0; j < populationSize; ++j) {
	    		std::cout << simpleGAChromosome[j][i] << "\t";
	    	}
	    	std::cout << std::endl;
	    }
    }

    // remove edges from the random population generated by assigning community IDS
    for (int i = 0; i < populationSize; ++i) {

    	// edge ID state maintenance
    	int edgeIDState[networkNumEdges];
    	
    	// before removing all edges are in
    	for (int k = 0; k < networkNumEdges; ++k)
    		edgeIDState[k] = 1;

    	for (int j = 0; j < networkNumVertices; ++j) {
	    	Edge *edgePtr;
	    	edgePtr = &gaSparseNetwork->vertices[j].firstEdge;

	    	

	    	if (edgePtr->next == 0) // no edges for this vertex
	    		continue;

	    	while(edgePtr->next != 0) {
	    		edgePtr = edgePtr->next;
	    		// Avoid RETAINSYMMETRIC
	    		if(j < edgePtr->target){

	    			if(simpleGAChromosome[i][j] != simpleGAChromosome[i][edgePtr->target]) {
	    				std::cout << j << " : " << edgePtr->target << std::endl;
	    				if(gaSparseNetwork->haveEdge(j, edgePtr->target))
	    					removeEdgeByID(networkNumVertices*(j)+(edgePtr->target));
	    				edgeIDState[getEdgeIDIndex(j, edgePtr->target)] = -1;
	    			}
	    		}
	    	}
	    }

	    // after one chromosome generation
	    for (int k = 0; k < networkNumEdges; ++k){
    		if(edgeIDState[k] == -1){
    			std::cout << originalEdgeIDS[k] << "\t";
    		}
	    }
	    std::cout << std::endl;

    }

}

GA::~GA() {
	delete [] chromosomes;
	delete [] originalEdgeIDS;
	delete [] simpleGAChromosome;
}

void GA::getOriginalEdgeIDS() {
  	for (int i = 0; i < networkNumEdges; ++i) {
      std::cout << originalEdgeIDS[i] << "\t";
    }
    std::cout <<  std::endl;
}

// generate random number in range (min, max)
int GA::generateRandomNumber(int min, int max) {
	return ((rand() % max) + min); 	
}

// max EdgeID can be (networkNumVertices)^2 for generating random edgeIDs
//num of edges removed can be a max upto (2, networkNumEdges/2)
void GA::generate_GA() {
	
}