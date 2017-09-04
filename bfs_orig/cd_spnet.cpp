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

double Chromosome::getFitness() { 
	return fitness;
}

double Chromosome::calculateFitness() {

	for (int i = 0; i < networkNumEdgesChr; ++i)
		if(edgeIDS[i] != -1){
			gaSparseNetworkChr->removeEdge(edgeIDS[i] / networkNumVerticesChr, edgeIDS[i] % networkNumVerticesChr);
		}

	double calc_fitness = gaSparseNetworkChr->q_calc("ga.out");

	for (int i = 0; i < networkNumEdgesChr; ++i)
		if(edgeIDS[i] != -1){
			gaSparseNetworkChr->addEdge(edgeIDS[i] / networkNumVerticesChr, edgeIDS[i] % networkNumVerticesChr, 1.0);
		}

	fitness = calc_fitness;

	return fitness;
}

// if chr_index = -1 then all else calcu fitness of chromosome[chr_index]
void GA::calculateFitness(int chr_index) {
	if(chr_index == -1)
		for (int i = 0; i < populationSize; ++i)
			chromosomes[i].calculateFitness();
	else
		chromosomes[chr_index].calculateFitness();
}

// need to filter which EdgeID exist in original network
// how to access network here?
int GA::removeEdgeByID(int edgeID) {
	return gaSparseNetwork->removeEdge(edgeID / networkNumVertices, edgeID % networkNumVertices);
}

int GA::removeEdgeByPosition(int v1, int v2) {
	return gaSparseNetwork->removeEdge(v1, v2);
}

int GA::addEdgeByEdgeID(int edgeID) {
	return gaSparseNetwork->addEdge(edgeID / networkNumVertices, edgeID % networkNumVertices, 1.0);
}

// get index by using binary search
int GA::getEdgeIDIndex(int edgeID) {
	return binarySearch(originalEdgeIDS, 0, networkNumEdges, edgeID);
}

GA::GA(Network &sparseNetwork, int popSize, int generations, int numNodes, int numEdges) {

	srand (time(NULL));
	// set class data memebers
	populationSize = popSize;
    networkNumVertices = numNodes;
    networkNumEdges = numEdges;
    gaSparseNetwork = &sparseNetwork;
    numGenerations = generations;

	if ((chromosomes = new Chromosome[popSize]) == NULL) // allocate memory to set of Chromosomes
    	fatal("memory not allocated");

    if ((simpleGAChromosome = new int*[popSize]) == NULL) // allocate memory to set of Basic GA Chromosomes
    	fatal("memory not allocated");

    // initalize size of basic GA chromosome to be equal to numVertices
    for (int i = 0; i < popSize; ++i) {
    	simpleGAChromosome[i] = new int[networkNumVertices]; // size of array equal to numVertices
    	chromosomes[i].edgeIDS = new int[numEdges]; // size of each chromosome maximum to be all the edges
    	chromosomes[i].gaSparseNetworkChr = &sparseNetwork; // store the reference for use in Chromosome class
    	chromosomes[i].networkNumVerticesChr = networkNumVertices;
    	chromosomes[i].networkNumEdgesChr = networkNumEdges;
    }

    if ((originalEdgeIDS = new int[numEdges]) == NULL)
    	fatal("memory not allocated");

    int edgePos = 0;

    // generate originalEdgeIDS for future use
    for (int i = 0; i < networkNumVertices; ++i)
    {
    	Edge *edgePtr;
    	edgePtr = &gaSparseNetwork->vertices[i].firstEdge;

    	if (edgePtr->next == NULL) // no edges for this vertex
    		continue;

    	
    	while(edgePtr->next != NULL) {
    		edgePtr = edgePtr->next;

    		// Avoid RETAINSYMMETRIC
    		if(i < edgePtr->target){
    			if(GA_DEBUG)
    				std::cout << "#"<< edgePos << " :-: " << networkNumVertices * (i) << " + " << (edgePtr->target) << std::endl;
    			originalEdgeIDS[edgePos++] = networkNumVertices*(i)+(edgePtr->target);
    		}
    	}
    	// std::cout << "NULL" << std::endl;
    }

    // sort originalEdgeIDS for binary search use
    std::sort(originalEdgeIDS, originalEdgeIDS + networkNumEdges);

    // original Edge ID Array Print
    if(GA_DEBUG){
	    std::cout << "original Edge ID Array : \n";
	    for (int k = 0; k < networkNumEdges; ++k){
		    		std::cout << originalEdgeIDS[k] << "\t";
			    }
		std::cout << "\n-=-=-=-=-=-=-=-=-=-=-=-=-=-" << std::endl;
	}

    // initialize each chromosome with -1 (empty)
    for (int i = 0; i < populationSize; ++i) {
    	for (int j = 0; j < networkNumEdges; ++j) {
    		chromosomes[i].edgeIDS[j] = -1;
    	}
    }

    // generate initial population
    for (int i = 0; i < populationSize; ++i) {
    	for (int j = 0; j < networkNumVertices; ++j) {
    		simpleGAChromosome[i][j] = generateRandomNumber(0, GA_NUM_COMMUNITY);
    		// originalNUmClusters - command line/ header
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

    // edge ID state maintenance
    int edgeIDState[networkNumEdges];
    	
    // before removing all edges are in
    for (int k = 0; k < networkNumEdges; ++k)
    	edgeIDState[k] = 1;

    // remove edges from the random population generated by assigning community IDS
    for (int i = 0; i < populationSize; ++i) {

    	if(GA_DEBUG) {
    		//before each chromosome generated
    		std::cout << "-=-=-=-=-=-=-=-=-=-=-=-=-=-" << std::endl;
	    	std::cout << " Edge State ID Array " << std::endl;
	    	std::cout << " Total Num Edges : " << gaSparseNetwork->getNumEdges() << std::endl;
	    	for (int k = 0; k < networkNumEdges; ++k){
	    		std::cout << edgeIDState[k] << "\t";
		    }
	    	std::cout << "\n-=-=-=-=-=-=-=-=-=-=-=-=-=-" << std::endl;
    	}

    	for (int j = 0; j < networkNumVertices; ++j) {
	    	Edge *edgePtr;
	    	edgePtr = &gaSparseNetwork->vertices[j].firstEdge;

	    	if (edgePtr->next == NULL) // no edges for this vertex
	    		continue;

	    	while(edgePtr->next != NULL) {
	    		edgePtr = edgePtr->next;
	    		// Avoid RETAINSYMMETRIC
	    		if(j < edgePtr->target){
	    			if(simpleGAChromosome[i][j] != simpleGAChromosome[i][edgePtr->target]) {
	    				// std::cout << "(" << j << "," << edgePtr->target << ") EdgeID :-> " << originalEdgeIDS[getEdgeIDIndex(j, edgePtr->target)] << std::endl; 
	    				// overhead check can be removed as we are sure to be within bounds
	    				if(gaSparseNetwork->haveEdge(j, edgePtr->target)){
	    					int edgeIDToRemove = networkNumVertices*(j)+(edgePtr->target);
	    					removeEdgeByID(edgeIDToRemove);
	    					// fix here
	    					edgeIDState[getEdgeIDIndex(edgeIDToRemove)] = -1;
	    					if(GA_DEBUG)
	    						std::cout << "Edge state update :-> " << edgeIDToRemove << " :-> " << getEdgeIDIndex(edgeIDToRemove) << std::endl;
	    				}
	    			}
	    		}
	    	}
	    }

	    if(GA_DEBUG) {
    		//afetr each chromosome generated
    		std::cout << "-=-=-=-=-=-=-=-=-=-=-=-=-=-" << std::endl;
	    	std::cout << " Edge State ID Array after removing edges " << std::endl;
	    	std::cout << " Total Num Edges : " << gaSparseNetwork->getNumEdges() << std::endl;
	    	for (int k = 0; k < networkNumEdges; ++k){
	    		std::cout << edgeIDState[k] << "\t";
		    }
	    	std::cout << "\n-=-=-=-=-=-=-=-=-=-=-=-=-=-" << std::endl;
    	}

	    // add edges back after generating chromosome
	    // after one chromosome generation
	    int chromosome_edgeID_pos = 0;
	    if(GA_DEBUG)
	    	std::cout << "Removed edges for chromosome : "<< (i+1) << " : \n";
	    for (int k = 0; k < networkNumEdges; ++k){
    		if(edgeIDState[k] == -1){
    			edgeIDState[k] = addEdgeByEdgeID(originalEdgeIDS[k]);
    			chromosomes[i].edgeIDS[chromosome_edgeID_pos++] = originalEdgeIDS[k];
    		}
	    }
	    if(GA_DEBUG)
	    	std::cout << "\n";

	    

    }

    if(GA_DEBUG) {
    	std::cout << "Actual Chromosomes : " << std::endl;
    	for (int i = 0; i < populationSize; ++i)
    		std::cout << "CHR#" << (i+1) << "\t";
    	std::cout << std::endl;
    	for (int i = 0; i < networkNumEdges; ++i) {
	    	for (int j = 0; j < populationSize; ++j) {
	    		std::cout << chromosomes[j].edgeIDS[i] << "\t";
	    	}
	    	std::cout << std::endl;
	    }

	    if(GA_DEBUG && GA_DEBUG_FILE){
	    	char* outputFile = "chromosomes.ga";
		    FILE *output;
			if ((output = fopen(outputFile, "w")) == NULL)
				fatal ("Unable to open output file");

			fprintf(output,"Actual Chromosomes\n");  

			for (int i = 0; i < populationSize; ++i)
				fprintf(output,"CHR#%d\t\t", i+1);
			fprintf(output,"\n");
			for (int i = 0; i < networkNumEdges; ++i) {
		    	for (int j = 0; j < populationSize; ++j) {
		    		fprintf(output,"%d\t", chromosomes[j].edgeIDS[i]);
		    		fprintf(output,"\t\t");		
		    	}
		    	fprintf(output,"\n");
		    }

			if ((output = fopen(outputFile, "a")) == NULL)
			 	fatal ("Unable to open output file");
  
			fclose(output);
		}
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

// taken from http://preshing.com/20121224/how-to-generate-a-sequence-of-unique-random-integers/
unsigned int permuteQPR(unsigned int x)
{
    static const unsigned int prime = 4294967291;
    if (x >= prime)
        return x;  // The 5 integers out of range are mapped to themselves.
    unsigned int residue = ((unsigned long long) x * x) % prime;
    return (x <= prime / 2) ? residue : prime - residue;
}

// generate random number in range (min, max)
int GA::generateRandomNumber(int min, int max) {
	return ((permuteQPR(rand()) % max) + min); 	
}

double GA::averageFitnessForPopulation() {
	double total_fitness = 0.0;
	for (int i = 0; i < populationSize; ++i)
	{
		chromosomes[i].calculateFitness();
		total_fitness += chromosomes[i].getFitness();
	}
	return (total_fitness/populationSize);
}

// max EdgeID can be (networkNumVertices)^2 for generating random edgeIDs
//num of edges removed can be a max upto (2, networkNumEdges/2)
void GA::generate_GA() {
	// for (int i = 0; i < populationSize; ++i)
	// {
	// 	chromosomes[i].calculateFitness();
	// }
	// for (int i = 0; i < populationSize; ++i)
	// {
	// 	std::cout << "Fitness for CHR#" << (i+1) << " :-> " <<chromosomes[i].getFitness() << std::endl;
	// }
	// double avg_fitness = averageFitnessForPopulation();
	// std::cout << "Average Fitness : " << avg_fitness << std::endl;

	// tournament selection
	int nextGenChromosomeState[populationSize];
	for (int i = 0; i < populationSize; ++i) {
		nextGenChromosomeState[i] = 0; // no chromosome is selected 0, if selected 1
	}

	for (int i = 0; i < GA_TOURNAMENT_SIZE; ++i) {
		int index = generateRandomNumber(0, populationSize);
		while(nextGenChromosomeState[index]) {
			index = generateRandomNumber(0, populationSize);
		}
		nextGenChromosomeState[index] = 1;
	}

	for (int i = 0; i < populationSize; ++i) {
		std::cout << nextGenChromosomeState[i] << "\t";
		if ((i+1)%5==0 && i>0)
		{
			std::cout << std::endl;
		}
	}
	std::cout << std::endl;




}