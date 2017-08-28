// -------------------------------------------------------------------------
// cd_spnet.h - Header file for Genetic Algorithm of sparse networks
//
// written by Aditya Karnam, August 2017
//
// ------------------------------------------------------------------------

#ifndef _CDSPNET_H
#define _CDSPNET_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "timer.h"
#include "network.h"

const int GA_QUIET = 1;  // set to one to eliminate output to screen
const int GA_VERBOSE = 0;  // set to one to display maximum output to screen
const int GA_DEBUG = 1; //set to one to display debugging for function

class Chromosome {
	friend class Population; // Population class allowed access to private functions
	friend class GA; // GA class allowed access to private functions
	private:
		int *edgeID; // variable number of edgeIDs for each chromosome
};

class GA {
	public:
		GA(Network&, int = 0, int = 0, int = 0, int = 0); // create GA with populationSize for each Individual
		~GA(); // destructor to manage garbage collection
		void generate_GA(); // test func
		int removeEdgeByID(int); // remove edge by EdgeID 
		void getOriginalEdgeIDS(); // returns the EdgeIDS for the network
		int removeEdgeByPosition(int, int); // removes edge by (v1,v2) position
		int getEdgeIDIndex(int, int); // get index of EdgeID from originalEdgeIDS
	private:
		Network *gaSparseNetwork; // generated network reference
		Chromosome *chromosomes; // array of individuals/chromosomes
		int populationSize; // number of Individuals in the populations
		int networkNumVertices; // store the number of vertices in the graph
		int networkNumEdges; // store the number of edges in the graph
		int numGenerations; // number of generations of GA
		int *originalEdgeIDS; // network EdgeIDs for lookup in Genetic Algo
		int **simpleGAChromosome; // basic chromosome to generate initial populations
		int addEdgeByEdgeID(int); // add edge by EdgeID to the network
		int generateRandomNumber(int, int); // generate random number between min and max range
};

#endif
