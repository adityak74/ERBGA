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
#include "timer.h"
#include "network.h"

const int GA_QUIET = 1;  // set to one to eliminate output to screen
const int GA_VERBOSE = 0;  // set to one to display maximum output to screen

class Individual {
	friend class Population; // Population class allowed access to private functions
	friend class GA; // GA class allowed access to private functions
	private:
		int *edgeID; // variable number of edgeIDs for each individual
};

class GA {
	public:
		GA(Network*, int = 0, int = 0, int = 0); // create GA with populationSize for each Individual
		~GA(); // destructor to manage garbage collection
		int generateRandomNumber(int, int); // generate random number between min and max range
		void generate_GA(); // test func
		int removeEdgeByID(int); // remove edge by EdgeID 
	private:
		Network *gaSparseNetwork; // generated network reference
		Individual *individuals; // array of individuals
		int populationSize; // number of Individuals in the populations
		int networkNumVertices; // store the number of vertices in the graph
		int networkNumEdges; // store the number of edges in the graph
};

#endif
