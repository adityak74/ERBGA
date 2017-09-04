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

const int GA_QUIET = 0;  // set to one to eliminate output to screen
const int GA_VERBOSE = 0;  // set to one to display maximum output to screen
const int GA_DEBUG = 0; // set to one to display debugging for function
const int GA_DEBUG_FILE = 1; // prints debug to file
const int GA_TOURNAMENT_SIZE = 3; // tournament size for the seleciton operator
const int GA_NUM_COMMUNITY = 3; // original community size to start with

class Chromosome {
	friend class Population; // Population class allowed access to private functions
	friend class GA; // GA class allowed access to private functions
	public:
		double getFitness(); // get the fitness of the chromosome
		double calculateFitness(); // calculates the fitness of the chromosome
 	private:
		int *edgeIDS; // variable number of edgeIDs for each chromosome
		double fitness; // store the fitness of the chromosome
		Network *gaSparseNetworkChr; // generated network reference
		int networkNumVerticesChr; // store the number of vertices in the graph
		int networkNumEdgesChr; // store the number of edges in the graph
};

class GA {
	friend class Chromosome;
	public:
		GA(Network&, int = 0, int = 0, int = 0, int = 0); // create GA with populationSize for each Individual
		~GA(); // destructor to manage garbage collection
		void generate_GA(); // test func
		int removeEdgeByID(int); // remove edge by EdgeID 
		void getOriginalEdgeIDS(); // returns the EdgeIDS for the network
		int removeEdgeByPosition(int, int); // removes edge by (v1,v2) position
		int getEdgeIDIndex(int); // get index of EdgeID from originalEdgeIDS
		void getFitness(); // Calculate the fitness of the chromosome "i"
		void calculateFitness(int = -1); // calculates the fitness of the chromosome
		double averageFitnessForPopulation();
	private:
		Network *gaSparseNetwork; // generated network reference
		Chromosome *chromosomes; // array of individuals/chromosomes
		int populationSize; // number of Individuals in the populations
		int numGenerations; // number of generations of GA
		int networkNumVertices; // store the number of vertices in the graph
		int networkNumEdges; // store the number of edges in the graph
		int *originalEdgeIDS; // network EdgeIDs for lookup in Genetic Algo
		int **simpleGAChromosome; // basic chromosome to generate initial populations
		int addEdgeByEdgeID(int); // add edge by EdgeID to the network
		int generateRandomNumber(int, int); // generate random number between min and max range
};

#endif
