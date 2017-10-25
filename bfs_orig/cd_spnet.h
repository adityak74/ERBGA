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
const int GA_DEBUG = 1; // set to one to display debugging for function
const int GA_DEBUG_FILE = 0; // prints debug to file
const int GA_TOURNAMENT_SIZE = 2; // tournament size for the seleciton operator
const int GA_NUM_COMMUNITY = 3; // original community size to start with
const double GA_CROSSOVER_RATE = 0.8; // max chr size used for crossover
const double GA_REPRODUCTION_RATE = 0.1; // rate of reproduction producing the offsprings
const double GA_MUTATION_RATE = 0.2; // rate of mutation producing the offspring

// macro for calculating the size from actual array to bit array size
#define ARRAY_SIZE(x) (x/8+(!!(x%8)))

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
		int length; // non negative chromsome length
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
		double averageFitnessForPopulation(); // calculates the average fitness of population
		char get_bit(char *array, int index); // get the bit value for chromosomeBitArr
		void toggle_bit(char *array, int index); // toggle the bit value for chromosomeBitArr
		int getBitAt(int chrIndex, int pos, int popState); // gets the bit value at position
		void initializeRates(); // initalize GA params , normalize if needed
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
		char ***chromosomesBitArr; // bit array for edge state in the chromosome
		double crossover_rate; // crossover produced offsprings
		double mutation_rate; // mutation produced offsprings
		double reproduction_rate; // reproduction produced offsprings
};

#endif
