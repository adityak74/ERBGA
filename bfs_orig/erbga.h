// -------------------------------------------------------------------------
// erbga.h - Header file for Genetic Algorithm of sparse networks
//
// written by Aditya Karnam, August 2017
//
// ------------------------------------------------------------------------

#ifndef _ERBGA_H
#define _ERBGA_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <algorithm>
#include <iostream>
#include <unistd.h>
#include <fstream>
#include <cstdarg>
#include <time.h>
#include "timer.h"
#include "network.h"

const int GA_QUIET = 0;  // set to one to eliminate output to screen
const int GA_VERBOSE = 0;  // set to one to display maximum output to screen
const int GA_DEBUG = 1; // set to one to display debugging for function
const int GA_DEBUG_FILE = 0; // prints debug to file
const int GA_DEBUG_L2 = 0; // level 2 debugging

// program defaults
const int GA_POPULATION_SIZE = 100; // default population size
const int GA_NUM_GENERATIONS = 250; // default number of iterations for GA


const double GA_RANDOM_POP_PERCENT = 0.75; // set a minimum of 75 percent here
const int GA_TOURNAMENT_SIZE = 7; // tournament size for the seleciton operator
const double GA_CROSSOVER_SIZE = 0.5; // percentage of chromosome used for crossover
const double GA_CROSSOVER_RATE = 0.5; // max chr size used for crossover
const double GA_ELITISM_RATE = 0.1; // rate of reproduction producing the offsprings
const double GA_MUTATION_RATE = 0.5; // rate of mutation producing the offspring (1 / PopSize) 
const double GA_GENE_REPAIR_PERCENT = 0.5; // nodes of higher degree are repaired
const double GA_GENE_REPAIR_CHANCE = 0.2; // percent chance of edge added back

const int GA_FITNESS_MODULARITY = 1; // use modularity as fitness, else use Qs

const int GA_CROSSOVER_TYPE = 2; // ( 0 - One point, 1 - One way, 2 - Uniform)

const int GA_CROSSOVER_PARENT_CHILD_BEST = 0; // keep the best in parent and child

const int GA_ONE_WAY_CROSSOVER = 0; // use one way crossover

const int GA_ONE_POINT_CROSSOVER = 0; // one point crossover

const double GA_TOL = 0.0000000001;

const std::string GA_LOG_FILE = "ga_run-"; // log filename
const std::string GA_POP_FILE = "ga_pop-"; // random filename
const std::string GA_BST_FILE = "ga_bst-"; // best solutions log
const std::string GA_BST_AVG_FITNESS_RUN = "ga_bst_avg_run-"; // reports best and mean fitness
const std::string GA_BST_GML = "ga_bst_final-"; // name for final GML
const std::string GA_POPULATION_DATA = "ga_pop_data-"; // name for final GML

// macro for calculating the size from actual array to bit array size
#define ARRAY_SIZE(x) (x/8+(!!(x%8)))

// chromosome mapping dependency structure for std::sort
typedef struct chromosome_map
{
	double fitness;
	int chromosome_index;
};

// map degree and node index
typedef struct node_degree_map
{
	int node_index;
	int degree;
};

class Chromosome
{
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
		int removeEdgeByID(long long int); // remove edge by EdgeID 
		void getOriginalEdgeIDS(); // returns the EdgeIDS for the network
		int removeEdgeByPosition(int, int); // removes edge by (v1,v2) position
		int getEdgeIDIndex(long long int); // get index of EdgeID from originalEdgeIDS
		void getFitness(); // Calculate the fitness of the chromosome "i"
		double calculateFitness(int = -1, int = 0); // calculates the fitness of the chromosome
		double averageFitnessForPopulation(); // calculates the average fitness of population
		void toggle_bit(char *array, int index); // toggle the bit value for chromosomeBitArr
		int get_bit(int chrIndex, int pos, int popState); // gets the bit value at position
		void toggle_bit(int chrIndex, int pos, int popState); // toggles the bit value
		void initializeRates(); // initalize GA params , normalize if needed
		void set_bit(int chrIndex, int pos, int popState); // set bit
		void unset_bit(int chrIndex, int pos, int popState); // unsets bit
		void mutate(int, int); // does mutation on the chromosome index
		void chromosome_g2p(int, int); // maps the bitArr to chromosomes
		void printPopData(int = 0); // prints the current population
		void printChromosome(int, int); // prints the chromosome to log file
		void move_chromosome_to_next_gen(int, int); // moves chromosome to next generation
		void set_data_name(char*); //set dataset name
		void printChromosomes(int); // print the whole chromosome at depth
		void repairGene(int, int); // repair the gene for high degree nodes
		void exportPopulationData(int); // write the bits to a backup file 

	  private:
		Network *gaSparseNetwork; // generated network reference
		Chromosome *chromosomes; // array of individuals/chromosomes
		int populationSize; // number of Individuals in the populations
		int numGenerations; // number of generations of GA
		int networkNumVertices; // store the number of vertices in the graph
		int networkNumEdges; // store the number of edges in the graph
		long long int *originalEdgeIDS; // network EdgeIDs for lookup in Genetic Algo
		int **simpleGAChromosome; // basic chromosome to generate initial populations
		int addEdgeByEdgeID(long long int); // add edge by EdgeID to the network
		int generateRandomNumber(int, int); // generate random number between min and max range
		char ***chromosomesBitArr; // bit array for edge state in the chromosome
		double crossover_rate; // crossover produced offsprings
		double mutation_rate; // mutation produced offsprings
		double reproduction_rate; // reproduction produced offsprings
		int numGenomeMutations; // number of genomes mutated in the chromosome
		int numGenesRepairSize; // number of nodes to be repaired after the breeding phase
		int next = 1;
		int prev = 0;
		char *dataset_name;
		node_degree_map *ndmap; // map node and degree to sort with context
};

#endif
