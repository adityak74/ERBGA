/****************************************************************************
*
*	cd_spnet.cpp:	Code for finding communities using Genetic Programming. 
*						Outputs optimal Qs value.
*
*                       Aditya Karnam
*                       Aug, 2017 - Mar, 2018
*						v 0.1 - Added Genetic operators
*						v 0.2 - Added Elitism Operator
*						v 0.3 - One Point/Two Point Crossover
*						v 0.4 - Added Gene Repair Operator
*
*
****************************************************************************/

#include "network.h"
#include "cd_spnet.h"

time_t now;
char the_date[256];

// binary search part
// A recursive binary search function. It returns location of x in
// given array arr[l..r] is present, otherwise -1
int binarySearch(long long int arr[], int l, int r, int x)
{
	if (r >= l)
	{
		int mid = l + (r - l) / 2;

		// If the element is present at the middle itself
		if (arr[mid] == x)
			return mid;

		// If element is smaller than mid, then it can only be present
		// in left subarray
		if (arr[mid] > x)
			return binarySearch(arr, l, mid - 1, x);

		// Else the element can only be present in right subarray
		return binarySearch(arr, mid + 1, r, x);
	}

	// We reach here when element is not present in array
	return -1;
}
// ends

bool chromosome_sorter(const chromosome_map &a, const chromosome_map &b)
{
	return a.fitness > b.fitness;
}

bool node_sorter(const node_degree_map &a, const node_degree_map &b)
{
	return a.degree > b.degree;
}

double Chromosome::getFitness()
{
	return fitness;
}

double Chromosome::calculateFitness()
{
	double calc_fitness = 0.0f;
	if (GA_FITNESS_MODULARITY)
		calc_fitness = gaSparseNetworkChr->modularity("ga.out");
	else
		calc_fitness = gaSparseNetworkChr->q_calc("ga.out");
	fitness = calc_fitness;
	return calc_fitness;
}

// if chr_index = -1 then all else calcu fitness of chromosome[chr_index]
double GA::calculateFitness(int chr_index, int depth)
{
	double fitness = 0.0;
	if (chr_index < 0)
	{
		fatal("Invalid chromosomeindex");
	}
	else
	{
		chromosome_g2p(chr_index, depth); // convert the chromosome from genotype to phenotype
		fitness = chromosomes[chr_index].calculateFitness();
	}
	return fitness;
}

void GA::chromosome_g2p(int chromosomeIndex, int depth)
{

	if (GA_DEBUG_L2)
		std::cout << "\n --- Removing Edges from chromosome : " << chromosomeIndex + 1 << std::endl;

	for (int i = 0; i < networkNumEdges; ++i)
	{
		if (get_bit(chromosomeIndex, i, depth))
		{
			removeEdgeByID(originalEdgeIDS[i]);
			if (GA_DEBUG_L2)
				std::cout << "\n Removed Edge Index : " << i << " originalEdgeIDS : " << originalEdgeIDS[i] << std::endl;
		}
		else
			addEdgeByEdgeID(originalEdgeIDS[i]);
	}

	if (GA_DEBUG_L2)
		std::cout << "\n --- Removed edges from : " << chromosomeIndex + 1 << std::endl;
}

int GA::removeEdgeByID(long long int edgeID)
{
	return gaSparseNetwork->removeEdge(edgeID / networkNumVertices, edgeID % networkNumVertices);
}

int GA::removeEdgeByPosition(int v1, int v2)
{
	return gaSparseNetwork->removeEdge(v1, v2);
}

int GA::addEdgeByEdgeID(long long int edgeID)
{
	return gaSparseNetwork->addEdge(edgeID / networkNumVertices, edgeID % networkNumVertices, 1.0);
}

// get index by using binary search
int GA::getEdgeIDIndex(long long int edgeID)
{
	return binarySearch(originalEdgeIDS, 0, networkNumEdges, edgeID);
}

// toggle bit at index in array
void GA::toggle_bit(char *array, int index)
{
	array[index / 8] ^= 1 << (index % 8);
}

void GA::printChromosomes(int depth)
{
	fprintf(stderr, "-=-=-=-=- Complete Population bit array \n");
	for (int i = 0; i < populationSize; ++i)
	{
		for (int j = 0; j < networkNumEdges; ++j)
		{
			std::cout << get_bit(i, j, 0) << "\t";
		}
		std::cout << std::endl;
	}
}

GA::GA(Network &sparseNetwork, int popSize, int generations, int numNodes, int numEdges)
{

	the_date[0] = '\0';
	now = time(NULL);
	strftime(the_date, 256, "-%F-%T-", gmtime(&now));
	strcat(the_date, std::to_string(now).c_str());
	strcat(the_date, ".log");

	// set class data memebers
	populationSize = popSize;
	networkNumVertices = numNodes;
	networkNumEdges = numEdges;
	gaSparseNetwork = &sparseNetwork;
	numGenerations = generations;

	if ((ndmap = new node_degree_map[networkNumVertices]) == NULL) // allocate memory to node degree mapping
		fatal("memory not allocated");

	// map node degrees for use in gene repair
	for (int i = 0; i < networkNumVertices; i++)
	{
		ndmap[i].node_index = i;
		ndmap[i].degree = gaSparseNetwork->getDegree(i);
	}

	std::sort(ndmap, ndmap + networkNumVertices, &node_sorter);

	int max_degree = ndmap[0].degree; // store the max degree

	if (GA_DEBUG)
	{
		std::cout << "--- Nodes sorted by degree ascending : " << std::endl;
		for (int i = 0; i < networkNumVertices; i++)
		{
			if (ndmap[i].degree > -1)
				std::cout << ndmap[i].node_index << " , " << ndmap[i].degree << std::endl;
		}
	}

	initializeRates();

	if ((chromosomes = new Chromosome[popSize]) == NULL) // allocate memory to set of Chromosomes
		fatal("memory not allocated");

	if ((simpleGAChromosome = new int *[popSize]) == NULL) // allocate memory to set of Basic GA Chromosomes
		fatal("memory not allocated");

	// initalize size of basic GA chromosome to be equal to numVertices
	for (int i = 0; i < popSize; ++i)
	{
		simpleGAChromosome[i] = new int[networkNumVertices]; // size of array equal to numVertices
		chromosomes[i].edgeIDS = new int[numEdges];			 // size of each chromosome maximum to be all the edges
		chromosomes[i].gaSparseNetworkChr = &sparseNetwork;  // store the reference for use in Chromosome class
		chromosomes[i].networkNumVerticesChr = networkNumVertices;
		chromosomes[i].networkNumEdgesChr = networkNumEdges;
	}

	if ((originalEdgeIDS = new long long int[numEdges]) == NULL)
		fatal("memory not allocated");

	// allocate space for bit locations
	// show edge state (1=removed, 0=in_network)
	int bit_arr_pop_size = (popSize);
	int bit_arr_num_edges = ARRAY_SIZE(numEdges);

	if ((chromosomesBitArr = new char **[bit_arr_pop_size]) == NULL) // allocate memory to set of Basic GA Chromosomes
		fatal("memory not allocated");

	for (int i = 0; i < bit_arr_pop_size; ++i)
	{
		if ((chromosomesBitArr[i] = new char *[bit_arr_num_edges]) == NULL) // allocate memory to set of Basic GA Chromosomes
			fatal("memory not allocated");
	}

	for (int i = 0; i < bit_arr_pop_size; ++i)
	{
		for (int j = 0; j < bit_arr_num_edges; ++j)
		{
			// this dimesnion should be 2 times the number of subpopulation
			if ((chromosomesBitArr[i][j] = new char[2]) == NULL) // allocate memory to set of Basic GA Chromosomes
				fatal("memory not allocated");
		}
	}

	for (int i = 0; i < 2; ++i)
	{
		for (int j = 0; j < bit_arr_num_edges; ++j)
		{
			for (int k = 0; k < bit_arr_pop_size; ++k)
			{
				chromosomesBitArr[k][j][i] = 0;
			}
		}
	}

	int edgePos = 0;

	// generate originalEdgeIDS for future use
	for (long long int i = 0; i < networkNumVertices; ++i)
	{
		Edge *edgePtr;
		edgePtr = &gaSparseNetwork->vertices[i].firstEdge;

		if (edgePtr->next == NULL) // no edges for this vertex
			continue;

		while (edgePtr->next != NULL)
		{
			edgePtr = edgePtr->next;

			// Avoid RETAINSYMMETRIC
			if (i < edgePtr->target)
			{
				if (GA_DEBUG_L2) {
					std::cout << "#" << edgePos << " :-: " << networkNumVertices * (i) << " + " << (edgePtr->target) << " = " << networkNumVertices * (i) + (edgePtr->target) << std::endl;
					std::cout << "\t---(" << (i + 1) << " , " << (edgePtr->target+1) << " )\n\n";
				}
				originalEdgeIDS[edgePos++] = networkNumVertices * (i) + (edgePtr->target);
			}
		}
		// std::cout << "NULL" << std::endl;
	}

	// sort originalEdgeIDS for binary search use
	std::sort(originalEdgeIDS, originalEdgeIDS + networkNumEdges);

	// original Edge ID Array Print
	if (GA_DEBUG_L2)
	{
		std::cout << "original Edge ID Array : \n";
		for (int k = 0; k < networkNumEdges; ++k)
		{
			std::cout << k << "\t";
		}
		std::cout << "\n";
		for (int k = 0; k < networkNumEdges; ++k)
		{
			std::cout << originalEdgeIDS[k] << "\t";
		}
		std::cout << "\n-=-=-=-=-=-=-=-=-=-=-=-=-=-" << std::endl;
	}

	// initialize each chromosome with -1 (empty)
	for (int i = 0; i < populationSize; ++i)
	{
		for (int j = 0; j < networkNumEdges; ++j)
		{
			chromosomes[i].edgeIDS[j] = -1;
		}
	}

	// Main POP CREATION
	// generate initial population
	for (int i = 0; i < populationSize; ++i)
	{

		int hits = 0;
		for (int j = 0; j < networkNumEdges; ++j)
		{
			// changed to random community size
			int r = generateRandomNumber(0, networkNumEdges);
			if (r < round(GA_RANDOM_POP_PERCENT * networkNumEdges))
			{
				set_bit(i, j, prev);
				//hits++;
				//std ::cout << "------==r for chr #" << i + 1 << " :- " << r << "\n";
			}
			else
				unset_bit(i, j, prev);
			// originalNUmClusters - command line/ header
		}
		//std ::cout << "------==hits for chr #" << i + 1 << " : " << hits << " , " << (float)hits / networkNumEdges << "\n";
	}

	// Alternate POP CREATION
	// initialize each chromosome with -1 (empty)
	// std::cout << "\n-=-=-=-=-=-=generating chromosomes : \n";
	// int kcluster = 0;
	// for (int i = 0; i < populationSize; ++i)
	// {
	// 	kcluster = generateRandomNumber(9, 10);
	// 	// std::cout << "-->" << kcluster << "\t";
	// 	for (int j = 0; j < networkNumVertices; ++j)
	// 	{
	// 		int clusterLabel = generateRandomNumber(1, kcluster+1);
	// 		simpleGAChromosome[i][j] = clusterLabel;
	// 	}
	// }

	// if(GA_DEBUG) {
	// 	std::cout << "\n-=-=-=-=-=-=original chromosomes : \n";
	// 	for (int i = 0; i < populationSize; ++i)
	// 	{
	// 		for (int j = 0; j < networkNumVertices; ++j)
	// 		{
	// 			std::cout << simpleGAChromosome[i][j] << "\t";
	// 		}
	// 		std::cout << std::endl;
	// 	}
	// }

	// // edge ID state maintenance
	// int edgeIDState[networkNumEdges];

	// // before removing all edges are in
	// for (int k = 0; k < networkNumEdges; ++k)
	// 	edgeIDState[k] = 1;

	// // remove edges from the random population generated by assigning community IDS
	// for (int i = 0; i < populationSize; ++i) {

	// 	if(GA_DEBUG_L2) {
	// 		//before each chromosome generated
	// 		std::cout << "-=-=-=-=-=-=-=-=-=-=-=-=-=-" << std::endl;
	//     	std::cout << " Edge State ID Array " << std::endl;
	//     	std::cout << " Total Num Edges : " << gaSparseNetwork->getNumEdges() << std::endl;
	//     	for (int k = 0; k < networkNumEdges; ++k){
	//     		std::cout << edgeIDState[k] << "\t";
	// 	    }
	//     	std::cout << "\n-=-=-=-=-=-=-=-=-=-=-=-=-=-" << std::endl;
	// 	}

	// 	for (int j = 0; j < networkNumVertices; ++j) {
	//     	Edge *edgePtr;
	//     	edgePtr = &gaSparseNetwork->vertices[j].firstEdge;

	//     	if (edgePtr->next == NULL) // no edges for this vertex
	//     		continue;

	//     	while(edgePtr->next != NULL) {
	//     		edgePtr = edgePtr->next;
	//     		// Avoid RETAINSYMMETRIC
	//     		if(j < edgePtr->target){
	//     			if(simpleGAChromosome[i][j] != simpleGAChromosome[i][edgePtr->target]) {
	// 					// std::cout << "(" << j << "," << edgePtr->target << ") EdgeID :-> " << originalEdgeIDS[getEdgeIDIndex(networkNumVertices * (j) + (edgePtr->target))] << std::endl;
	// 					// overhead check can be removed as we are sure to be within bounds
	//     				if(gaSparseNetwork->haveEdge(j, edgePtr->target)){
	//     					int edgeIDToRemove = networkNumVertices*(j)+(edgePtr->target);
	//     					// removeEdgeByID(edgeIDToRemove);
	//     					// fix here
	//     					edgeIDState[getEdgeIDIndex(edgeIDToRemove)] = -1;
	//     					if(GA_DEBUG_L2)
	//     						std::cout << "Edge state update :-> " << edgeIDToRemove << " :-> " << getEdgeIDIndex(edgeIDToRemove) << std::endl;
	//     				}
	//     			}
	//     		}
	//     	}
	//     }

	//     if(GA_DEBUG_L2) {
	// 		//afetr each chromosome generated
	// 		std::cout << "-=-=-=-=-=-=-=-=-=-=-=-=-=-" << std::endl;
	//     	std::cout << " Edge State ID Array after removing edges " << std::endl;
	//     	std::cout << " Total Num Edges : " << gaSparseNetwork->getNumEdges() << std::endl;
	//     	for (int k = 0; k < networkNumEdges; ++k) {
	//     		std::cout << edgeIDState[k] << "\t";
	// 	    }
	//     	std::cout << "\n-=-=-=-=-=-=-=-=-=-=-=-=-=-" << std::endl;
	// 	}

	//     int chromosome_edgeID_pos = 0;
	//     if(GA_DEBUG_L2)
	//     	std::cout << "Removed edges for chromosome : "<< (i+1) << " : \n";

	//     // add edges back after generating chromosome
	//     // after one chromosome generation
	//     for (int k = 0; k < networkNumEdges; ++k){
	// 		if(edgeIDState[k] == -1){
	// 			edgeIDState[k] = 1;
	// 			// chromosomes[i].edgeIDS[chromosome_edgeID_pos++] = originalEdgeIDS[k];
	// 			set_bit(i, k, 0);
	// 		}
	//     }
	//     // assign length to the chromosome
	//     chromosomes[i].length = chromosome_edgeID_pos;

	//     if(GA_DEBUG_L2)
	//     	std::cout << "\n";

	// }

	// initial population using Basic Chromosome (Community ID based)
	if (GA_DEBUG_L2)
	{
		fprintf(stdout, "-=-=-=-=-Initial population bit array\n");
		for (int i = 0; i < populationSize; ++i)
		{
			for (int j = 0; j < networkNumVertices; ++j)
			{
				std::cout << get_bit(i, j, prev) << "\t";
			}
			std::cout << std::endl;
		}
	}

	if (GA_DEBUG_L2)
	{
		std::cout << "BitArr Dimensions : " << bit_arr_pop_size << " : " << bit_arr_num_edges << std::endl;
		std::cout << "SIZE OF BIT ARRAY : " << sizeof(chromosomesBitArr) << std::endl;
		std::cout << "--- BIT ARRAY ---\n\n";
		for (int i = 0; i < 2; ++i)
		{
			std::cout << "Depth (" << (i + 1) << ") :" << std::endl;
			for (int j = 0; j < bit_arr_pop_size; ++j)
			{
				for (int k = 0; k < bit_arr_num_edges; ++k)
				{
					int index = 0;
					while (index < 8)
					{
						std::cout << (k * 8) + index << " : (" << j << "," << k << "," << i << "," << index << ") -> ";
						std::cout << (1 & (chromosomesBitArr[j][k][i] >> (index))) << "\n";
						index++;
					}
				}
				std::cout << std::endl;
			}
			std::cout << "End of Depth\n\n";
		}
	}

	if (GA_DEBUG_L2)
	{
		std::cout << "Actual Chromosomes : " << std::endl;
		for (int i = 0; i < populationSize; ++i)
			std::cout << "CHR#" << (i + 1) << "\t";
		std::cout << std::endl;
		for (int i = 0; i < networkNumEdges; ++i)
		{
			for (int j = 0; j < populationSize; ++j)
			{
				std::cout << chromosomes[j].edgeIDS[i] << "\t";
			}
			std::cout << std::endl;
		}

		if (GA_DEBUG && GA_DEBUG_FILE)
		{
			char *outputFile = "chromosomes.ga";
			FILE *output;
			if ((output = fopen(outputFile, "w")) == NULL)
				fatal("Unable to open output file");

			fprintf(output, "Actual Chromosomes\n");

			for (int i = 0; i < populationSize; ++i)
				fprintf(output, "CHR#%d\t\t", i + 1);
			fprintf(output, "\n");
			for (int i = 0; i < networkNumEdges; ++i)
			{
				for (int j = 0; j < populationSize; ++j)
				{
					fprintf(output, "%d\t", chromosomes[j].edgeIDS[i]);
					fprintf(output, "\t\t");
				}
				fprintf(output, "\n");
			}

			if ((output = fopen(outputFile, "a")) == NULL)
				fatal("Unable to open output file");

			fclose(output);
		}
	}
}

GA::~GA()
{
	delete[] chromosomesBitArr;
	delete[] chromosomes;
	delete[] originalEdgeIDS;
	delete[] simpleGAChromosome;
	delete[] ndmap;
}

void GA::getOriginalEdgeIDS()
{
	for (int i = 0; i < networkNumEdges; ++i)
	{
		std::cout << originalEdgeIDS[i] << "\t";
	}
	std::cout << std::endl;
}

// taken from http://preshing.com/20121224/how-to-generate-a-sequence-of-unique-random-integers/
unsigned int permuteQPR(unsigned int x)
{
	static const unsigned int prime = 4294967291;
	if (x >= prime)
		return x; // The 5 integers out of range are mapped to themselves.
	unsigned int residue = ((unsigned long long)x * x) % prime;
	return (x <= prime / 2) ? residue : prime - residue;
}

// generate random number in range (min, max)
int GA::generateRandomNumber(int min, int max)
{
	return ((rand() % (max - min)) + min);
}

double GA::averageFitnessForPopulation()
{
	double total_fitness = 0.0;
	for (int i = 0; i < populationSize; ++i)
	{
		//chromosomes[i].calculateFitness();
		total_fitness += calculateFitness(i);
	}
	return (total_fitness / populationSize);
}

// params : chromosomeIndex, edgeIndex, depth(0-prev, 1-next)
int GA::get_bit(int chrIndex, int pos, int popState)
{
	int return_val = (1 & chromosomesBitArr[chrIndex][pos / 8][popState] >> (pos - (8 * (pos / 8))));
	return return_val;
}

// toggle the bit at chromosomeIndex, edgeIndex, depth(0-prev, 1-next)
void GA::toggle_bit(int chrIndex, int pos, int popState)
{
	chromosomesBitArr[chrIndex][pos / 8][popState] ^= 1 << (pos - (8 * (pos / 8)));
}

// sets bit at chromosomeIndex, edgeIndex, depth(0-prev, 1-next)
void GA::set_bit(int chrIndex, int pos, int popState)
{
	chromosomesBitArr[chrIndex][pos / 8][popState] |= 1 << (pos - (8 * (pos / 8)));
}

// unsets bit at chromosomeIndex, edgeIndex, depth(0-prev, 1-next)
void GA::unset_bit(int chrIndex, int pos, int popState)
{
	chromosomesBitArr[chrIndex][pos / 8][popState] &= ~(1 << (pos - (8 * (pos / 8))));
}

void GA::initializeRates()
{

	crossover_rate = GA_CROSSOVER_RATE;
	mutation_rate = GA_MUTATION_RATE;
	reproduction_rate = GA_REPRODUCTION_RATE;
	numGenomeMutations = round(mutation_rate * networkNumEdges);
	numGenesRepairSize = round(GA_GENE_REPAIR_PERCENT * networkNumVertices);

	// double total_rates_sum = crossover_rate + mutation_rate + reproduction_rate;
	// double delta = (1 - total_rates_sum) / 3;
	// crossover_rate += delta;
	// mutation_rate += delta;
	// reproduction_rate += delta;

	if (GA_DEBUG)
	{
		std::cout << "- - - GA RATES - - - " << std::endl;
		std::cout << "NORMALIZED CROSSOVER RATE : " << crossover_rate << std::endl;
		std::cout << "NORMALIZED MUTATION RATE : " << mutation_rate << std::endl;
		std::cout << "NORMALIZED REPRODUCTION RATE : " << reproduction_rate << std::endl;
		std::cout << "TOURNAMENT SIZE : " << GA_TOURNAMENT_SIZE << std::endl;
		std::cout << "NUM MUTATION SIZE : " << numGenomeMutations << std::endl;
		std::cout << "NUM GENE REPAIR SIZE : " << numGenesRepairSize << std::endl;
		std::cout << "- - - - - - - - - - -" << std::endl;
	}
}

void GA::mutate(int chromosomeIndex, int popState)
{
	int mutation_site = -1;
	// int curMutations = -1;
	// int rchoice = generateRandomNumber(0, 2);
	// if ( rchoice )
	// 	curMutations = generateRandomNumber(numGenomeMutations/2, numGenomeMutations);
	// else
	// 	curMutations = generateRandomNumber(0 ,numGenomeMutations / 2);
	for (int i = 0; i < networkNumEdges; i++)
	{
		int rchoice = generateRandomNumber(1, 100);
		if (rchoice < (int)(GA_MUTATION_RATE * 100))
			toggle_bit(chromosomeIndex, i, popState);
	}
	if (GA_DEBUG_L2)
	{
		std::cout << "MUTATION SITE FOR CHR# " << chromosomeIndex << " -> " << mutation_site << std::endl;
	}
}

void swap(int &a, int &b)
{
	int temp = a;
	a = b;
	b = temp;
}

void GA::printPopData(int depth)
{
	std::fstream file; // declare an object of fstream class
	char filename[256] = {0};
	strcpy(filename, GA_POP_FILE.c_str());
	strcat(filename, dataset_name);
	strcat(filename, the_date);
	file.open(filename, std::ios::out | std::ios::app); // open file in append mode

	for (int i = 0; i < populationSize; ++i)
	{
		for (int j = 0; j < networkNumEdges; ++j)
		{
			file << get_bit(i, j, depth) << "\t";
		}
		file << "(" << calculateFitness(i, depth) << ")";
		file << std::endl;
	}
	file << std::endl;
	file << "----------------------";
	file << std::endl;
	file.close();
}

void GA::printChromosome(int chr_index, int depth)
{
	std::fstream file; // declare an object of fstream class
	char filename[256] = {0};
	strcpy(filename, GA_BST_FILE.c_str());
	strcat(filename, dataset_name);
	strcat(filename, the_date);
	file.open(filename, std::ios::out | std::ios::app); // open file in append mode
	file << "BEST CHR INDEX : " << chr_index << std::endl;
	for (int j = 0; j < networkNumEdges; ++j)
	{
		file << get_bit(chr_index, j, depth) << "\t";
	}
	file << "(" << calculateFitness(chr_index, depth) << ")";
	file << std::endl;

	file << std::endl;
	file << "----------------------";
	file << std::endl;
	file.close();
}

void GA::move_chromosome_to_next_gen(int chromosomeIndex, int populationIndex)
{
	for (int i = 0; i < networkNumEdges; ++i)
	{
		if (get_bit(chromosomeIndex, i, prev))
			set_bit(populationIndex, i, next);
		else
			unset_bit(populationIndex, i, next);
	}
}

void GA::set_data_name(char *name)
{
	dataset_name = (char *)malloc(sizeof(name));
	strcpy(dataset_name, name);
}

int nearestEvenInt(int to)
{
	return (to % 2 == 0) ? to : (to + 1);
}

void GA::repairGene(int chromosomeIndex, int depth)
{

	if (GA_DEBUG_L2)
	{
		std::cout << "Before Repair : \n";
		for (int n = 0; n < networkNumEdges; n++)
		{
			std::cout << get_bit(chromosomeIndex, n, depth) << "\t";
		}
		std::cout << "\n";
	}

	for (int p = 0; p < networkNumEdges; p++)
	{
		addEdgeByEdgeID(originalEdgeIDS[p]);
	}

	for (int ni = 0; ni < numGenesRepairSize; ni++)
	{
		// std :: cout << "\n\t------rep perc for node # " << ndmap[ni].node_index << " : " << repairPercent << "\n";
		// std ::cout << "\n Cdegree : " << gaSparseNetwork->getDegree(ndmap[ni].node_index) << " , Odegree : " << ndmap[ni].degree << "\n\n";
		int nodeEndPoints[ndmap[ni].degree];
		gaSparseNetwork->getNodeEndPoints(ndmap[ni].node_index, *nodeEndPoints);
		chromosome_g2p(chromosomeIndex, depth);
		double repairPercent = 1 - (double)gaSparseNetwork->getDegree(ndmap[ni].node_index) / ndmap[ni].degree;
		// std ::cout << "\n\t------repair percent :" << repairPercent << "  \n\n";
		for (int p = 0; p < networkNumEdges; p++)
			addEdgeByEdgeID(originalEdgeIDS[p]);

		for (int pts = 0; pts < ndmap[ni].degree; pts++)
		{
			long long int start = (ndmap[ni].node_index < nodeEndPoints[pts]) ? ndmap[ni].node_index : nodeEndPoints[pts];
			long long int end = (ndmap[ni].node_index < nodeEndPoints[pts]) ? nodeEndPoints[pts] : ndmap[ni].node_index;
			long long int edgeID = (networkNumVertices * start) + end;
			long long int position = getEdgeIDIndex(edgeID);
			// std ::cout << "\n\t edgeID : " << edgeID << " , position " << position << " , " << start << " , " << end << "\n";
			if (get_bit(chromosomeIndex, position, depth))
			{
				// std ::cout << "\n\t edge : " << edgeID << " , " << position << " added\n";
				int r = generateRandomNumber(0, 100);
				if (r < GA_GENE_REPAIR_CHANCE * 100)
				{
					unset_bit(chromosomeIndex, position, depth);
				}
			}
		}
	}

	if (GA_DEBUG_L2)
	{
		std::cout << "After Repair: \n";
		for (int n = 0; n < networkNumEdges; n++)
		{
			std::cout << get_bit(chromosomeIndex, n, depth) << "\t";
		}
		std::cout << "\n";
	}
}

// max EdgeID can be (networkNumVertices)^2 for generating random edgeIDs
// num of edges removed can be a max upto (2, networkNumEdges/2)
void GA::generate_GA()
{
	int crossover_discards = 0, numCrossovers = 0, numMutations = 0;

	std::fstream file, pop_file; // declare an object of fstream class
	int currentGeneration = 0;
	int minCrossoverSize = round(GA_CROSSOVER_SIZE * networkNumEdges);
	// int minCrossoverSize = 0;
	int mutation_pop_size = (int)(GA_MUTATION_RATE * populationSize);
	int numEliteChromosomes = nearestEvenInt((int)(GA_REPRODUCTION_RATE * populationSize));

	char filename[256] = {0};
	strcpy(filename, GA_LOG_FILE.c_str());
	strcat(filename, dataset_name);
	strcat(filename, the_date);
	file.open(filename, std::ios::out | std::ios::app); // open file in append mode

	file << "------ GA RUN PARAMS ------" << std::endl;
	file << "DATASET NAME : " << dataset_name << std::endl;
	file << "GENERATIONS : " << numGenerations << std::endl;
	file << "POPULATION SIZE : " << populationSize << std::endl;
	file << "GRAPH NUM VERTICES : " << networkNumVertices << std::endl;
	file << "GRAPH NUM EDGES : " << networkNumEdges << std::endl;
	file << "MINIMUM CROSSOVER SIZE PERCENT : " << GA_CROSSOVER_SIZE << std::endl;
	file << "ELITISM (INDIVIDUALS) : " << numEliteChromosomes << std::endl;
	file << "MIN CROSSOVER SIZE (INDIVIDUALS) : " << minCrossoverSize << std::endl;
	file << "FITNESS : " << ((GA_FITNESS_MODULARITY == 1) ? "MODULARITY" : "Qs") << std::endl;
	file << "------ GA RUN PARAMS ------" << std::endl;

	memset(filename, 0, sizeof(filename));

	strcpy(filename, GA_POP_FILE.c_str());
	strcat(filename, dataset_name);
	strcat(filename, the_date);
	pop_file.open(filename, std::ios::out | std::ios::app); // open file in append mode
	pop_file << "originalEdgeIDS : " << std::endl;
	for (int i = 0; i < networkNumEdges; ++i)
	{
		pop_file << originalEdgeIDS[i] << "\t";
	}
	pop_file << std::endl
			 << std::endl;
	pop_file.close();

	// printPopData(0);
	double threshold_fitness = 1.0 - GA_TOL;
	double max_fitness_generation = 0.0;
	int max_fitness_chr = -1;
	chromosome_map cmap[populationSize], tcmap[GA_TOURNAMENT_SIZE];

	// print the chromosome array
	// printChromosomes(0);

	if (GA_DEBUG)
	{

		double max_fitness = calculateFitness(0, prev);
		int max_fitness_chr = 0, i = 0;
		double total_fitness = max_fitness;
		for (i = 1; i < populationSize; ++i)
		{
			int numSingletons = 0;
			total_fitness += calculateFitness(i, prev);
			
			// print the clustering scheme
			// std::cout << "\n\n------Clustering Scheme for Chromosome : " << i + 1 << std::endl;
			// for (int j = 0; j < networkNumVertices; j++)
			// {
			// 	std::cout << gaSparseNetwork->globalClusterNum[j] << "\t";
			// 	if (gaSparseNetwork->globalClusterNum[j] == -1)
			// 		numSingletons += 1;
			// }
			// std::cout << "\n\n Number of clusters : " << *(gaSparseNetwork->kcluster) << std::endl;
			// std::cout << "\n\n Singletons : " << numSingletons << std::endl;

			std::cout << std::endl;

			if (max_fitness < chromosomes[i].getFitness())
			{
				max_fitness_chr = i;
				max_fitness = chromosomes[i].getFitness();
			}
		}
		file << "averageFitnessForGen #" << currentGeneration << " : " << (double)total_fitness / populationSize << std::endl;
		file << "Best Chr# " << max_fitness_chr << " -> Fitness : " << max_fitness << std::endl;
		printChromosome(max_fitness_chr, next);
		max_fitness_generation = max_fitness;

		std::fstream test_file; // declare an object of fstream class
		char filename[256] = {0};
		strcpy(filename, GA_BST_AVG_FITNESS_RUN.c_str());
		strcat(filename, dataset_name);
		strcat(filename, the_date);
		test_file.open(filename, std::ios::out | std::ios::app); // open file in append mode
		test_file << max_fitness_generation << "\t" << (double)total_fitness / populationSize << std::endl;
		test_file.close();
	}

	while (currentGeneration < numGenerations && max_fitness_generation < threshold_fitness)
	{

		std::cout << "\n---GENERATION #" << currentGeneration + 1 << "---" << std::endl;

		int populationIndex = 0;

		// skip first generation if needed
		if (currentGeneration >= 0)
		{
			// ELITISM
			for (int i = 0; i < populationSize; ++i)
			{
				cmap[i].chromosome_index = i;
				cmap[i].fitness = calculateFitness(i, prev);
			}
			std::sort(cmap, cmap + populationSize, &chromosome_sorter);

			// move elite chromosomes
			for (int i = 0; i < numEliteChromosomes; ++i)
			{
				move_chromosome_to_next_gen(cmap[i].chromosome_index, populationIndex++);
			}
		}

		while (populationIndex < populationSize)
		{

			// tournament selection

			int parentsForCrossover[2]; // 2 parents for crossover
			int nextGenChromosomeState[populationSize];

			for (int ip = 0; ip < 2; ip++)
			{
				for (int i = 0; i < populationSize; ++i)
					nextGenChromosomeState[i] = 0; // no chromosome is selected 0, if selected 1

				// bias the selection to have atleast one of the elite chromosomes
				int bias_chr_index = generateRandomNumber(0, numEliteChromosomes);
				nextGenChromosomeState[bias_chr_index] = 1;
				tcmap[0].chromosome_index = bias_chr_index;
				tcmap[0].fitness = calculateFitness(bias_chr_index, prev);

				for (int i = 1; i < GA_TOURNAMENT_SIZE; ++i)
				{
					// srand (time(NULL));
					int index = generateRandomNumber(0, populationSize);
					while (nextGenChromosomeState[index])
					{
						index = generateRandomNumber(0, populationSize);
					}

					nextGenChromosomeState[index] = 1;

					tcmap[i].chromosome_index = index;
					tcmap[i].fitness = calculateFitness(index, prev);
				}

				std::sort(tcmap, tcmap + GA_TOURNAMENT_SIZE, &chromosome_sorter);
				parentsForCrossover[ip] = tcmap[0].chromosome_index;
			}

			if (GA_DEBUG_L2)
			{
				for (int i = 0; i < populationSize; ++i)
				{
					std::cout << nextGenChromosomeState[i] << "\t";
					if ((i + 1) % 5 == 0 && i > 0)
					{
						std::cout << std::endl;
					}
				}
				std::cout << std::endl;

				file << "TOURNAMENT SELECTION (nextGenChromosomeState) : " << std::endl;
				for (int i = 0; i < populationSize; ++i)
				{
					file << nextGenChromosomeState[i] << "\t";
				}
				file << std::endl;

				std::cout << "\n---TOURNAMENT SORT : \n";
				for (int i = 0; i < GA_TOURNAMENT_SIZE; i++)
				{
					std::cout << tcmap[i].chromosome_index << " , " << tcmap[i].fitness << "\n";
				}
				std::cout << "\n---TOURNAMENT SORT END \n";
			}

			if (GA_DEBUG_L2)
				for (int i = 0; i < 2; ++i)
				{
					std::cout << "parentsForCrossover #" << i << " :-> " << parentsForCrossover[i] << std::endl;
					std::cout << "parentsForCrossoverFitness #" << i << " :-> " << chromosomes[parentsForCrossover[i]].getFitness() << std::endl;

					file << "parentsForCrossover #" << i << " :-> " << parentsForCrossover[i] << std::endl;
					file << "parentsForCrossoverFitness #" << i << " :-> " << chromosomes[parentsForCrossover[i]].getFitness() << std::endl;
				}

			// CROSSOVER

			// int r2 = generateRandomNumber(1, 101);
			// if ( r2 < (GA_CROSSOVER_RATE*100) ) {
			if (1)
			{

				int randomBitVector[networkNumEdges];
				int populationCounter = 0;

				if (GA_DEBUG_L2)
				{
					for (int i = 0; i < 2; ++i)
					{
						std::cout << "PARENTS # " << parentsForCrossover[i] << std::endl;
						for (int j = 0; j < networkNumEdges; ++j)
						{
							std::cout << get_bit(parentsForCrossover[i], j, 0) << std::endl;
						}
						std::cout << "---" << std::endl;
					}
				}

				// int numCrossoverSites = generateRandomNumber(1, networkNumEdges);
				// int numCrossoverSites = generateRandomNumber(minCrossoverSize, networkNumEdges);
				// int numCrossoverSites = generateRandomNumber(0, minCrossoverSize);
				// int numCrossoverSites = minCrossoverSize;
				int numCrossoverSites = -1;
				int rchoice = generateRandomNumber(0, 2);
				if (rchoice)
					numCrossoverSites = generateRandomNumber(1, networkNumEdges / 2);
				else
					numCrossoverSites = generateRandomNumber(networkNumEdges / 2, networkNumEdges);

				if (GA_DEBUG_L2)
					std::cout << "\t-----CROSSOVER SITES : " << numCrossoverSites << "\n";

				if (GA_ONE_POINT_CROSSOVER)
				{
					int firstPoint = generateRandomNumber(1, networkNumEdges / 2);
					int secondPoint = generateRandomNumber(networkNumEdges / 2, networkNumEdges);
					if (GA_DEBUG)
					{
						std::cout << "\n\t-----ONE WAY PT : " << firstPoint << "\n";
						std::cout << "\n\t-----TWO WAY PT : " << secondPoint << "\n";
					}
					for (int i = 0; i < networkNumEdges; i++)
					{
						// perform crossover at sites
						//if (i <= firstPoint && i >= secondPoint)
						if (i >= numCrossoverSites)
						{
							// swap values at i
							if (get_bit(parentsForCrossover[1], i, prev))
								set_bit(populationIndex, i, next);
							else
								unset_bit(populationIndex, i, next);

							if (get_bit(parentsForCrossover[0], i, prev))
								set_bit(populationIndex + 1, i, next);
							else
								unset_bit(populationIndex + 1, i, next);
						}
						else
						{
							if (get_bit(parentsForCrossover[0], i, prev))
								set_bit(populationIndex, i, next);
							else
								unset_bit(populationIndex, i, next);

							if (get_bit(parentsForCrossover[1], i, prev))
								set_bit(populationIndex + 1, i, next);
							else
								unset_bit(populationIndex + 1, i, next);
						}
					}
				}
				else
				{

					for (int i = 0; i < networkNumEdges; ++i)
						randomBitVector[i] = 0;

					while (numCrossoverSites)
					{
						int crossoverSite = generateRandomNumber(0, networkNumEdges);
						while (randomBitVector[crossoverSite])
							crossoverSite = generateRandomNumber(0, networkNumEdges);
						randomBitVector[crossoverSite] = 1;
						numCrossoverSites--;
					}

					for (int i = 0; i < networkNumEdges; ++i)
					{
						if (GA_ONE_WAY_CROSSOVER)
						{
							// perform crossover at sites
							if (randomBitVector[i])
							{
								// swap values at i
								if (get_bit(parentsForCrossover[0], i, prev))
									set_bit(populationIndex, i, next);
								else
									unset_bit(populationIndex, i, next);

								if (get_bit(parentsForCrossover[0], i, prev))
									set_bit(populationIndex + 1, i, next);
								else
									unset_bit(populationIndex + 1, i, next);
							}
							else
							{
								if (get_bit(parentsForCrossover[0], i, prev))
									set_bit(populationIndex, i, next);
								else
									unset_bit(populationIndex, i, next);

								if (get_bit(parentsForCrossover[1], i, prev))
									set_bit(populationIndex + 1, i, next);
								else
									unset_bit(populationIndex + 1, i, next);
							}
						}
						else
						{
							// perform crossover at sites
							if (randomBitVector[i])
							{
								// swap values at i
								if (get_bit(parentsForCrossover[1], i, prev))
									set_bit(populationIndex, i, next);
								else
									unset_bit(populationIndex, i, next);

								if (get_bit(parentsForCrossover[0], i, prev))
									set_bit(populationIndex + 1, i, next);
								else
									unset_bit(populationIndex + 1, i, next);
							}
							else
							{
								if (get_bit(parentsForCrossover[0], i, prev))
									set_bit(populationIndex, i, next);
								else
									unset_bit(populationIndex, i, next);

								if (get_bit(parentsForCrossover[1], i, prev))
									set_bit(populationIndex + 1, i, next);
								else
									unset_bit(populationIndex + 1, i, next);
							}
						}
					}
				}

				if (GA_CROSSOVER_PARENT_CHILD_BEST)
				{
					double p1_fitness = calculateFitness(parentsForCrossover[0], prev);
					double p2_fitness = calculateFitness(parentsForCrossover[1], prev);
					double c1_fitness = calculateFitness(populationIndex, next);
					double c2_fitness = calculateFitness(populationIndex + 1, next);

					if (p1_fitness > c1_fitness)
					{
						for (int i = 0; i < networkNumEdges; ++i)
						{
							if (get_bit(parentsForCrossover[0], i, prev))
								set_bit(populationIndex, i, next);
							else
								unset_bit(populationIndex, i, next);
						}
						if (GA_DEBUG)
							std::cout << "---P1 Retained\n";
						crossover_discards++;
					}
					else
					{
						numCrossovers++;
					}

					if (p2_fitness > c2_fitness)
					{
						for (int i = 0; i < networkNumEdges; ++i)
						{
							if (get_bit(parentsForCrossover[1], i, prev))
								set_bit(populationIndex + 1, i, next);
							else
								unset_bit(populationIndex + 1, i, next);
						}
						if (GA_DEBUG)
							std::cout << "---P2 Retained\n";
						crossover_discards++;
					}
					else
					{
						numCrossovers++;
					}
				}

				if (GA_DEBUG_L2)
				{
					std::cout << "Printing bits : " << std::endl;
					for (int i = 0; i < 2; ++i)
					{
						for (int j = 0; j < networkNumEdges; ++j)
						{
							std::cout << get_bit(i, j, next) << "\t";
						}
						std::cout << std::endl;
					}
				}
			}

			// MUTATION
			mutate(populationIndex, next);
			mutate(populationIndex + 1, next);

			// REPAIR GENE
			repairGene(populationIndex, next);
			repairGene(populationIndex + 1, next);

			populationIndex += 2;
			printf("Finished PopINdex : %d\n", populationIndex);
		}
		printf("Finished Gen # : %d\n", currentGeneration);
		// END OF CROSSOVER AND MUTATION
		// print population data after each run
		// printPopData(next);

		// STATS
		if (GA_DEBUG)
		{
			double max_fitness = calculateFitness(0, next);
			max_fitness_chr = 0;
			int i = 0;
			double total_fitness = max_fitness;
			for (i = 1; i < populationSize; ++i)
			{
				calculateFitness(i, next);
				total_fitness += chromosomes[i].getFitness();
				if (max_fitness < chromosomes[i].getFitness())
				{
					max_fitness_chr = i;
					max_fitness = chromosomes[i].getFitness();
				}
			}
			file << "averageFitnessForGen #" << currentGeneration << " : " << (double)total_fitness / populationSize << std::endl;
			file << "Best Chr# " << max_fitness_chr << " -> Fitness : " << max_fitness << std::endl;
			printChromosome(max_fitness_chr, next);
			
			max_fitness_generation = max_fitness;

			std::fstream test_file; // declare an object of fstream class
			char filename[256] = {0};
			strcpy(filename, GA_BST_AVG_FITNESS_RUN.c_str());
			strcat(filename, dataset_name);
			strcat(filename, the_date);
			test_file.open(filename, std::ios::out | std::ios::app); // open file in append mode
			test_file << max_fitness_generation << "\t" << (double)total_fitness / populationSize << std::endl;
			test_file.close();

			chromosome_g2p(max_fitness_chr, next);

			char new_filename[256] = {0};
			strcpy(new_filename, GA_BST_GML.c_str());
			strcat(new_filename, the_date);
			strcat(new_filename, dataset_name);
			gaSparseNetwork->printEdges(new_filename);
		}

		std::cout << "\n---END OF GENERATION #" << currentGeneration + 1 << "---" << std::endl;

		prev = !prev;
		next = !next;
		currentGeneration++;
	}

	file.close();

	chromosome_g2p(max_fitness_chr, next);
	printChromosomes(next);
	char new_filename[256] = {0};
	strcpy(new_filename, GA_BST_GML.c_str());
	strcat(new_filename, the_date);
	strcat(new_filename, dataset_name);
	gaSparseNetwork->printEdges(new_filename);

	std::cout << "\n CROSSOVER DISCARDS : " << crossover_discards << std::endl;
	std::cout << "\n CROSSOVER USEFUL : " << numCrossovers << std::endl;

	// EVALUATE

	if (GA_DEBUG_L2)
	{
		int bit_arr_pop_size = (populationSize);
		int bit_arr_num_edges = ARRAY_SIZE(networkNumEdges);
		std::cout << "BitArr Dimensions : " << bit_arr_pop_size << " : " << bit_arr_num_edges << std::endl;
		std::cout << "SIZE OF BIT ARRAY : " << sizeof(chromosomesBitArr) << std::endl;
		std::cout << "--- BIT ARRAY ---\n\n";
		for (int i = 0; i < 2; ++i)
		{
			std::cout << "Depth (" << (i + 1) << ") :" << std::endl;
			for (int j = 0; j < bit_arr_pop_size; ++j)
			{
				for (int k = 0; k < bit_arr_num_edges; ++k)
				{
					int index = 0;
					while (index < 8)
					{
						std::cout << (k * 8) + index << " : (" << j << "," << k << "," << i << "," << index << ") -> ";
						std::cout << (1 & (chromosomesBitArr[j][k][i] >> (index))) << "\n";
						index++;
					}
				}
				std::cout << std::endl;
			}
			std::cout << "End of Depth\n\n";
		}
	}
}
