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

	double calc_fitness = gaSparseNetworkChr->q_calc("ga.out");
	fitness = calc_fitness;
	return calc_fitness;
}

// if chr_index = -1 then all else calcu fitness of chromosome[chr_index]
double GA::calculateFitness(int chr_index, int depth) {
	double fitness = 0.0;
	if(chr_index < 0) {
		fatal("Invalid chromosomeindex");
	}
	else {
		chromosome_g2p(chr_index, depth);
		fitness = chromosomes[chr_index].calculateFitness();
	}
	return fitness;
}

void GA::chromosome_g2p(int chromosomeIndex, int depth) {
	for (int i = 0; i < networkNumEdges; ++i)
	{
		if( get_bit ( chromosomeIndex, i, depth) )
			removeEdgeByID( originalEdgeIDS[i] );
		else
			addEdgeByEdgeID( originalEdgeIDS[i] );
	}
}

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

// toggle bit at index in array
void GA::toggle_bit(char *array, int index) {
    array[index / 8] ^= 1 << (index % 8);
}

GA::GA(Network &sparseNetwork, int popSize, int generations, int numNodes, int numEdges) {

	// set class data memebers
	populationSize = popSize;
    networkNumVertices = numNodes;
    networkNumEdges = numEdges;
    gaSparseNetwork = &sparseNetwork;
    numGenerations = generations;

    initializeRates();

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

    // allocate space for bit locations
    // show edge state (0=removed, 1=in_network)		
    int bit_arr_pop_size = (popSize);
    int bit_arr_num_edges = ARRAY_SIZE(numEdges);

    if ((chromosomesBitArr = new char**[bit_arr_pop_size]) == NULL) // allocate memory to set of Basic GA Chromosomes
    	fatal("memory not allocated");

    for (int i = 0; i < bit_arr_pop_size; ++i)
    {
    	if ((chromosomesBitArr[i] = new char*[bit_arr_num_edges]) == NULL) // allocate memory to set of Basic GA Chromosomes
    		fatal("memory not allocated");
    }

    for (int i = 0; i < bit_arr_pop_size; ++i)
    {
    	for (int j = 0; j < bit_arr_num_edges; ++j)
    	{
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
	    	for (int k = 0; k < networkNumEdges; ++k) {
	    		std::cout << edgeIDState[k] << "\t";
		    }
	    	std::cout << "\n-=-=-=-=-=-=-=-=-=-=-=-=-=-" << std::endl;
    	}

	    
	    int chromosome_edgeID_pos = 0;
	    if(GA_DEBUG)
	    	std::cout << "Removed edges for chromosome : "<< (i+1) << " : \n";

	    // add edges back after generating chromosome
	    // after one chromosome generation
	    for (int k = 0; k < networkNumEdges; ++k){
    		if(edgeIDState[k] == -1){
    			edgeIDState[k] = addEdgeByEdgeID(originalEdgeIDS[k]);
    			chromosomes[i].edgeIDS[chromosome_edgeID_pos++] = originalEdgeIDS[k];
    			set_bit(i, k, 0);
    		}
	    }
	    // assign length to the chromosome
	    chromosomes[i].length = chromosome_edgeID_pos;

	    if(GA_DEBUG)
	    	std::cout << "\n";

    }

    if(GA_DEBUG_L2) {
    	std::cout << "BitArr Dimensions : " << bit_arr_pop_size << " : " << bit_arr_num_edges << std::endl;
    	std::cout << "SIZE OF BIT ARRAY : " << sizeof(chromosomesBitArr) << std::endl;
    	std::cout << "--- BIT ARRAY ---\n\n";
    	for (int i = 0; i < 2; ++i)
	    {
	    	std::cout << "Depth (" << (i+1) << ") :" << std::endl;
	    	for (int j = 0; j < bit_arr_pop_size; ++j)
		    {
		    	for (int k = 0; k < bit_arr_num_edges; ++k)
		    	{
		    		int index = 0;
		    		while(index < 8){
		    			std::cout << (k*8)+index << " : (" << j << "," <<  k << "," << i << "," << index << ") -> ";
		    			std::cout << (1 & (chromosomesBitArr[j][k][i] >> (index))) << "\n";
		    			index++;
		    		}
		    	}
		    	std::cout << std::endl;
		    }
		    std::cout << "End of Depth\n\n";		
	    }
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
	delete [] chromosomesBitArr;
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
	
	return (((rand()) % max) + min); 	
}

double GA::averageFitnessForPopulation() {
	double total_fitness = 0.0;
	for (int i = 0; i < populationSize; ++i)
	{
		//chromosomes[i].calculateFitness();
		total_fitness += calculateFitness(i);
	}
	return (total_fitness/populationSize);
}

// params : chromosomeIndex, edgeIndex, depth(0-prev, 1-next)
int GA::get_bit(int chrIndex, int pos, int popState) {
	int return_val = (1 & chromosomesBitArr[chrIndex][pos/8][popState] >> (pos - (8 * (pos/8))) );
	return return_val;
}

// toggle the bit at chromosomeIndex, edgeIndex, depth(0-prev, 1-next)
void GA::toggle_bit(int chrIndex, int pos, int popState) {
    chromosomesBitArr[chrIndex][pos/8][popState] ^= 1 << (pos - (8 * (pos/8)));
}

// sets bit at chromosomeIndex, edgeIndex, depth(0-prev, 1-next)
void GA::set_bit(int chrIndex, int pos, int popState) {
    chromosomesBitArr[chrIndex][pos/8][popState] |= 1 << (pos - (8 * (pos/8)));
}

// unsets bit at chromosomeIndex, edgeIndex, depth(0-prev, 1-next)
void GA::unset_bit(int chrIndex, int pos, int popState) {
    chromosomesBitArr[chrIndex][pos/8][popState] &= ~(1 << (pos - (8 * (pos/8))));
}

void GA::initializeRates() {

	crossover_rate = GA_CROSSOVER_RATE;
	mutation_rate = GA_MUTATION_RATE;
	reproduction_rate = GA_REPRODUCTION_RATE;

	double total_rates_sum = crossover_rate + mutation_rate + reproduction_rate;
	double delta = (1 - total_rates_sum) / 3;
	crossover_rate += delta;
	mutation_rate += delta;
	reproduction_rate += delta;

	if(GA_DEBUG) {
		std::cout << "- - - GA RATES - - - " << std::endl;
		std::cout << "NORMALIZED CROSSOVER RATE : " << crossover_rate << std::endl;
		std::cout << "NORMALIZED MUTATION RATE : " << mutation_rate << std::endl;
		std::cout << "NORMALIZED REPRODUCTION RATE : " << reproduction_rate << std::endl;
		std::cout << "- - - - - - - - - - -" << std::endl;
	}

}

void GA::mutate(int chromosomeIndex, int popState) {
	int mutation_site = generateRandomNumber(0, networkNumEdges);
	toggle_bit( chromosomeIndex, mutation_site, popState);

	if(GA_DEBUG_L2) {
		std::cout << "MUTATION SITE FOR CHR# " << chromosomeIndex << " -> " << mutation_site << std::endl;
	}
}

void swap(int &a, int &b) {
	int temp = a;
	a = b;
	b = temp;
}

void GA::printPopData(int depth) {
	std :: fstream file; // declare an object of fstream class
    file.open(GA_POP_FILE.c_str(), std :: ios :: out | std :: ios :: app); // open file in append mode
   	
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

void GA::printChromosome(int chr_index, int depth) {
	std :: fstream file; // declare an object of fstream class
    file.open(GA_BST_FILE.c_str(), std :: ios :: out | std :: ios :: app); // open file in append mode
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

// max EdgeID can be (networkNumVertices)^2 for generating random edgeIDs
// num of edges removed can be a max upto (2, networkNumEdges/2)
void GA::generate_GA() {

	if( remove( GA_LOG_FILE.c_str() ) != 0 )
    	perror( "Error deleting file" );

    if( remove( GA_POP_FILE.c_str() ) != 0 )
    	perror( "Error deleting file" );

    if( remove( GA_BST_FILE.c_str() ) != 0 )
    	perror( "Error deleting file" );

	std :: fstream file, pop_file; // declare an object of fstream class
    file.open(GA_LOG_FILE.c_str(), std :: ios :: out | std :: ios :: app); // open file in append mode

    pop_file.open(GA_POP_FILE.c_str(), std :: ios :: out | std :: ios :: app); // open file in append mode
    pop_file << "originalEdgeIDS : " << std::endl;
    for (int i = 0; i < networkNumEdges; ++i)
    {
    	pop_file << originalEdgeIDS[i] << "\t";
    }
    pop_file << std::endl << std::endl;
    pop_file.close();

	int currentGeneration = 0;
	int minCrossoverSize = (int) GA_CROSSOVER_SIZE_RATE * networkNumEdges;
	int mutation_pop_size = (int)(mutation_rate * populationSize);

	if(GA_DEBUG) {
		double total_fitness = 0.0;
		for (int i = 0; i < populationSize; ++i) {
			calculateFitness(i, prev);
			total_fitness += chromosomes[i].getFitness();
		}
		file << "averageFitnessForGen #0 : " << total_fitness/populationSize << std::endl;
	}

	printPopData(0);
	double threshold_fitness = 0.95;
	double max_fitness_generation = 0.0;

	while(currentGeneration < numGenerations && max_fitness_generation <= threshold_fitness) {

		int populationIndex = 0;

		while( populationIndex != populationSize ) {

			// tournament selection

			int parentsForCrossover[2]; // 2 parents for crossover
			int nextGenChromosomeState[populationSize];
				
			for (int i = 0; i < populationSize; ++i)
				nextGenChromosomeState[i] = 0; // no chromosome is selected 0, if selected 1

			int start_index = 999999999;
			for (int i = 0; i < GA_TOURNAMENT_SIZE; ++i) {
				// srand (time(NULL));
				int index = generateRandomNumber(0, populationSize);
				while(nextGenChromosomeState[index]) {
					index = generateRandomNumber(0, populationSize);
				}
				if(start_index>index)
					start_index = index;
				nextGenChromosomeState[index] = 1;
			}

			if(GA_DEBUG) {
				for (int i = 0; i < populationSize; ++i) {
					std::cout << nextGenChromosomeState[i] << "\t";
					if ((i+1)%5==0 && i>0)
					{
						std::cout << std::endl;
					}
				}
				std::cout << std::endl;
				std::cout << "Min Index : " << start_index << std::endl;

				file << "TOURNAMENT SELECTION (nextGenChromosomeState) : " << std::endl;
				for (int i = 0; i < populationSize; ++i)
				{
					file << nextGenChromosomeState[i] << "\t";
				}
				file << std::endl;
			}

			// start_index helps to jump to the first minimum to skip few iterations
			double max_fitness = 0.00 - GA_TOL;
			int parent_chr_1 = -1;
			int parent_chr_2 = -1;

			int next_start_index = start_index + 1;
			while( nextGenChromosomeState[next_start_index] != 1 )
				next_start_index += 1;

			double f_c1 = calculateFitness(start_index, prev);
			double f_c2 = calculateFitness(next_start_index, prev);

			parent_chr_1 = start_index;
			parent_chr_2 = next_start_index;
			// make the first chr to be the highest
			if(f_c1 < f_c2) {
				swap(parent_chr_1, parent_chr_2);
				double f_temp = f_c1;
				f_c1 = f_c2;
				f_c2 = f_temp;
			}

			int iter = (parent_chr_1 < parent_chr_2) ? parent_chr_2 : parent_chr_1;

			for (int i = iter+1; i < populationSize; i++) {
				if(nextGenChromosomeState[i]) {
					double chr_fitness = calculateFitness(i, prev);
					if(chr_fitness > f_c1) {
						parent_chr_2 = parent_chr_1;
						f_c2 = f_c1;
						parent_chr_1 = i;
						f_c1 = chr_fitness;
					} else if ( chr_fitness > f_c2 ) {
						parent_chr_2 = i;
						f_c2 = chr_fitness;
					}
				}
			}

			parentsForCrossover[0] = parent_chr_1;
			parentsForCrossover[1] = parent_chr_2;
			

			if( GA_DEBUG )
				for (int i = 0; i < 2; ++i)
				{
					std::cout << "parentsForCrossover #" << i << " :-> " << parentsForCrossover[i] << std::endl;
					std::cout << "parentsForCrossoverFitness #" << i << " :-> " << chromosomes[parentsForCrossover[i]].getFitness() << std::endl;

					file << "parentsForCrossover #" << i << " :-> " << parentsForCrossover[i] << std::endl;
					file << "parentsForCrossoverFitness #" << i << " :-> " << chromosomes[parentsForCrossover[i]].getFitness() << std::endl;
				}

			// CROSSOVER
			
			int randomBitVector[networkNumEdges];
			int populationCounter = 0;

			if(GA_DEBUG_L2) {
				for (int i = 0; i < 2; ++i)
				{
					std::cout << "PARENTS # " << parentsForCrossover[i] << std::endl;
					for (int j = 0; j < networkNumEdges; ++j)
					{
						std::cout << get_bit( parentsForCrossover[i], j, 0 ) << std::endl;
					}
					std::cout << "---" << std::endl;
				}
			}

			int numCrossoverSites = generateRandomNumber(minCrossoverSize, networkNumEdges - 1);

			for (int i = 0; i < networkNumEdges; ++i)
				randomBitVector[i] = 0;

			// while(numCrossoverSites) {
			// 	int crossoverSite = generateRandomNumber(0, networkNumEdges);
			// 	while(randomBitVector[crossoverSite])
			// 		crossoverSite = generateRandomNumber(0, networkNumEdges);
			// 	randomBitVector[crossoverSite] = 1;
			// 	numCrossoverSites--;
			// }

			for (int i = 0; i < networkNumEdges; ++i)
			{
				// perform crossover at sites
				if(randomBitVector[i]) {
					// swap values at i
					if( get_bit( parentsForCrossover[1], i, prev ) )
						set_bit ( populationIndex, i, next );
					else
						unset_bit( populationIndex, i, next );

					if( get_bit( parentsForCrossover[0], i, prev ) )
						set_bit ( populationIndex+1, i, next );
					else
						unset_bit( populationIndex+1, i, next );

				} else {
					if( get_bit( parentsForCrossover[0], i, prev ) )
						set_bit ( populationIndex, i, next );
					else
						unset_bit( populationIndex, i, next );

					if( get_bit( parentsForCrossover[1], i, prev ) )
						set_bit ( populationIndex+1, i, next );
					else
						unset_bit( populationIndex+1, i, next );
				}
			}

			if(GA_DEBUG_L2) {
				std::cout << "Printing bits : " << std::endl;
		    	for (int i = 0; i < 2; ++i)
		    	{
		    		for (int j = 0; j < networkNumEdges; ++j)
					{
						std::cout << get_bit( i, j, next ) << "\t";
					}
					std::cout << std::endl;	
		    	}
		    }

		    populationIndex += 2;

		    // MUTATE SOME INDIVIDUAL
		    int mutated_chr = 0;
			while(mutated_chr < mutation_pop_size) {
				int mutation_chr = generateRandomNumber(0, populationIndex);
				mutate(mutation_chr, next);
				mutated_chr++;
				if(GA_DEBUG_L2)
					std::cout << "MUTATION CHR#: " << mutation_chr << std::endl;

			}
		    
		}

		
	    

		
		// END OF CROSSOVER AND MUTATION
		// print population data after each run
		printPopData(next);
		
		// STATS
		if(GA_DEBUG) {
			double max_fitness = 0.0;
			int max_fitness_chr = -1, i = 0;
			double total_fitness = 0.0;
			for (i = 0; i < populationSize; ++i) {
				calculateFitness(i, next);
				total_fitness += chromosomes[i].getFitness();
				if( max_fitness < chromosomes[i].getFitness() ) {
					max_fitness_chr = i;
					max_fitness = chromosomes[i].getFitness();
				}
			}
			file << "averageFitnessForGen #" << currentGeneration << " : " << (double)total_fitness/populationSize << std::endl;
			file << "Best Chr# " << max_fitness_chr << " -> Fitness : " << max_fitness << std::endl;
			printChromosome(max_fitness_chr, next);
			max_fitness_generation = max_fitness;

			std :: fstream test_file; // declare an object of fstream class
    		test_file.open("ga_bst_test.log", std :: ios :: out | std :: ios :: app); // open file in append mode
    		test_file << max_fitness_generation << std::endl;
    		test_file.close();
		}

		prev = !prev;
		next = !next;
		currentGeneration++;

	}

	file.close();

	// EVALUATE

	if(GA_DEBUG_L2) {
		double total_fitness = 0.0;
		for (int i = 0; i < populationSize; ++i) {
			std::cout << "Fitness for CHR#" << (i+1) << " :-> " << chromosomes[i].getFitness() << std::endl;
			total_fitness += chromosomes[i].getFitness();
		}
		std::cout << "averageFitnessForPopulation : " << total_fitness/populationSize << std::endl;
	}

	if(GA_DEBUG_L2) {
		int bit_arr_pop_size = (populationSize);
    	int bit_arr_num_edges = ARRAY_SIZE(networkNumEdges);
    	std::cout << "BitArr Dimensions : " << bit_arr_pop_size << " : " << bit_arr_num_edges << std::endl;
    	std::cout << "SIZE OF BIT ARRAY : " << sizeof(chromosomesBitArr) << std::endl;
    	std::cout << "--- BIT ARRAY ---\n\n";
    	for (int i = 0; i < 2; ++i)
	    {
	    	std::cout << "Depth (" << (i+1) << ") :" << std::endl;
	    	for (int j = 0; j < bit_arr_pop_size; ++j)
		    {
		    	for (int k = 0; k < bit_arr_num_edges; ++k)
		    	{
		    		int index = 0;
		    		while(index < 8){
		    			std::cout << (k*8)+index << " : (" << j << "," <<  k << "," << i << "," << index << ") -> ";
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