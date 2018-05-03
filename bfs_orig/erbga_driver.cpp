/****************************************************************************
*
*	erbga_driver.cpp:	Code for exploring a network using breadth-first 
*                       search.  Outputs connected components.
*
*                       Sharlee Climer
*                       October/November 2007
*
*      Modified to allow duplicate edges, count them separately, 
*      and record highest weight from the duplicates.
*      Sharlee Climer, July 2009
*
*
*      Modified by Aditya Karnam for (Efficient Reduced-Bias Genetic 
*      Algorithm for Generic Community Detection Objectives) ERBGA. 
*      July 2017 - March 2018.
*      Modified to assign node mapping inside the Network class. 
*      Added function for running GA. 
*
****************************************************************************/

#include "erbga_driver.h"
#include "network.h"
#include "erbga.h"

char* input_fname = "";
char* output_fname = "";

int helpFlag = 0;

int popSize = 0;
int generations = 0;
float random_pop_rate = GA_RANDOM_POP_PERCENT, 
        elitism_rate = GA_ELITISM_RATE,
        mutation_rate = GA_MUTATION_RATE,
        gene_repair_rate = GA_GENE_REPAIR_CHANCE;
int tournament_size = GA_TOURNAMENT_SIZE, 
        gene_repair_size = GA_GENE_REPAIR_PERCENT, 
        crossover_type = 0;
std::string population_checkpoint_file = "";

void showProgramVersion() {
    printf("ERBGA -- v0.4.\n");
    
    printf("Developed by \n\tAditya Karnam Gururaj Rao\n\tJuly 2017 - May 2018\n");
    printf("Changelog--\n");
    printf("\tv0.1 - Sept 2017 - Added Genetic Operators\n");
    printf("\tv0.2 - Oct 2017 - Added Elitism Operator\n");
    printf("\tv0.3 - One Point/One Way/Uniform Crossover Operators\n");
    printf("\tv0.4 - Added Gene Repair Operator\n");
    printf("\tv0.5 - Refactoring and Command Line args support.\n");

    printf("\tDeveloped on the base code \033[1mNetwork\033[0m Library\n");
    printf("\tShoutout to Sharlee Climer for her guidance.\n\n\n");
    
}
void showHelpMessage() {
    showProgramVersion();
	  printf("--help: Prints this help message.\n");
    printf("--popsize: Allows you to set the population size of the GA.\n");
    printf("\tThis is a required argument.\n");
    printf("--gen: Set the number of epochs for the GA.\n");
    printf("\tThis is a required argument.\n");
    printf("--in: Input GML filename.\n");
    printf("\tThis is a required argument.\n");
    printf("--out: BFS Output filename.\n");
    printf("\tThis is a required argument.\n");
    printf("--rpop: Set the random population rate.\n");
    printf("\tOptional. Ranges 0 through 1.\n");
    printf("--tsize: Set the tournament size for breeding population of GA.\n");
    printf("\tOptional. Ranges 0 through Population Size.\n");
    printf("--elirate: Set the Elitism for breeding population of GA.\n");
    printf("\tOptional. Ranges 0 through 1.\n");
    printf("--mutrate: Set the Mutation rate for breeding population of GA.\n");
    printf("\tOptional. Ranges 0 through 1.\n");
    printf("--grepsize: Set the Gene Repair size for breeding population of GA.\n");
    printf("\tOptional. Ranges 0 through Population Size.\n");
    printf("--greprate: Set the Gene Repair rate for breeding population of GA.\n");
    printf("\tOptional. Ranges 0 through 1.\n");
    printf("--crosstype: Set the Crossover Type for breeding population of GA.\n");
    printf("\tOptional. Ranges 0,1 or 2.\n");
    printf("\t\t0. One Point Crossover.\n");
    printf("\t\t1. One Way Crossover.\n");
    printf("\t\t2. Uniform Crossover.\n");
    printf("--popcheckpt: Population Checkpoint file to resume GA.\n");
    printf("\tOptional. Checkpoint filename.\n");
}

// validates the output from the command line args
// defaults to values in the header files
void readCommandLineInputs(int argc, char **argv) {

  int nonoptarg = 0;
  int c;
  opterr = 0;
  int option_index = 0;
  
  while((c = getopt_long (argc, argv, "a:b:c:d:efghijklm", long_options, &option_index)) != -1){
    switch (c)
    {
      case 0:
          /* If this option set a flag, do nothing else now. */
          if (long_options[option_index].flag != 0)
            break;
          printf ("option %s", long_options[option_index].name);
          if (optarg)
            printf (" with arg %s", optarg);
          printf ("\n");
          break;
      case 'a': 
        popSize = atoi(optarg);
        break;
      case 'b': 
        generations = atoi(optarg);
        break;
      case 'c': 
        input_fname = optarg;
        break;
      case 'd': 
        output_fname = optarg;
        break;
      case 'e': 
        random_pop_rate = atof(optarg);
        break;
      case 'f':
        tournament_size = atoi(optarg);
        break;
      case 'g':
        elitism_rate = atof(optarg);
        break;
      case 'h':
        mutation_rate = atof(optarg);
        break;  
      case 'i':
        gene_repair_size = atoi(optarg);
        break;
      case 'j':
        gene_repair_rate = atof(optarg);
        break;
      case 'k':
        crossover_type = atoi(optarg);
        break;
      case 'l':
        population_checkpoint_file = optarg;
        break;
      case 'm':
        // Help section
        helpFlag = 1;
        break;
      case '?':
        if (optopt == 'a' || optopt == 'b' || optopt == 'c' || optopt == 'd')
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr, "Unknown option character `\\x%x'.\n",optopt);
        abort ();
        break;
      default:
        abort ();
    }
  }

  if(!helpFlag) {
    if ( input_fname == "" ) {
      fatal("Option --in requires a string argement. Example filename.gml");
    }
    
    if ( output_fname == "" ) {
      fatal("Option --out requires a string argement. Example out.bfs");
    }

    if ( popSize == 0 ) {
      fatal("Option --popsize requires an integer argement. Example 50");
    }

    if ( generations == 0 ) {
      fatal("Option --gen requires an integer argement. Example 100");
    }
  } else {
    showHelpMessage();
  }

}

int main(int argc, char **argv)
{
  // Randomizer for GA
  srand(time(NULL));
  // if (argc != 3)
  //   fatal("Usage:\n   erbga input.gml output.bfs");

  readCommandLineInputs(argc, argv);

  FILE *input;
  FILE *output;
  timer t;

  t.start("Timer started");

  int testFile = system("test ! -f comp1.gml");
  if (testFile != 0)
    fatal("Component files already exist in this directory");

  if (((input = fopen(input_fname, "r")) == NULL) || ((output = fopen(output_fname, "w")) == NULL))
    fatal("File could not be opened.\n");

  if (DESCRIPTIVE_OUTPUT)
    fprintf(output, "%s\n", input_fname);

  

  long long int min = 10000000000;  // hold min node number
  long long int max = -10000000000; // hold max node number
  char string[50];
  int startOne; // 1 if start node number is 1, 0 if start number is 0

  while (1)
  { // throw away everything until get to "graph" declaration
    if (feof(input))
      fatal("no graph to read in input file");

    fscanf(input, "%s", string); // read in string from file

    if (strncmp(string, "graph", 5) == 0)
      break; // stop when get to "graph" declaration
  }

  // find number of nodes in graph
  long long int numNodes = 0; // number of nodes
  long long int numEdges = 0; // number of edges

  while (1)
  { // throw away everything until get to "edge" declaration
    if (feof(input))
      fatal("no edges to read in input file");

    fscanf(input, "%s", string); // read in string from file

    if (strncmp(string, "id", 2) == 0)
    {
      fscanf(input, "%s", string); // read in id value
      long int num = atoi(string);

      if (min > num)
        min = num; // update minimum node number
      if (max < num)
        max = num; // update maximum node number
      numNodes++;  // one more node found
    }

    if (strncmp(string, "edge", 4) == 0)
      break; // stop when read in first "edge"
  }

  std::cout << "\n Number of nodes read from the GML file : " << numNodes << std::endl;
  std::cout << "\nNode numbers range from " << min << " to " << max << std::endl;
  std::cout << "Reading in graph from " << input_fname << "...\n"
            << std::endl;

  // Create the SparseNetwork
  Network sparseNet(numNodes, DIRECTED, min, max); // if DIRECTED = 0, undirected, min, max

  if (DESCRIPTIVE_OUTPUT)
    fprintf(output, "Node numbers range from %d to %d\n\n", min, max);
  fclose(output);

  startOne = min;
  rewind(input); // start at beginning of input file again

  if (startOne < 0)
    fatal("Node numbers can not be negative");
  if ((startOne != 0) && (startOne != 1))
    warning("Smallest node number should be 0 or 1");

  // record node id numbers
  int ptr = 0; // pointer for filling nodeNumbers array

  while (1)
  { // throw away everything until get to "graph" declaration
    if (feof(input))
      fatal("no graph to read in input file");

    fscanf(input, "%s", string); // read in string from file

    if (strncmp(string, "graph", 5) == 0)
      break; // stop when get to "graph" declaration
  }

  while (1)
  { // record node id numbers
    if (feof(input))
      fatal("no edges to read in input file");

    fscanf(input, "%s", string); // read in string from file

    if (strncmp(string, "id", 2) == 0)
    {
      fscanf(input, "%s", string); // read in id value
      int num = atoi(string);

      if ((num < min) || (num > max))
        fatal("error while reading in node numbers");

      // assign ID to nodeNum
      sparseNet.assignID(ptr++, num);
      //std::cout << ptr << ", " << id[ptr-1]<< std::endl;
    }

    if (strncmp(string, "edge", 4) == 0)
      break; // stop when read in first "edge"
  }

  if (ptr != numNodes)
    fatal("Error reading in node numbers");

  int dupEdges = 0; // record number of duplicate edges

  while (!feof(input))
  {                              // read edges until reach end of file
    fscanf(input, "%s", string); // read in string from file

    if (!strncmp(string, "[", 1) == 0) // check for correct character
      fatal("No '[' after edge declaration");

    fscanf(input, "%s", string);            // read in string from file
    if (!strncmp(string, "source", 6) == 0) // check for correct character
      fatal("No 'source' declaration");

    fscanf(input, "%s", string); // read in start node from file
    int source = atoi(string);

    if ((source < min) || (source > max))
      fatal("Invalid node number");

    fscanf(input, "%s", string);            // read in string from file
    if (!strncmp(string, "target", 6) == 0) // check for correct character
      fatal("No 'target' declaration");

    fscanf(input, "%s", string); // read in target node from file
    int target = atoi(string);

    if ((target < min) || (target > max))
      fatal("Invalid node number");

    fscanf(input, "%s", string); // read in string from file

    float weight = 1.0;
    if (strncmp(string, "weight", 6) == 0)
    {                              // check to see if weight specified
      fscanf(input, "%s", string); // read in string from file
      weight = atof(string);
      fscanf(input, "%s", string); // read in string from file
    }

    // (source, target, weight)
    if (!sparseNet.addEdgeFromGML(source, target, weight))
      dupEdges++;
    else
      numEdges++; // count number of edges

    if (numEdges % 10000000 == 0) // message every 10 million edges
      std::cout << numEdges / 1000000 << " million edges read" << std::endl;

    if (!strncmp(string, "]", 1) == 0) // check for correct character
      fatal("No ']' after edge declaration");

    if (feof(input))
      break; // break is at end of file

    fscanf(input, "%s", string); // read in string from file

    if (strncmp(string, "]", 1) == 0) // check for end of graph
      break;

    if (!strncmp(string, "edge", 4) == 0) // check for correct character
      fatal("No 'edge' declaration");
  }

  fclose(input);

  if (sparseNet.getNumEdges() != numEdges)
    fatal("error recording edges in network");

  // set and save the initial state of network for further use
  sparseNet.setGlobalNetworkGE();

  std::cout << "\nFinding components and printing them to compX.gml files...\n"
            << std::endl;

  std::cout << numEdges << " edges explored" << std::endl;
  std::cout << dupEdges << " duplicate edges not counted in edge count" << std::endl;

  // check if all the nodes have invID lookups set
  if (DEBUG)
  {
    for (int i = min; i < max; ++i)
    {
      std::cout << "GETID(" << i << ") : " << sparseNet.getID(i) << std::endl;
    }
  }

  if ((output = fopen(output_fname, "a")) == NULL)
    fatal("File could not be opened.\n");
  fprintf(output, "%f\n", t.timeVal());
  fclose(output);

  exit(0);

  // GA crossover, mutation and other operators could be modified in the erbga.h header file before the run. defaults will picked from the header until unless passed into the command line.

  // Genetic Algorithm Starts here
  GA sparseGA(sparseNet, popSize, generations, numNodes, numEdges);
  sparseGA.set_data_name(input_fname);
  std::cout << "\nGenetic Algo Start : " << std::endl
            << "-----------------------" << std::endl;
  sparseGA.generate_GA();

  t.stop("Timer stopped");
  std::cout << t << " seconds" << std::endl;
  return 0;
}
