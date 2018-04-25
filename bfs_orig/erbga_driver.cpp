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

int main(int argc, char **argv)
{
  // Randomizer for GA
  srand(time(NULL));
  if (argc != 3)
    fatal("Usage:\n   erbga input.gml output.bfs");

  opterr = 0;

  FILE *input;
  FILE *output;
  timer t;

  t.start("Timer started");

  int testFile = system("test ! -f comp1.gml");
  if (testFile != 0)
    fatal("Component files already exist in this directory");

  if (((input = fopen(argv[1], "r")) == NULL) || ((output = fopen(argv[2], "w")) == NULL))
    fatal("File could not be opened.\n");

  if (DESCRIPTIVE_OUTPUT)
    fprintf(output, "%s\n", argv[1]);

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
  std::cout << "Reading in graph from " << argv[1] << "...\n"
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

  if ((output = fopen(argv[2], "a")) == NULL)
    fatal("File could not be opened.\n");
  fprintf(output, "%f\n", t.timeVal());
  fclose(output);

  // GA params
  int popSize = 100;
  int generations = 500;

  // GA crossover, mutation and other operators could be modified in the erbga.h header file before the run. defaults will picked from the header until unless passed into the command line.

  // Genetic Algorithm Starts here
  GA sparseGA(sparseNet, popSize, generations, numNodes, numEdges);
  sparseGA.set_data_name(argv[1]);
  std::cout << "\nGenetic Algo Start : " << std::endl
            << "-----------------------" << std::endl;
  sparseGA.generate_GA();

  t.stop("Timer stopped");
  std::cout << t << " seconds" << std::endl;
  return 0;
}