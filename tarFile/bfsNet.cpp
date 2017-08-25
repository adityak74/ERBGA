/****************************************************************************
*
*	bfsNet.cpp:	Code for exploring a network using breadth-first 
*                       search.  Outputs connected components.
*
*                       Sharlee Climer
*                       October/November 2007
*
*      Modified to allow duplicate edges, count them separately, 
*      and record highest weight from the duplicates.
*      Sharlee Climer, July 2009
*
****************************************************************************/
  

#include "bfsNet.h"
#include "network.h"


int main(int argc, char ** argv)
{
  if (argc != 3)
    fatal("Usage:\n  bfs input.gml output.bfs"); 

  FILE *input;
  FILE *output;
  timer t;

  t.start("Timer started");

  int testFile = system("test ! -f comp1.gml");
  if(testFile != 0)
	fatal("Component files already exist in this directory");
  
  if (((input = fopen(argv[1], "r")) == NULL) || ((output = fopen(argv[2], "w"))
 == NULL))
    fatal("File could not be opened.\n");

  if(DESCRIPTIVE_OUTPUT)
    fprintf(output, "%s\n", argv[1]);

  int min = 1000000000;   // hold min node number
  int max = -1000000000;  // hold max node number
  char string[50];
  int startOne; // 1 if start node number is 1, 0 if start number is 0

  while (1) { // throw away everything until get to "graph" declaration
    if(feof(input)) fatal("no graph to read in input file");

    fscanf(input, "%s", string); // read in string from file

    if(strncmp(string, "graph", 5) == 0) 
      break; // stop when get to "graph" declaration
  }

  // find number of nodes in graph
  int numNodes = 0; // number of nodes
  int numEdges = 0; // number of edges

  while (1) { // throw away everything until get to "edge" declaration
    if(feof(input)) fatal("no edges to read in input file");

    fscanf(input, "%s", string); // read in string from file

    if(strncmp(string, "id", 2) == 0) { 
      fscanf(input, "%s", string); // read in id value
      int num = atoi(string);

      if(min > num) min = num; // update minimum node number
      if(max < num) max = num; // update maximum node number
      numNodes++; // one more node found
    }

    if(strncmp(string, "edge", 4) == 0) 
      break; // stop when read in first "edge"
  }

  std::cout << "\nNode numbers range from " << min <<" to " << max << std::endl
;
  std::cout << "Reading in graph from " << argv[1] << "...\n" << std::endl;

  if (DESCRIPTIVE_OUTPUT)
    fprintf(output, "Node numbers range from %d to %d\n\n", min, max);
  fclose(output); 

  startOne = min;
  rewind(input); // start at beginning of input file again

  if(startOne < 0) fatal("Node numbers can not be negative");
  if ((startOne != 0) && (startOne != 1)) 
    warning("Smallest node number should be 0 or 1");

  // record node id numbers
  int ptr = 0; // pointer for filling nodeNumbers array
  int *id; // hold the node IDs as given in input file
  if ((id = new int[numNodes]) == NULL)
    fatal("memory not allocated");

  while (1) { // throw away everything until get to "graph" declaration
    if(feof(input)) fatal("no graph to read in input file");

    fscanf(input, "%s", string); // read in string from file

    if(strncmp(string, "graph", 5) == 0) 
      break; // stop when get to "graph" declaration
  }

  while (1) { // record node id numbers 
    if(feof(input)) fatal("no edges to read in input file");

    fscanf(input, "%s", string); // read in string from file
  
    if(strncmp(string, "id", 2) == 0) { 
      fscanf(input, "%s", string); // read in id value
      int num = atoi(string);

      if((num < min) || (num > max))
	fatal("error while reading in node numbers");

      id[ptr++] = num; // record node id number
      //std::cout << ptr << ", " << id[ptr-1]<< std::endl;
    }

    if(strncmp(string, "edge", 4) == 0) 
      break; // stop when read in first "edge"
  }

  if (ptr != numNodes)
    fatal("Error reading in node numbers");

  int *idInv; // invert the ID numbers for easy look-up
  if ((idInv = new int[max+1]) == NULL)
    fatal("memory not allocated");

  //std::cout << "id: " ;
  for (int i = 0; i < max+1; i++)
    idInv[i] = -1; // initialize values
  for (int i = 0; i < numNodes; i++) {
    //std::cout << id[i] << " ";
    if(idInv[id[i]] != -1)
      fatal("Error recording index for ID number");
    idInv[id[i]] = i; // record index for given id number
  }
  //std::cout << std::endl;

  if (0) {
    std::cout << "idInv: ";
    for (int i = 0; i < max+1; i++)
      std::cout << idInv[i] << " ";
    std::cout << std::endl;
  }

  // create network with numNodes vertices
  Network sparseNet(numNodes, DIRECTED); // if DIRECTED = 0, undirected
  int dupEdges = 0; // record number of duplicate edges

  while (!feof(input)) { // read edges until reach end of file
    fscanf(input, "%s", string); // read in string from file

    if(!strncmp(string, "[", 1) == 0) // check for correct character
      fatal("No '[' after edge declaration");

    fscanf(input, "%s", string); // read in string from file
    if(!strncmp(string, "source", 6) == 0) // check for correct character
      fatal("No 'source' declaration");

    fscanf(input, "%s", string); // read in start node from file
    int source = atoi(string);

    if((source < min) || (source > max))
      fatal("Invalid node number");

    //std::cout << "source: " << idInv[source];

    fscanf(input, "%s", string); // read in string from file
    if(!strncmp(string, "target", 6) == 0) // check for correct character
      fatal("No 'target' declaration");

    fscanf(input, "%s", string); // read in target node from file
    int target = atoi(string);

    if((target < min) || (target > max))
      fatal("Invalid node number");

    //std::cout << ", target: " << idInv[target];

    fscanf(input, "%s", string); // read in string from file

    float weight = 1.0;
    if(strncmp(string, "weight", 6) == 0) {// check to see if weight specified
      fscanf(input, "%s", string); // read in string from file
      weight = atof(string);
      fscanf(input, "%s", string); // read in string from file
    }

    //std::cout << ", weight: " << weight << std::endl;
    
    if(!sparseNet.addEdge(idInv[source],idInv[target],weight))
	  dupEdges++;
      //warning("Duplicate edge in input");
    
    else
      numEdges++; // count number of edges

	if(numEdges % 10000000 == 0) // message every 10 million edges
	  std::cout << numEdges / 1000000 << " million edges read" << std::endl;

    if(!strncmp(string, "]", 1) == 0) // check for correct character
      fatal("No ']' after edge declaration");

    //fprintf(output, "%d\t%d\t%f\n", source, target, weight);

    if (feof(input)) break; // break is at end of file

    fscanf(input, "%s", string); // read in string from file

    if(strncmp(string, "]", 1) == 0) // check for end of graph
      break;

    if(!strncmp(string, "edge", 4) == 0) // check for correct character
      fatal("No 'edge' declaration");
  }

  //std::cout << std::endl;
 
  fclose(input);

  if (sparseNet.getNumEdges() != numEdges)
    fatal("error recording edges in network");

  std::cout << "\nFinding components and printing them to compX.gml files...\n" << std::endl;

  sparseNet.bfs(argv[2]);

  std::cout << numEdges << " edges explored" << std::endl;
  std::cout << dupEdges << " duplicate edges not counted in edge count" << std::endl;
  //if (dupEdges > 0)
  //std::cout << "Highest edge weight recorded for duplicate edges" << std::endl;

  t.stop("Timer stopped");
  std::cout << t << " seconds" << std::endl;

  if ((output = fopen(argv[2], "a")) == NULL)
    fatal("File could not be opened.\n");
  fprintf(output,"%f\n",t.timeVal());
  fclose(output);

  delete [] id;
  delete [] idInv;

  return 1;
}

