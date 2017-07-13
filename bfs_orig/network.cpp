/***********************************************************
/
/   Implementation of network data structure
/
/   Sharlee Climer 
/   June 2007
/
/   Breadth-first search added October/November 2007
/   Print out component size summary (PRINTCOMPSIZES)
/    added January 2011
/
/**********************************************************/



#include "network.h"

using namespace std;


Edge::Edge(int endpt, double wt) // create a new edge with default weight of 1
{  
  if (weight < 0-TOL) fatal("Negative edge weight");
  if (endpt < 0)  fatal("Invalid endpoint for edge");

  target = endpt; // assign endpoint
  weight = wt; // assign weight
}

int Edge::getTarget() // get other edgepoint  
{
  return target;
}
  
double Edge::getWeight() // get weight of edge 
{
  return weight;  
}

Vertex::Vertex() // create a new vertex
{
  //if (indxNum < 0)  fatal("Vertex indices must be non-negative");
  //index = indxNum; // assign index of node
  //ID = GML; // assign GML ID number
  //label = lbl; // assign GML label
  degree = 0; // no incident edges yet
  firstEdge.next = 0; // first edge is pointed to by firstEdge.next
}

Vertex::~Vertex() // destructor
{
  /* delete edges; 
  Edge *edgePtr, *followPtr; // pointers to move through linked list of edges
  edgePtr = &firstEdge; // point to first edge
  followPtr = edgePtr;  

  while (edgePtr->next != 0) {// follow until find last edge
    edgePtr = edgePtr->next; // pointer points at next edge in list
    delete followPtr;
    followPtr = edgePtr;
  }
  Need to fix!!!!!!!!!!!!!!!!!!!!!!!!
  */
}

//int Vertex::getID() // get GML ID number of vertex

//note: firstEdge is a placeholder that points to true first edge
int Vertex::addEdge(int endpt, double wt) // add an edge to vertex, return 1 if successful
{
  if (endpt < 0) fatal("invalid endpoint for edge");
  if (wt < 0-TOL) fatal("invalid weight for edge");

  Edge *edgePtr; // pointer to move through linked list of edges
  edgePtr = &firstEdge; // point to first edge

  while (edgePtr->next != 0) // follow until find last edge
    edgePtr = edgePtr->next; // pointer points at next edge in list

  Edge *newEdge; // newEdge created
  if((newEdge = new Edge) == NULL)
    fatal("memory not allocated for edge");

  edgePtr->next = newEdge; // add new edge to list
  newEdge->target = endpt; // assign properties to new edge
  newEdge->weight = wt;
  newEdge->next = 0; // mark new edge as last on list

  return 1;
}

double Vertex::getWeight(int endpoint) // return weight of edge
{
  Edge *edgePtr; // pointer to move through linked list of edges
  edgePtr = &firstEdge; // point to first edge

  for (int i = 0; i < endpoint; i++)
    edgePtr = edgePtr->next; // pointer points at next edge in list

  return edgePtr->weight;
}

void Vertex::changeWeight(int endpoint, double newWeight) // change weight of edge
{
  Edge *edgePtr; // pointer to move through linked list of edges
  edgePtr = &firstEdge; // point to first edge

  for (int i = 0; i < endpoint; i++)
	edgePtr = edgePtr->next; // pointer points at next edge in list

  edgePtr->weight = newWeight;
}

int Vertex::haveEdge(int endpoint) // return 1 if edge exists, 0 if not
{
  Edge *edgePtr; // pointer to move through linked list of edges
  edgePtr = &firstEdge; // point to first edge

  if (edgePtr->next == 0) // no edges for this vertex
    return 0;

  while (edgePtr->next != 0) { // follow until find last edge
    edgePtr = edgePtr->next; // pointer points at next edge in list
    if (edgePtr->target == endpoint)
      return 1;
  }

  return 0;
}

int Vertex::getDegree() // get out-degree (full degree for undirected)
{
  return degree;
}

//char* Vertex::getLabel() // get GML label of vertex (NULL if none specified)

Edge* Vertex::getEdges() // get copy of edges eminating from vertex
{
  Edge* edge1; // FIX
  return edge1;
}

void Vertex::printEdges(int node, char *outputFile) // print edges eminating from vertex
{
  Edge *edgePtr; // pointer to move through linked list of edges
  edgePtr = &firstEdge; // point to first edge

  if (edgePtr->next == 0) return; // no edges from this node

  FILE *output; 
  if ((output = fopen(outputFile, "a")) == NULL)
      fatal ("Unable to open output file");

  while (edgePtr->next != 0) { // follow until find last edge
    edgePtr = edgePtr->next; // pointer points at next edge in list

    if(!RETAINSYMMETRIC || (node < edgePtr->target)) // print only once
      fprintf(output,"\n\tedge\n\t[\n\tsource %d\n\ttarget %d\n\t]", node+1, edgePtr->target+1);
  }

  fclose(output);
}

Network::Network(int num, int dir) // create a network with num vertices
{
  if (num < 1)  fatal("Network requires at least 1 vertex");
  numVertices = num;
  directed = dir; // 0 if undirected, 1 if directed network
  if((vertices = new Vertex[numVertices]) == NULL)
    fatal("memory not allocated"); // allocate memory

  numEdges = 0; // no edges added yet
}


Network::~Network() // destructor
{
  delete [] vertices; 
}

int Network::getNvertices() // get number of vertices in network
{
  return numVertices;
}

int Network::getNumEdges() // get number of vertices in network
{
  return numEdges;
}

int Network::getDegree(int vertex) // get degree of node
{
  int degree = vertices[vertex].getDegree();
  return degree;
}

int Network::isDirected() // return 0 if undirected, 1 if directed network
{
  return directed;
}

int Network::addEdge(int v1, int v2, double weight) // return 1 if successfully add edge
{
  if ((v1 > numVertices-1) || (v2 > numVertices-1))
    fatal("Attempt to add edge to non-existent node");
  if ((v1 < 0) || (v2 < 0))
    fatal("Attempt to add edge to negative numbered node");

  if (!DIRECTED) {
    if (v1 > v2) { // order v1 and v2
      int temp = v1;
      v1 = v2;
      v2 = temp;
    }
  }

  //std::cout << "v1 = " << v1 << ", v2 = " << v2 << std::endl;

  if(vertices[v1].haveEdge(v2)) { // edge already exist
	//	if (weight > vertices[v1].getWeight(v2))
	//  vertices[v1].changeWeight(v2, weight); // use highest weight for edge
    return 0; // signal no new edge added
  }

  vertices[v1].addEdge(v2, weight); // add edge

  if(RETAINSYMMETRIC) { // retain both (i,j) and (j,i)
    if(vertices[v2].haveEdge(v1)) // edge already exists, return failure
      return 0;
    vertices[v2].addEdge(v1, weight);
  }

  numEdges++; // update number of edges
  vertices[v1].degree++; // update degree of vertices
  if (!DIRECTED)
    vertices[v2].degree++; // increase degree for both vertices if undirected

  return 1;
}

void Network::printEdges(char *outputFile) // print all edges in network
{
  FILE *output;
  if ((output = fopen(outputFile, "w")) == NULL)
      fatal ("Unable to open output file");

  fprintf(output,"Graph with %d nodes and %d edges.\ngraph\n[", numVertices, numEdges);  

  for (int i = 0; i < numVertices; i++)
    fprintf(output,"\n  node\n  [\n    id %d\n  ]", i+1);

  fclose(output);

  for (int i = 0; i < numVertices; i++)
    vertices[i].printEdges(i, outputFile);

  if ((output = fopen(outputFile, "a")) == NULL)
      fatal ("Unable to open output file");

  fprintf(output,"\n]\n");  
  fclose(output);
}

double Network::getEdgeWeight(int v1, int v2) // get weight of edge between endpoints
{
  if ((v1 > numVertices-1) || (v2 > numVertices-1)) 
    fatal("Attempt to get weight of edge to non-existent node");
  if ((v1 < 0) || (v2 < 0))
    fatal("Attempt to get weight of edge to negative numbered node");
  
  if (!DIRECTED) {
    if (v1 > v2) { // order v1 and v2
      int temp = v1;
      v1 = v2;
      v2 = temp;
    }
  }

  if(!vertices[v1].haveEdge(v2)) // edge is not in graph
	fatal("Attempt to get weight of edge that is not in graph");

  double wt = vertices[v1].getWeight(v2); // get weight
  if (wt < 0-TOL)
	fatal ("Negative edge weight");

  return wt;
}

int Network::haveEdge(int v1, int v2) // return 1 if edge is in graph
{
  if ((v1 > numVertices-1) || (v2 > numVertices-1)) 
    fatal("Attempt to check edge to non-existent node");
  if ((v1 < 0) || (v2 < 0))
    fatal("Attempt to check edge to negative numbered node");
  
  if (!DIRECTED) {
    if (v1 > v2) { // order v1 and v2
      int temp = v1;
      v1 = v2;
      v2 = temp;
    }
  }

  if(vertices[v1].haveEdge(v2)) // edge is in graph
    return 1;
  else 
    return 0;
}

//Removed BFS function