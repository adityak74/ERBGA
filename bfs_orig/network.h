/***********************************************************
/
/   Header file for network data structure
/
/   Sharlee Climer 
/   June 2007
/
/   Breadth-first search added October/November 2007
/
/**********************************************************/

#ifndef NETWORK_H
#define NETWORK_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string.h>
#include <time.h>
#include "bfsNet.h"

const int RETAINSYMMETRIC = 1;// 1 to retain both (i,j) and (j,i) (use for BFS)
const int DESCRIPTIVE_OUTPUT = 0; // 0 for just cluster membership numbers
const int PRINTCOMPSIZES = 1; // 1 to print summary of component sizes out to "saveCompSizes.txt"
const int BFS_GML = 0; // 1 to output BFS components to .gml files
const int BFS_WG2 = 0; // 1 to output BFS components to .wg2 files
const double MIN_COMPLETE = 1.1; // print comp. only if less than this value (0.5)
const int MAX_NUM_COMPS = 500000; // max number of non-singleton components
const double TOL = 0.000001; // tolerance

class Edge
{
  friend class Vertex; // Vertex class allowed access to private data
  friend class Network; 
private: 
  Edge(int = 0, double = 1); // create a new edge with default weight of 1
  int getTarget(); // get other edgepoint  
  double getWeight(); // get weight of edge 
  int target; // other edgepoint
  double weight; // weight 
  Edge *next; // pointer to next edge eminating from source vertex
};

class Vertex
{
  friend class Network; // Network class allowed access to private functions
 private:
  Vertex(); // create a new vertex
  ~Vertex(); // destructor
  int index; // index number of node
  int degree; // out-degree of vertex (full degree for undirected)
  Edge firstEdge; // first edge eminating from vertex
  //int ID; // GML ID number of vertex
  //char *label; // GML label of vertex (NULL if none specified)

  int addEdge(int, double=1); // add an edge to vertex, return 1 if successful
  int haveEdge(int); // check to see if edge already exists, return 0 if not
  double getWeight(int); // return weight of edge
  void changeWeight(int, double); // change weight of edge
  int getDegree(); // get out-degree of vertex (full degree for undirected)
  Edge* getEdges(); // get copy of edges eminating from vertex
  void printEdges(int, char *); // print edges eminating from vertex
  //int getID(); // get GML ID number of vertex (-1 if none specified)
  //char *getLabel(); // get GML label of vertex (NULL if none specified)
};

class Network
{
 public:
  Network(int = 0, int = 0); // create a network with 0 vertices
  ~Network(); // destructor
  int getNvertices(); // get number of vertices in network
  int getNumEdges(); // get number of edges in network
  int getDegree(int); // get degree of vertex
  int isDirected(); // return 0 if undirected, 1 if directed network
  int addEdge(int, int, double=1); // return 1 if edge successfully added
  int haveEdge(int, int); // return 1 if edge is in graph
  void printEdges(char *); // print all edges in network
  void bfs(char *); // breadth-first search outputs connected components
  void q_calc(char *); // clone of bfs with q_value calculation
  // void printAllEdges(char *); // print all the edges of the Network
  double getEdgeWeight(int, int); // return weight of edge, given endpoints
  void removeEdge(int, int); // remove edge (x,y) from the list

 private:
  int numVertices; // number of vertices
  int directed; // 0 if undirected, 1 if directed
  int numEdges; // number of edges
  Vertex *vertices; // vertices in network
};

#endif
