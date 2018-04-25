/***********************************************************
/
/   Header file for network data structure
/
/   Sharlee Climer 
/   June 2007
/
/   Breadth-first search added October/November 2007
/
/   Aditya Karnam
/   July 2017 - March 2018
/   Added Remove Edge
/   Added original degree to save initial degree states for Objective func
/   Added function to return the list of endpoints for a given vertices
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
#include "erbga_driver.h"

const int RETAINSYMMETRIC = 1;    // 1 to retain both (i,j) and (j,i) (always use 1 when using BFS)
const int DESCRIPTIVE_OUTPUT = 0; // 0 for just cluster membership numbers
const int PRINTCOMPSIZES = 0;     // 1 to print summary of component sizes out to "saveCompSizes.txt"
const int BFS_GML = 0;            // 1 to output BFS components to .gml files
const int BFS_WG2 = 0;            // 1 to output BFS components to .wg2 files
const double MIN_COMPLETE = 1.1;  // print comp. only if less than this value (0.5)
const int MAX_NUM_COMPS = 500000; // max number of non-singleton components
const double TOL = 0.000001;      // tolerance
const int NETWORK_API_DEBUG = 0;  // debug messages for API functions

class Edge
{
  friend class GA;
  friend class Vertex; // Vertex class allowed access to private data
  friend class Network;

private:
  Edge(int = 0, double = 1.0); // create a new edge with default weight of 1
  int getTarget();             // get other edgepoint
  double getWeight();          // get weight of edge
  int target;                  // other edgepoint
  double weight;               // weight
  Edge *next;                  // pointer to next edge eminating from source vertex
};

class Vertex
{
  friend class Network; // Network class allowed access to private functions
  friend class GA;

private:
  Vertex();            // create a new vertex
  ~Vertex();           // destructor
  int index;           // index number of node
  int degree;          // out-degree of vertex (full degree for undirected)
  int original_degree; // original degree of the node in the original network
  double degree_rate;  // degree rate defines the deviation of the degree of node from the max degree node
  Edge firstEdge;      // first edge eminating from vertex
                       // change to edgeListPlaceholder.next
                       //int ID; // GML ID number of vertex
                       //char *label; // GML label of vertex (NULL if none specified)
public:
  int addEdge(int, double = 1.0); // add an edge to vertex, return 1 if successful
  int haveEdge(int);              // check to see if edge already exists, return 0 if not
  double getWeight(int);          // return weight of edge
  void changeWeight(int, double); // change weight of edge
  int getDegree();                // get out-degree of vertex (full degree for undirected)
  void getEdges();                // get copy of edges eminating from vertex
  void printEdges(int, char *);   // print edges eminating from vertex
  int removeEdge(int end);        // remove edge from the edge list eminating from vertex
  //int getID(); // get GML ID number of vertex (-1 if none specified)
  //char *getLabel(); // get GML label of vertex (NULL if none specified)
  int getOriginalDegree();  // returns original_degree
  void getEndPoints(int &); // get the list of end points
};

class Network
{
  friend class GA;

public:
  Network(int = 0, int = 0, int = 0, int = 0); // create a network with 0 vertices, minimumIndex = 0 and maxIndex = 0
  ~Network();                                  // destructor
  int getNvertices();                          // get number of vertices in network
  int getNumEdges();                           // get number of edges in network
  int getDegree(int);                          // get degree of vertex
  int isDirected();                            // return 0 if undirected, 1 if directed network
  int addEdge(int, int, double = 1);           // add edge for class internal return 1 if edge successfully added
  int haveEdge(int, int);                      // return 1 if edge is in graph
  void printEdges(char *);                     // print all edges in network
  void bfs(char *);                            // breadth-first search outputs connected components
  double q_calc(char *);                       // clone of bfs with q_value calculation
  // void printAllEdges(char *); // print all the edges of the Network
  double getEdgeWeight(int, int);             // return weight of edge, given endpoints
  int removeEdge(int, int);                   // remove edge (x,y) from the list return bool
  void assignID(int, int);                    // creates the ID and invID for easy lookup for (vertex, actualVertexID)
  int getID(int);                             // return the ID for easy lookup
  int addEdgeFromGML(int, int, double = 1.0); // add edge from GML file return 1 if edge successfully added
  int *globalClusterNum;                      // stores the cluster numbering after each BFS run
  void setGlobalNetworkGE();                  // set the original vertives and edges values from read graph
  double modularity(char *);                  // calculates the modularity of the network
  int getOriginalDegree(int);                 // gets the original degree of vertex
  void getNodeEndPoints(int, int &);
  double average_degree;

private:
  int numVertices;         // number of vertices
  int directed;            // 0 if undirected, 1 if directed
  int numEdges;            // number of edges
  int originalNumVertices; // original graph nodes count
  int originalNumEdges;    // original graph edges count
  Vertex *vertices;        // vertices in network
  int *id;                 // hold the node IDs as given in input network
  int *invID;              // inverted ID numbers for easy look-up
  int *kcluster;           // number of clusters in the network
};

#endif
