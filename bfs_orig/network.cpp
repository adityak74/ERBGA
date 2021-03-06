// *********************************************************
//
//    Implementation of network data structure
//
//    Sharlee Climer 
//    June 2007
//
//    Breadth-first search added October/November 2007
//    Print out component size summary (PRINTCOMPSIZES)
//     added January 2011
//
//
//    Modified by Aditya Karnam
//
//    Add New Edge to the front of the list
//    Moved ID and invID from BFS to Network
//    assignID and getID functions added
//    Added remove edge
//    Fixed garbage collection
//    Added Qs    
//    Added Modularity
// 
//    July 2017 - March 2018
//
// *********************************************************

#include "network.h"

using namespace std;

// Edge Constructor
// create a new edge with default weight of 1
Edge::Edge(int endpt, double wt)
{
    if (wt < 0 - TOL) // fixed was weight instead of wt
        fatal("Negative edge weight");

    if (endpt < 0)
        fatal("Invalid endpoint for edge");

    target = endpt; // assign endpoint
    weight = wt; // assign weight
}

// getTarget() returns int
// get other edgepoint
int Edge::getTarget()
{
    return target;
}

// getWeight() returns double
// get weight of edge
double Edge::getWeight()
{
    return weight;
}

// Vertex Constructor
// create a new vertex
Vertex::Vertex()
{
    //if (indxNum < 0)  fatal("Vertex indices must be non-negative");
    //index = indxNum; // assign index of node
    //ID = GML; // assign GML ID number
    //label = lbl; // assign GML label
    degree = 0; // no incident edges yet
    original_degree = 0; // placeholder for holding the degree of the node in original network
    firstEdge.next = NULL; // first edge is pointed to by firstEdge.next
}

// Vertex destructor
Vertex::~Vertex()
{
    Edge *edgePtr, *followPtr; // pointers to move through linked list of edges
    edgePtr = firstEdge.next; // point to first edge (head of the list)
    followPtr = edgePtr;

    while (edgePtr != NULL) { // follow until find last edge
        edgePtr = edgePtr->next; // pointer points at next edge in list
        delete followPtr;
        followPtr = edgePtr;
    }
}

//int Vertex::getID() // get GML ID number of vertex

// Vertex addEdge() returns int
// note: firstEdge is a placeholder that points to true first edge
// add an edge to vertex, return 1 if successful
int Vertex::addEdge(int endpt, double wt)
{
    if (endpt < 0)
        fatal("invalid endpoint for edge");
    if (wt < 0 - TOL)
        fatal("invalid weight for edge");

    Edge* edgePtr; // pointer to move through linked list of edges
    edgePtr = &firstEdge; // point to first edge

    Edge* newEdge; // newEdge created
    if ((newEdge = new Edge) == NULL)
        fatal("memory not allocated for edge");

    newEdge->target = endpt; // assign properties to new edge
    newEdge->weight = wt;
    newEdge->next = edgePtr->next; // add new edge in the beginning of the list
    firstEdge.next = newEdge; // Update the firstEdge to point to the newEdge

    return 1;
}

// Vertex::getWeight() returns double
// return weight of edge
double Vertex::getWeight(int endpoint)
{
    Edge* edgePtr; // pointer to move through linked list of edges
    edgePtr = &firstEdge; // point to first edge

    for (int i = 0; i < endpoint; i++)
        edgePtr = edgePtr->next; // pointer points at next edge in list

    return edgePtr->weight;
}

// Vertex::changeWeight()
// update weight of edge
void Vertex::changeWeight(int endpoint, double newWeight)
{
    Edge* edgePtr; // pointer to move through linked list of edges
    edgePtr = &firstEdge; // point to first edge

    for (int i = 0; i < endpoint; i++)
        edgePtr = edgePtr->next; // pointer points at next edge in list

    edgePtr->weight = newWeight;
}

// Vertex::haveEdge() returns int
// return 1 if edge exists, 0 if not
int Vertex::haveEdge(int endpoint)
{
    Edge* edgePtr; // pointer to move through linked list of edges
    edgePtr = &firstEdge; // point to first edge

    if (edgePtr->next == NULL) // no edges for this vertex
        return 0;

    while (edgePtr->next != NULL) { // follow until find last edge
        edgePtr = edgePtr->next; // pointer points at next edge in list
        if (edgePtr->target == endpoint)
            return 1;
    }

    return 0;
}

// Vertex::getDegree() returns int
// get out-degree (full degree for undirected)
int Vertex::getDegree()
{
    return degree;
}

// Vertex::getOriginalDegree() returns int
// get out-degree (full degree for undirected)
int Vertex::getOriginalDegree()
{
    return original_degree;
}

//char* Vertex::getLabel() // get GML label of vertex (NULL if none specified)

// Vertex::getEdges() returns pointer to Edge list
// get copy of edges eminating from vertex
void Vertex::getEdges()
{
    Edge* edgePtr; // pointer to move through linked list of edges
    edgePtr = &firstEdge; // point to first edge

    if (edgePtr->next == NULL)
        return; // no edges from this node

    while (edgePtr->next != NULL) { // follow until find last edge
        edgePtr = edgePtr->next; // pointer points at next edge in list

        // if(!RETAINSYMMETRIC || (node < edgePtr->target)) // print only once
    }
}

// Vertex::printEdges
// print edges eminating from vertex
void Vertex::printEdges(int node, char* outputFile)
{
    Edge* edgePtr; // pointer to move through linked list of edges
    edgePtr = &firstEdge; // point to first edge

    if (edgePtr->next == NULL)
        return; // no edges from this node

    FILE* output;
    if ((output = fopen(outputFile, "a")) == NULL)
        fatal("Unable to open output file");

    while (edgePtr->next != NULL) { // follow until find last edge
        edgePtr = edgePtr->next; // pointer points at next edge in list

        if (!RETAINSYMMETRIC || (node < edgePtr->target)) // print only once
            fprintf(output, "\n\tedge\n\t[\n\tsource %d\n\ttarget %d\n\t]", node + 1, edgePtr->target + 1);
    }

    fclose(output);
}

// Network Constructor
// create a network with num vertices, DIRECTED/INDIRECT, min and max nodeID
Network::Network(int numNodes, int dir, int min, int max)
{
    if (numNodes < 1)
        fatal("Network requires at least 1 vertex");
    numVertices = numNodes;
    directed = dir; // 0 if undirected, 1 if directed network
    if ((vertices = new Vertex[numVertices]) == NULL)
        fatal("memory not allocated"); // allocate memory

    numEdges = 0; // no edges added yet

    if ((id = new int[numNodes]) == NULL) // allocate memory to node IDs
        fatal("memory not allocated");
    if ((invID = new int[max + 1]) == NULL) // allote memory to invIDs
        fatal("memory not allocated");
    for (int i = 0; i < max + 1; i++)
        invID[i] = -1; // initialize values for invID some values can be -1
    if ((kcluster = new int) == NULL) // allocate memory to node IDs
        fatal("memory not allocated");
    // allocate memory for edgeIDs change to numEdges
}

// Network::assignID()
// assigns id and invID to each node
void Network::assignID(int index, int nodeIdNum)
{
    id[index] = nodeIdNum; // record node id number
    invID[id[index]] = index; // record index for given id number
}

// Network::getID() returns int
// returns the lookup index for the vertex
int Network::getID(int vertex)
{
    return invID[vertex];
}

// ~Network() destructor
Network::~Network()
{
    delete[] id;
    delete[] invID;
    delete[] vertices;
    delete[] globalClusterNum;
}

// Network::getNvertices() returns int
// get number of vertices in network
int Network::getNvertices()
{
    return numVertices;
}

// Network::getNumEdges() returns int
// get number of vertices in network
int Network::getNumEdges()
{
    return numEdges;
}

// Network::getDegree() returns int
// get degree of node
int Network::getDegree(int vertex)
{
    int degree = vertices[vertex].getDegree();
    return degree;
}

// Network::getOriginalDegree() returns int
// get degree of node
int Network::getOriginalDegree(int vertex)
{
    int degree = vertices[vertex].getOriginalDegree();
    return degree;
}

// Network::isDirected()
// return 0 if undirected, 1 if directed network
int Network::isDirected()
{
    return directed;
}

// Network::addEdge() returns int
// return 1 if successfully add edge from v1 to v2 with weight
int Network::addEdgeFromGML(int v1, int v2, double weight)
{
    // Lookup the actual ID for v1,v2
    v1 = invID[v1];
    v2 = invID[v2];

    if ((v1 > numVertices - 1) || (v2 > numVertices - 1))
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

    if (vertices[v1].haveEdge(v2)) { // edge already exist
        //	if (weight > vertices[v1].getWeight(v2))
        //  vertices[v1].changeWeight(v2, weight); // use highest weight for edge
        return 0; // signal no new edge added
    }

    vertices[v1].addEdge(v2, weight); // add edge

    if (RETAINSYMMETRIC) { // retain both (i,j) and (j,i)
        if (vertices[v2].haveEdge(v1)) // edge already exists, return failure
            return 0;
        vertices[v2].addEdge(v1, weight);
    }

    numEdges++; // update number of edges
    vertices[v1].degree++; // update degree of vertices
    if (!DIRECTED)
        vertices[v2].degree++; // increase degree for both vertices if undirected

    return 1;
}

// Network::addEdge() returns int
// return 1 if successfully add edge from v1 to v2 with weight
int Network::addEdge(int v1, int v2, double weight)
{
    // Lookup the actual ID for v1,v2
    // v1 = invID[v1];
    // v2 = invID[v2];

    if ((v1 > numVertices - 1) || (v2 > numVertices - 1))
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

    if (vertices[v1].haveEdge(v2)) { // edge already exist
        //  if (weight > vertices[v1].getWeight(v2))
        //  vertices[v1].changeWeight(v2, weight); // use highest weight for edge
        return 0; // signal no new edge added
    }

    vertices[v1].addEdge(v2, weight); // add edge

    if (RETAINSYMMETRIC) { // retain both (i,j) and (j,i)
        if (vertices[v2].haveEdge(v1)) // edge already exists, return failure
            return 0;
        vertices[v2].addEdge(v1, weight);
    }

    numEdges++; // update number of edges
    vertices[v1].degree++; // update degree of vertices
    if (!DIRECTED)
        vertices[v2].degree++; // increase degree for both vertices if undirected

    return 1;
}

int Vertex::removeEdge(int end)
{
    Edge *edgePtr, *prevEdgePtr;
    int edgeFound = 0;

    // since this is initiliazed in stack no pointer reference.
    edgePtr = &firstEdge; // fix for first and last element deletion and decerement the degrees
    prevEdgePtr = edgePtr; // prevEdge points to the first Edge for initial setup

    if (edgePtr->next == NULL) // no edges for this vertex
        return 0;

    while (edgePtr->next != NULL) { // follow until find last edge
        prevEdgePtr = edgePtr;
        edgePtr = edgePtr->next; // pointer points at next edge in list

        // std::cout << "DEBUG | at here :->" << edgePtr->target << " |" <<std::endl;
        if (edgePtr->target == end) {
            edgeFound = 1;
            break;
        }
    }
    if (edgeFound) {
        // std::cout << "DEBUG | edge found :->" << edgePtr->target << " |" <<std::endl;
        prevEdgePtr->next = edgePtr->next;
        if (edgePtr != NULL)
            delete edgePtr;
        return 1;
    }
    else {
        fatal("No such Edge in the List");
        return 0;
    }
}

// Network::removeEdge() returns int
// removes edge from v1 to v2, if success returns 1 else 0
int Network::removeEdge(int v1, int v2)
{

    int start = (v1 < v2) ? v1 : v2;
    int end = (v1 < v2) ? v2 : v1;

    if ((start > numVertices - 1) || (end > numVertices - 1))
        fatal("Attempt to add edge to non-existent node");
    if ((start < 0) || (end < 0))
        fatal("Attempt to add edge to negative numbered node");

    if (!vertices[v1].haveEdge(v2)) {
        return 0; // no such edge in the list
    }

    int ret = vertices[start].removeEdge(end);
    if (RETAINSYMMETRIC) {
        if (vertices[end].haveEdge(start)) {
            // std::cout << "--- Removing Symmetric Edge : (" << end << "," << start << ")" << std::endl;
            vertices[end].removeEdge(start);
        }
    }

    if (ret) {
        vertices[start].degree--; // update degree of vertices
        numEdges--; // update number of edges
        if (!DIRECTED)
            vertices[end].degree--; // increase degree for both vertices if undirected
        return 1;
    }
}

// Network::printEdges()
// prints the network to a file in GML format
void Network::printEdges(char* outputFile) // print all edges in network
{
    FILE* output;
    if ((output = fopen(outputFile, "w")) == NULL)
        fatal("Unable to open output file");

    fprintf(output, "Graph with %d nodes and %d edges.\ngraph\n[", numVertices, numEdges);

    for (int i = 0; i < numVertices; i++)
        fprintf(output, "\n  node\n  [\n    id %d\n  ]", i + 1);

    fclose(output);

    for (int i = 0; i < numVertices; i++) {
        vertices[i].printEdges(i, outputFile);
        vertices[i].getEdges();
    }

    if ((output = fopen(outputFile, "a")) == NULL)
        fatal("Unable to open output file");

    fprintf(output, "\n]\n");
    fclose(output);
}

void Network::setGlobalNetworkGE() {
    originalNumEdges = numEdges;
    originalNumVertices = numVertices;

    // assign the degree of the original_degree of every node
    // this will be used for calculating modularity

    int max_degree = 0, total_degree = 0;
    double average_degree = 0.0f;

    for (int i = 0; i < numVertices; i++)
    {
        vertices[i].original_degree = vertices[i].degree;
        total_degree += vertices[i].degree;
        if ( max_degree < vertices[i].degree )
            max_degree = vertices[i].degree;

        if ( NETWORK_API_DEBUG )
            std :: cout << "====Vertex #" << i+1 << ", degree : " << vertices[i].degree << "\n";
    }

    // set the average degree
    average_degree = (double) total_degree / numVertices;

    if ( NETWORK_API_DEBUG ) {
        std :: cout << "====max degree : " << max_degree << "\n";
        std :: cout << "====average degree : " << average_degree << "\n";
    }

    for (int i = 0; i < numVertices; i++)
    {
        vertices[i].degree_rate = (double)(max_degree - vertices[i].degree) / max_degree + (double)(abs(average_degree - vertices[i].degree)) / max_degree;
        
        if ( !NETWORK_API_DEBUG )
            std :: cout << "====Vertex #" << i+1 << ", degree_rate : " << vertices[i].degree_rate << " max degree : " << max_degree<< "\n";
    }
    
}

    // Network::getEdgeWeight() returns double
    // get weight of edge between endpoints
    double
    Network::getEdgeWeight(int v1, int v2)
{
    if ((v1 > numVertices - 1) || (v2 > numVertices - 1))
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

    if (!vertices[v1].haveEdge(v2)) // edge is not in graph
        fatal("Attempt to get weight of edge that is not in graph");

    double wt = vertices[v1].getWeight(v2); // get weight
    if (wt < 0 - TOL)
        fatal("Negative edge weight");

    return wt;
}

void Vertex::getEndPoints(int& endPointsList) {
    Edge *edgePtr;        // pointer to move through linked list of edges
    edgePtr = &firstEdge; // point to first edge

    if (edgePtr->next == NULL)
        return; // no edges from this node
    int index = 0;
    int* p = &endPointsList;
    while (edgePtr->next != NULL)
    {                            // follow until find last edge
        edgePtr = edgePtr->next; // pointer points at next edge in list
        *(p+(index++)) = edgePtr->target;
    }
    if (NETWORK_API_DEBUG) {
        std :: cout << "EPTS index : " << index << std :: endl;
    }
}

// returns the endPointsList reference
void Network::getNodeEndPoints(int node, int& endPointList) {
    if ((node > numVertices - 1))
        fatal("Attempt to check edge to non-existent node");
    if ((node < 0))
        fatal("Attempt to check edge to negative numbered node");

    vertices[node].getEndPoints(endPointList);
}

// Network::haveEdge() returns int
// return 1 if edge is in graph
int Network::haveEdge(int v1, int v2)
{
    if ((v1 > numVertices - 1) || (v2 > numVertices - 1))
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

    if (vertices[v1].haveEdge(v2)) // edge is in graph
        return 1;
    else
        return 0;
}

// breadth-first search, outputs connected components
void Network::bfs(char* outputFile)
{
    FILE* output;

    int testFile = system("test ! -f comp1.gml");
    if (testFile != 0)
        fatal("Component files already exist in this directory");

    if ((output = fopen(outputFile, "w")) == NULL)
        fatal("Unable to open output file");

    int ptr; // point to current index of component array
    int startPtr, endPtr; // pointers for queue
    int *visited, *component, *queue, *clusterNum;
    int numSingle = 0; // count nodes with degree zero
    int numTwo = 0; // count number of components with 2 vertices
    int k = 0; // current cluster number
    int totalEdges = 0; // total edges in graph
    int numComp3orMore = 0; // number of components with 3 or more vertices
    FILE* outputGML; // output components to .gml files
    FILE* outputWG2; // output components to .wg2 files
    FILE* outputCompnn; // output node numbering for components
    int randomVal; // random value to place at end of comp file names
    int numNonSingle = 0; // number of non-singleton components
    int printNum = 0; // number of cluster printed out
    int maxSizeComp = 0; // max size of component printed to gml file
    int maxCompAll = 0; // max size of component overall
    int num[10000]; // hold number of components with each size
    int numNotCliques = 0; // hold number of components that aren't cliques

    // globalCusterNum - accessed across Network class

    for (int i = 0; i < 10000; i++)
        num[i] = 0;

    int compSize[MAX_NUM_COMPS]; // hold non-singleton component sizes
    double compDensity[MAX_NUM_COMPS]; // hold non-singleton component densities
    int compPtr = 0; // pointer for these two arrays

    if (BFS_WG2 && !RETAINSYMMETRIC)
        fatal("Need to retain symmetric edges for .wg2 output");

    // srand(time(NULL)); // assign random value between 1 - 1000
    randomVal = rand() % 1000 + 1;

    if (((visited = new int[numVertices]) == NULL) || ((component = new int[numVertices]) == NULL))
        fatal("memory not allocated"); // allocate memory
    if (((queue = new int[numVertices]) == NULL) || ((clusterNum = new int[numVertices]) == NULL))
        fatal("memory not allocated"); // allocate memory
    if (((globalClusterNum = new int[numVertices]) == NULL))
        fatal("memory not allocated"); // allocate memory

    for (int i = 0; i < numVertices; i++)
        visited[i] = 0; // initialize values to not visited

    for (int i = 0; i < numVertices; i++) {
        int degree = getDegree(i);

        if (degree == 0) {
            numSingle++;
            clusterNum[i] = -1; // designate node as singleton
            globalClusterNum[i] = -1;
            continue;
        }

        if (visited[i] == 0) { // unvisited vertex
            visited[i] = 1;
            numNonSingle++; // another nonsingleton cluster
            ptr = 0; // start a new component
            startPtr = 0; // start a new queue
            queue[0] = i; // first vertex in queue
            endPtr = 1;

            while (startPtr != endPtr) { // while queue is not empty

                int node = queue[startPtr++]; // remove first vertex on queue
                component[ptr++] = node; // this vertex is in component
                visited[node] = 1; // mark as visited

                Edge* edgePtr; // pointer to move through linked list of edges
                edgePtr = &(vertices[node].firstEdge); // point to first edge

                while (edgePtr->next != NULL) { // follow until find last edge

                    edgePtr = edgePtr->next; // pointer points at next edge in list
                    int endpt = edgePtr->target; // find adjacent vertex

                    if (!visited[endpt]) {
                        visited[endpt] = 1;
                        queue[endPtr++] = endpt; // add to queue
                    }
                }
            }

            long int numberEdges = 0; // count edges in component
            for (int j = 0; j < ptr; j++)
                numberEdges += getDegree(component[j]);
            numberEdges /= 2; // count (i,j) and (j,i) only once
            //totalEdges += numberEdges; // total edges in graph

            long double complete = (double)numberEdges; // compute completeness of cluster
            complete /= (long double)ptr * (long double)(ptr - 1);
            complete *= (double)(2.0);
            if (complete < 0 - TOL)
                fatal("Negative density of component computed");

            for (int j = 0; j < ptr; j++) { // assign cluster number
                clusterNum[component[j]] = k;
                globalClusterNum[component[j]] = k;
            }
            k++; // move to next cluster number

            if (compPtr == MAX_NUM_COMPS)
                fatal("MAX_NUM_COMPS exceeded.  Change to larger value in network.h");

            compSize[compPtr] = ptr; // record component size
            compDensity[compPtr++] = complete; // record density

            if (ptr == 2) // only 2 vertices in component
                numTwo++;

            else {
                numComp3orMore++; // one more with 3 or more vertices

                if (complete < 0.99999) // not a clique
                    numNotCliques++;

                for (int k = 0; k < 10000; k++)
                    if (k == ptr) {
                        num[k]++;
                        break;
                    }

                if (1) {
                    std::cout << numNonSingle << ": " << ptr << " vertices, " << numberEdges << " edges (";
                    std::cout << complete << " complete)" << std::endl;
                }

                if (0)
                    std::cout << ptr << "  " << complete << std::endl;

                if (DESCRIPTIVE_OUTPUT) {
                    fprintf(output, "%d vertices, %d edges (%.4f complete)\n", ptr, numberEdges, complete);
                    for (int j = 0; j < ptr; j++) {
                        std::cout << component[j] << " ";
                        fprintf(output, "%d ", component[j]);
                        if ((j + 1) % 10 == 0) // 10 vertices per line
                            fprintf(output, "\n");
                    }

                    std::cout << std::endl;
                    fprintf(output, "\n");
                }
            }

            if (ptr > maxCompAll)
                maxCompAll = ptr; // keep track of largest component overall

            if ((ptr > 2) && (BFS_GML || BFS_WG2)) // print out nontrivial components
                if (complete < MIN_COMPLETE + TOL) { // print only if not too complete
                    const char base[] = "comp"; // make up new filename for current comp.
                    const char suffix1[] = ".gml";
                    const char suffix2[] = ".nn";
                    const char suffix3[] = ".wg2";
                    char filename[50]; // .gml file
                    char filenn[50]; // node numbering file
                    char fileWG2[50]; //.wg2 file
                    printNum++; // increase to next cluster printed
                    //sprintf(filename, "%s%d_%d%s", base, randomVal, printNum, suffix1);
                    //sprintf(filenn, "%s%d_%d%s", base, randomVal, printNum, suffix2);
                    //sprintf(fileWG2, "%s%d_%d%s", base, randomVal, printNum, suffix3);

                    sprintf(filename, "%s%d%s", base, printNum, suffix1);
                    sprintf(filenn, "%s%d%s", base, printNum, suffix2);
                    sprintf(fileWG2, "%s%d%s", base, printNum, suffix3);

                    if (BFS_GML)
                        if ((outputGML = fopen(filename, "w")) == NULL)
                            fatal("Component output file could not be opened");

                    if ((outputCompnn = fopen(filenn, "w")) == NULL)
                        fatal("Component node numbering file could not be opened");

                    if (BFS_WG2)
                        if ((outputWG2 = fopen(fileWG2, "w")) == NULL)
                            fatal("Component output file could not be opened");

                    fprintf(outputCompnn, "%d nodes\n", ptr);

                    for (int j = 0; j < ptr; j++) // output node numbering
                        fprintf(outputCompnn, "%d ", component[j]);
                    fprintf(outputCompnn, "\n");
                    fclose(outputCompnn);

                    if (BFS_GML) {
                        fprintf(outputGML, "Graph with %d nodes and %d edges.\ngraph\n[", ptr, numberEdges);

                        for (int j = 0; j < ptr; j++)
                            fprintf(outputGML, "\n  node\n  [\n    id %d\n  ]", j + 1);
                        fclose(outputGML);

                        if ((outputGML = fopen(filename, "a")) == NULL)
                            fatal("Unable to open component output file");
                    }

                    if (ptr > maxSizeComp)
                        maxSizeComp = ptr; // keep track of largest component printed

                    int* nn; // array to hold node numbering
                    if ((nn = new int[numVertices]) == NULL)
                        fatal("memory not allocated"); // allocate memory

                    for (int j = 0; j < ptr; j++)
                        nn[j] = -1; // initialize
                    for (int j = 0; j < ptr; j++)
                        nn[component[j]] = j; // insert new node numbers

                    Edge* edgePtr; // pointer to move through linked list of edges

                    for (int j = 0; j < ptr; j++) { // find edges in component
                        int node = component[j];
                        edgePtr = &(vertices[node].firstEdge); // point first edge

                        if (nn[node] < 0)
                            fatal("Node numbering array error");
                        if (BFS_WG2)
                            fprintf(outputWG2, "%d %d, X %d [ ", j, j, vertices[node].degree);

                        while (edgePtr->next != NULL) { // follow until find last edge
                            edgePtr = edgePtr->next; // pointer points at next edge in list
                            if (nn[edgePtr->target] < 0)
                                fatal("node numbering array error");

                            if (BFS_GML)
                                if (!RETAINSYMMETRIC || (node < edgePtr->target)) // print once
                                    fprintf(outputGML, "\n\tedge\n\t[\n\tsource %d\n\ttarget %d\n\tweight %f\n\t]", nn[node] + 1, nn[edgePtr->target] + 1, edgePtr->weight);

                            if (BFS_WG2)
                                fprintf(outputWG2, "(%d, 1.0) ", nn[edgePtr->target]);
                        }

                        if (BFS_WG2)
                            fprintf(outputWG2, "] \n");
                    }

                    //for (int j = 0; j < ptr; j++)
                    //vertices[component[j]].printEdges(component[j], filename);

                    if (BFS_GML) {
                        fprintf(outputGML, "\n]\n");
                        fclose(outputGML);
                    }

                    if (BFS_WG2)
                        fclose(outputWG2);
                }
        }
    }

    //for (int j = 0; j < numVertices; j++)
    //std::cout << visited[j] << std::endl;

    //if (totalEdges != numEdges) fatal("Error in number of edges found in bfs");

    if (!DESCRIPTIVE_OUTPUT) {
        fprintf(output, "%d nodes %d clusters %d edges\n", numVertices, k, numEdges);
        for (int i = 0; i < numVertices; i++)
            fprintf(output, "%d ", clusterNum[i]);
        fprintf(output, "\n\n");
    }

    std::cout << std::endl;
    std::cout << numSingle << " vertices with degree zero" << std::endl;
    std::cout << numTwo << " components with only 2 vertices" << std::endl;
    std::cout << numComp3orMore << " components with 3 or more vertices" << std::endl;

    fprintf(output, "%d singletons, %d components with only 2 vertices,\n%d components with 3 or more vertices\n", numSingle, numTwo, numComp3orMore);

    // compute statistics about nonsingleton component sizes and densities
    float avgdata = 0.0; // first compute for component sizes
    float sdata = 0.0;
    float confUpData, confLowData;
    float tempData;
    float median;
    int size = compPtr;

    if (size > 1) {
        for (int i = 0; i < size; i++) // sort data
            for (int j = size - 1; j > i; j--)
                if (compSize[j] < compSize[j - 1]) {
                    tempData = compSize[j];
                    compSize[j] = compSize[j - 1];
                    compSize[j - 1] = (int)tempData;
                }

        if ((size % 2) == 0) // an even number of data
            median = 0.5 * (compSize[size / 2] + compSize[size / 2 - 1]);

        else // odd number of data
            median = compSize[size / 2];

        for (int i = 0; i < size; i++) // find mean of data
            avgdata += compSize[i];
        avgdata /= (float)size;

        for (int i = 0; i < size; i++) // find 95% confidence interval
            sdata += (compSize[i] - avgdata) * (compSize[i] - avgdata);
        sdata /= (float)((size - 1) * size);
        sdata = sqrt(sdata);
    }

    else {
        median = avgdata = compSize[0];
        sdata = 0;
    }

    confUpData = avgdata + (2.0 * sdata);
    confLowData = avgdata - (2.0 * sdata);

    std::cout << "\nSizes of non-singleton components:" << std::endl;
    std::cout << "\tmax = " << compSize[size - 1];
    std::cout << ", min = " << compSize[0];
    std::cout << ", median = " << median;
    std::cout << ", average = " << avgdata << std::endl;
    std::cout << "\t95% confidence interval: " << confUpData << ", " << confLowData << std::endl;

    float meanSize = avgdata; // save in case PRINTCOMPSIZES

    avgdata = 0.0; // repeat for component densities
    sdata = 0.0;
    if (size > 1) {
        for (int i = 0; i < size; i++) // sort data
            for (int j = size - 1; j > i; j--)
                if (compDensity[j] < compDensity[j - 1]) {
                    tempData = compDensity[j];
                    compDensity[j] = compDensity[j - 1];
                    compDensity[j - 1] = tempData;
                }

        if ((size % 2) == 0) // an even number of data
            median = 0.5 * (compDensity[size / 2] + compDensity[size / 2 - 1]);

        else // odd number of data
            median = compDensity[size / 2];

        for (int i = 0; i < size; i++) // find mean of data
            avgdata += compDensity[i];
        avgdata /= (float)size;

        for (int i = 0; i < size; i++) // find 95% confidence interval
            sdata += (compDensity[i] - avgdata) * (compDensity[i] - avgdata);
        sdata /= (float)((size - 1) * size);
        sdata = sqrt(sdata);
    }

    else {
        median = avgdata = compDensity[0];
        sdata = 0;
    }

    confUpData = avgdata + (2.0 * sdata);
    confLowData = avgdata - (2.0 * sdata);

    std::cout << "\nDensities of non-singleton components:" << std::endl;
    std::cout << "\tmax = " << compDensity[size - 1];
    std::cout << ", min = " << compDensity[0];
    std::cout << ", median = " << median;
    std::cout << ", average = " << avgdata << std::endl;
    std::cout << "\t95% confidence interval: " << confUpData << ", " << confLowData << std::endl;
    std::cout << std::endl;

    if (1) {
        num[1] = numSingle;
        num[2] = numTwo;
        for (int i = 1; i < 10000; i++)
            if (num[i] > 0)
                //std::cout << num[i] << std::endl;
                std::cout << "\t" << num[i] << " components with " << i << " vertices" << std::endl;
        //std::cout << num3 << " " << num4 << " " << num5 << " " << num6 << " " << num7;
        //std::cout << " " << num8 << " " << num9 << " " << num10 << std::endl;

        if (0)
            for (int i = 1; i < 64; i++)
                std::cout << num[i] << std::endl;
    }

    if (BFS_GML || BFS_WG2) {
        fprintf(output, "%d vertices in largest component printed to separate file\n", maxSizeComp);
        std::cout << "\n" << maxSizeComp << " vertices in largest component printed to separate file" << std::endl;
    }

    fprintf(output, "%d vertices in largest component overall\n", maxCompAll);
    std::cout << maxCompAll << " vertices in largest component overall" << std::endl;

    fprintf(output, "%d components are not cliques\n", numNotCliques);
    std::cout << numNotCliques << " components are not cliques" << std::endl;

    fclose(output);

    if (PRINTCOMPSIZES) { // print out summary of component sizes to "saveCompSizes.txt"
        FILE* compsize;

        // check if new file, if so print out header
        if ((compsize = fopen("saveCompSizes.txt", "r")) == NULL) {
            if ((compsize = fopen("saveCompSizes.txt", "w")) == NULL)
                fatal("'saveCompSizes.txt' file could not be opened.\n");

            fprintf(compsize, "NumSingle NumDouble Num3orMore NumNotCliques MaxSize AvgSize\n");
            fclose(compsize);
        }

        if ((compsize = fopen("saveCompSizes.txt", "a")) == NULL)
            fatal("'saveCompSizes.txt' file could not be opened.\n");

        fprintf(compsize, "%d %d %d %d %d %.2f\n", numSingle, numTwo, numComp3orMore, numNotCliques, maxCompAll, meanSize);
        fclose(compsize);
    }

    delete[] visited;
    delete[] component;
    delete[] queue;
    delete[] clusterNum;
}

// breadth-first search and Q_Value calculation, outputs connected components
// and respected Q_Value of communities and total q_value
double Network::q_calc(char* outputFile)
{
    double q_value = 0.0f;
    FILE* output;

    int testFile = system("test ! -f comp1.gml");
    if (testFile != 0)
        fatal("Component files already exist in this directory");

    if ((output = fopen(outputFile, "w")) == NULL)
        fatal("Unable to open output file");

    int ptr; // point to current index of component array
    int startPtr, endPtr; // pointers for queue
    int *visited, *component, *queue, *clusterNum;
    int numSingle = 0; // count nodes with degree zero
    int numTwo = 0; // count number of components with 2 vertices
    int k = 0; // current cluster number
    kcluster = &k;
    // long int totalEdges = 0; // total edges in graph
    int numComp3orMore = 0; // number of components with 3 or more vertices
    FILE* outputGML; // output components to .gml files
    FILE* outputWG2; // output components to .wg2 files
    FILE* outputCompnn; // output node numbering for components
    int randomVal; // random value to place at end of comp file names
    int numNonSingle = 0; // number of non-singleton components
    int printNum = 0; // number of cluster printed out
    int maxSizeComp = 0; // max size of component printed to gml file
    int maxCompAll = 0; // max size of component overall
    int num[10000]; // hold number of components with each size
    int numNotCliques = 0; // hold number of components that aren't cliques

    for (int i = 0; i < 10000; i++)
        num[i] = 0;

    int compSize[MAX_NUM_COMPS]; // hold non-singleton component sizes
    double compDensity[MAX_NUM_COMPS]; // hold non-singleton component densities
    int compPtr = 0; // pointer for these two arrays

    if (BFS_WG2 && !RETAINSYMMETRIC)
        fatal("Need to retain symmetric edges for .wg2 output");

    // srand(time(NULL)); // assign random value between 1 - 1000
    randomVal = rand() % 1000 + 1;

    if (((visited = new int[numVertices]) == NULL) || ((component = new int[numVertices]) == NULL))
        fatal("memory not allocated"); // allocate memory
    if (((queue = new int[numVertices]) == NULL) || ((clusterNum = new int[numVertices]) == NULL))
        fatal("memory not allocated"); // allocate memory

    for (int i = 0; i < numVertices; i++)
        visited[i] = 0; // initialize values to not visited

    for (int i = 0; i < numVertices; i++) {
        int degree = getDegree(i);

        if (degree == 0) {
            numSingle++;
            clusterNum[i] = -1; // designate node as singleton
            continue;
        }

        if (visited[i] == 0) { // unvisited vertex
            visited[i] = 1;
            numNonSingle++; // another nonsingleton cluster
            ptr = 0; // start a new component
            startPtr = 0; // start a new queue
            queue[0] = i; // first vertex in queue
            endPtr = 1;

            while (startPtr != endPtr) { // while queue is not empty
                int node = queue[startPtr++]; // remove first vertex on queue
                component[ptr++] = node; // this vertex is in component
                visited[node] = 1; // mark as visited

                Edge* edgePtr; // pointer to move through linked list of edges
                edgePtr = &(vertices[node].firstEdge); // point to first edge

                while (edgePtr->next != NULL) { // follow until find last edge
                    edgePtr = edgePtr->next; // pointer points at next edge in list
                    int endpt = edgePtr->target; // find adjacent vertex

                    if (!visited[endpt]) {
                        visited[endpt] = 1;
                        queue[endPtr++] = endpt; // add to queue
                    }
                }
            }

            long int numberEdges = 0; // count edges in component
            for (int j = 0; j < ptr; j++)
                numberEdges += getDegree(component[j]);
            numberEdges /= 2; // count (i,j) and (j,i) only once
            // totalEdges += numberEdges; // total edges in graph

            long double complete = (double)numberEdges; // compute completeness of cluster
            complete /= (long double)ptr * (long double)(ptr - 1);
            complete *= (double)(2.0);
            if (complete < 0 - TOL)
                fatal("Negative density of component computed");

            for (int j = 0; j < ptr; j++) { // assign cluster number
                clusterNum[component[j]] = k;
            }
            k++; // move to next cluster number

            if (compPtr == MAX_NUM_COMPS)
                fatal("MAX_NUM_COMPS exceeded.  Change to larger value in network.h");

            compSize[compPtr] = ptr; // record component size
            compDensity[compPtr++] = complete; // record density

            if (ptr == 2) // only 2 vertices in component
                numTwo++;

            else {
                numComp3orMore++; // one more with 3 or more vertices

                if (complete < 0.99999) // not a clique
                    numNotCliques++;

                for (int k = 0; k < 10000; k++)
                    if (k == ptr) {
                        num[k]++;
                        break;
                    }

                // At this point we have the vertices and edges in a community.
                // Calcluating the q_value here

                if(0) {
                    std ::cout << "--------------" << numberEdges << "," << originalNumEdges << "," << ptr << "," << originalNumVertices << std ::endl;
                }

                float ci = (float)numberEdges / originalNumEdges;
                float ri = ((float)ptr * (ptr - 1)) / (originalNumVertices * (originalNumVertices - 1));

                if(0) {
                    std :: cout << "ci = " << ci << " , ri = " << ri << std :: endl;
                }

                q_value += ci - ri;

                if (1) {
                    std::cout << numNonSingle << ": " << ptr << " vertices, " << numberEdges << " edges (";
                    std::cout << complete << " complete), q_value : " << q_value << std::endl;
                }

                if (0)
                    std::cout << ptr << "  " << complete << std::endl;

                if (DESCRIPTIVE_OUTPUT) {
                    fprintf(output, "%d vertices, %d edges (%.4f complete)\n", ptr, numberEdges, complete);
                    for (int j = 0; j < ptr; j++) {
                        std::cout << component[j] << " ";
                        fprintf(output, "%d ", component[j]);
                        if ((j + 1) % 10 == 0) // 10 vertices per line
                            fprintf(output, "\n");
                    }

                    std::cout << std::endl;
                    fprintf(output, "\n");
                }
            }

            if (ptr > maxCompAll)
                maxCompAll = ptr; // keep track of largest component overall

            if ((ptr > 2) && (BFS_GML || BFS_WG2)) // print out nontrivial components
                if (complete < MIN_COMPLETE + TOL) { // print only if not too complete
                    const char base[] = "comp"; // make up new filename for current comp.
                    const char suffix1[] = ".gml";
                    const char suffix2[] = ".nn";
                    const char suffix3[] = ".wg2";
                    char filename[50]; // .gml file
                    char filenn[50]; // node numbering file
                    char fileWG2[50]; //.wg2 file
                    printNum++; // increase to next cluster printed
                    //sprintf(filename, "%s%d_%d%s", base, randomVal, printNum, suffix1);
                    //sprintf(filenn, "%s%d_%d%s", base, randomVal, printNum, suffix2);
                    //sprintf(fileWG2, "%s%d_%d%s", base, randomVal, printNum, suffix3);

                    sprintf(filename, "%s%d%s", base, printNum, suffix1);
                    sprintf(filenn, "%s%d%s", base, printNum, suffix2);
                    sprintf(fileWG2, "%s%d%s", base, printNum, suffix3);

                    if (BFS_GML)
                        if ((outputGML = fopen(filename, "w")) == NULL)
                            fatal("Component output file could not be opened");

                    if ((outputCompnn = fopen(filenn, "w")) == NULL)
                        fatal("Component node numbering file could not be opened");

                    if (BFS_WG2)
                        if ((outputWG2 = fopen(fileWG2, "w")) == NULL)
                            fatal("Component output file could not be opened");

                    fprintf(outputCompnn, "%d nodes\n", ptr);

                    for (int j = 0; j < ptr; j++) // output node numbering
                        fprintf(outputCompnn, "%d ", component[j]);
                    fprintf(outputCompnn, "\n");
                    fclose(outputCompnn);

                    if (BFS_GML) {
                        fprintf(outputGML, "Graph with %d nodes and %d edges.\ngraph\n[", ptr, numberEdges);

                        for (int j = 0; j < ptr; j++)
                            fprintf(outputGML, "\n  node\n  [\n    id %d\n  ]", j + 1);
                        fclose(outputGML);

                        if ((outputGML = fopen(filename, "a")) == NULL)
                            fatal("Unable to open component output file");
                    }

                    if (ptr > maxSizeComp)
                        maxSizeComp = ptr; // keep track of largest component printed

                    int* nn; // array to hold node numbering
                    if ((nn = new int[numVertices]) == NULL)
                        fatal("memory not allocated"); // allocate memory

                    for (int j = 0; j < ptr; j++)
                        nn[j] = -1; // initialize
                    for (int j = 0; j < ptr; j++)
                        nn[component[j]] = j; // insert new node numbers

                    Edge* edgePtr; // pointer to move through linked list of edges

                    for (int j = 0; j < ptr; j++) { // find edges in component
                        int node = component[j];
                        edgePtr = &(vertices[node].firstEdge); // point first edge

                        if (nn[node] < 0)
                            fatal("Node numbering array error");
                        if (BFS_WG2)
                            fprintf(outputWG2, "%d %d, X %d [ ", j, j, vertices[node].degree);

                        while (edgePtr->next != NULL) { // follow until find last edge
                            edgePtr = edgePtr->next; // pointer points at next edge in list
                            if (nn[edgePtr->target] < 0)
                                fatal("node numbering array error");

                            if (BFS_GML)
                                if (!RETAINSYMMETRIC || (node < edgePtr->target)) // print once
                                    fprintf(outputGML, "\n\tedge\n\t[\n\tsource %d\n\ttarget %d\n\tweight %f\n\t]", nn[node] + 1, nn[edgePtr->target] + 1, edgePtr->weight);

                            if (BFS_WG2)
                                fprintf(outputWG2, "(%d, 1.0) ", nn[edgePtr->target]);
                        }

                        if (BFS_WG2)
                            fprintf(outputWG2, "] \n");
                    }

                    //for (int j = 0; j < ptr; j++)
                    //vertices[component[j]].printEdges(component[j], filename);

                    if (BFS_GML) {
                        fprintf(outputGML, "\n]\n");
                        fclose(outputGML);
                    }

                    if (BFS_WG2)
                        fclose(outputWG2);
                }
        }
    }

    //for (int j = 0; j < numVertices; j++)
    //std::cout << visited[j] << std::endl;

    //if (totalEdges != numEdges) fatal("Error in number of edges found in bfs");

    if (!DESCRIPTIVE_OUTPUT) {
        fprintf(output, "%d nodes %d clusters %d edges\n", numVertices, k, numEdges);
        for (int i = 0; i < numVertices; i++)
            fprintf(output, "%d ", clusterNum[i]);
        fprintf(output, "\n\n");
    }

    std::cout << std::endl;
    std::cout << numSingle << " vertices with degree zero" << std::endl;
    std::cout << numTwo << " components with only 2 vertices" << std::endl;
    std::cout << numComp3orMore << " components with 3 or more vertices" << std::endl;

    std::cout << "Total Qs value for the given network : " << q_value << std::endl;
    std::cout << "Total communities in the network : " << k << std::endl;

    fprintf(output, "%d singletons, %d components with only 2 vertices,\n%d components with 3 or more vertices\n", numSingle, numTwo, numComp3orMore);

    // compute statistics about nonsingleton component sizes and densities
    float avgdata = 0.0; // first compute for component sizes
    float sdata = 0.0;
    float confUpData, confLowData;
    float tempData;
    float median;
    int size = compPtr;

    if (size > 1) {
        for (int i = 0; i < size; i++) // sort data
            for (int j = size - 1; j > i; j--)
                if (compSize[j] < compSize[j - 1]) {
                    tempData = compSize[j];
                    compSize[j] = compSize[j - 1];
                    compSize[j - 1] = (int)tempData;
                }

        if ((size % 2) == 0) // an even number of data
            median = 0.5 * (compSize[size / 2] + compSize[size / 2 - 1]);

        else // odd number of data
            median = compSize[size / 2];

        for (int i = 0; i < size; i++) // find mean of data
            avgdata += compSize[i];
        avgdata /= (float)size;

        for (int i = 0; i < size; i++) // find 95% confidence interval
            sdata += (compSize[i] - avgdata) * (compSize[i] - avgdata);
        sdata /= (float)((size - 1) * size);
        sdata = sqrt(sdata);
    }

    else {
        median = avgdata = compSize[0];
        sdata = 0;
    }

    confUpData = avgdata + (2.0 * sdata);
    confLowData = avgdata - (2.0 * sdata);

    std::cout << "\nSizes of non-singleton components:" << std::endl;
    std::cout << "\tmax = " << compSize[size - 1];
    std::cout << ", min = " << compSize[0];
    std::cout << ", median = " << median;
    std::cout << ", average = " << avgdata << std::endl;
    std::cout << "\t95% confidence interval: " << confUpData << ", " << confLowData << std::endl;

    float meanSize = avgdata; // save in case PRINTCOMPSIZES

    avgdata = 0.0; // repeat for component densities
    sdata = 0.0;
    if (size > 1) {
        for (int i = 0; i < size; i++) // sort data
            for (int j = size - 1; j > i; j--)
                if (compDensity[j] < compDensity[j - 1]) {
                    tempData = compDensity[j];
                    compDensity[j] = compDensity[j - 1];
                    compDensity[j - 1] = tempData;
                }

        if ((size % 2) == 0) // an even number of data
            median = 0.5 * (compDensity[size / 2] + compDensity[size / 2 - 1]);

        else // odd number of data
            median = compDensity[size / 2];

        for (int i = 0; i < size; i++) // find mean of data
            avgdata += compDensity[i];
        avgdata /= (float)size;

        for (int i = 0; i < size; i++) // find 95% confidence interval
            sdata += (compDensity[i] - avgdata) * (compDensity[i] - avgdata);
        sdata /= (float)((size - 1) * size);
        sdata = sqrt(sdata);
    }

    else {
        median = avgdata = compDensity[0];
        sdata = 0;
    }

    confUpData = avgdata + (2.0 * sdata);
    confLowData = avgdata - (2.0 * sdata);

    std::cout << "\nDensities of non-singleton components:" << std::endl;
    std::cout << "\tmax = " << compDensity[size - 1];
    std::cout << ", min = " << compDensity[0];
    std::cout << ", median = " << median;
    std::cout << ", average = " << avgdata << std::endl;
    std::cout << "\t95% confidence interval: " << confUpData << ", " << confLowData << std::endl;
    std::cout << std::endl;

    if (1) {
        num[1] = numSingle;
        num[2] = numTwo;
        for (int i = 1; i < 10000; i++)
            if (num[i] > 0)
                //std::cout << num[i] << std::endl;
                std::cout << "\t" << num[i] << " components with " << i << " vertices" << std::endl;
        //std::cout << num3 << " " << num4 << " " << num5 << " " << num6 << " " << num7;
        //std::cout << " " << num8 << " " << num9 << " " << num10 << std::endl;

        if (0)
            for (int i = 1; i < 64; i++)
                std::cout << num[i] << std::endl;
    }

    if (BFS_GML || BFS_WG2) {
        fprintf(output, "%d vertices in largest component printed to separate file\n", maxSizeComp);
        std::cout << "\n" << maxSizeComp << " vertices in largest component printed to separate file" << std::endl;
    }

    fprintf(output, "%d vertices in largest component overall\n", maxCompAll);
    std::cout << maxCompAll << " vertices in largest component overall" << std::endl;

    fprintf(output, "%d components are not cliques\n", numNotCliques);
    std::cout << numNotCliques << " components are not cliques" << std::endl;

    fclose(output);

    if (PRINTCOMPSIZES) { // print out summary of component sizes to "saveCompSizes.txt"
        FILE* compsize;

        // check if new file, if so print out header
        if ((compsize = fopen("saveCompSizes.txt", "r")) == NULL) {
            if ((compsize = fopen("saveCompSizes.txt", "w")) == NULL)
                fatal("'saveCompSizes.txt' file could not be opened.\n");

            fprintf(compsize, "NumSingle NumDouble Num3orMore NumNotCliques MaxSize AvgSize\n");
            fclose(compsize);
        }

        if ((compsize = fopen("saveCompSizes.txt", "a")) == NULL)
            fatal("'saveCompSizes.txt' file could not be opened.\n");

        fprintf(compsize, "%d %d %d %d %d %.2f\n", numSingle, numTwo, numComp3orMore, numNotCliques, maxCompAll, meanSize);
        fclose(compsize);
    }

    delete[] visited;
    delete[] component;
    delete[] queue;
    delete[] clusterNum;

    return q_value;
}

// breadth-first search and Q_Value calculation, outputs connected components
// and respected Q_Value of communities and total q_value
double Network::modularity(char *outputFile)
{
    double modularity_value = 0.0f;
    FILE *output;

    int testFile = system("test ! -f comp1.gml");
    if (testFile != 0)
        fatal("Component files already exist in this directory");

    if ((output = fopen(outputFile, "w")) == NULL)
        fatal("Unable to open output file");

    int ptr;              // point to current index of component array
    int startPtr, endPtr; // pointers for queue
    int *visited, *component, *queue, *clusterNum;
    int numSingle = 0; // count nodes with degree zero
    int numTwo = 0;    // count number of components with 2 vertices
    int k = 0;         // current cluster number
    kcluster = &k; // point kcluster to k
    // long int totalEdges = 0; // total edges in graph
    int numComp3orMore = 0; // number of components with 3 or more vertices
    FILE *outputGML;        // output components to .gml files
    FILE *outputWG2;        // output components to .wg2 files
    FILE *outputCompnn;     // output node numbering for components
    int randomVal;          // random value to place at end of comp file names
    int numNonSingle = 0;   // number of non-singleton components
    int printNum = 0;       // number of cluster printed out
    int maxSizeComp = 0;    // max size of component printed to gml file
    int maxCompAll = 0;     // max size of component overall
    int num[10000];         // hold number of components with each size
    int numNotCliques = 0;  // hold number of components that aren't cliques

    for (int i = 0; i < 10000; i++)
        num[i] = 0;

    int compSize[MAX_NUM_COMPS];       // hold non-singleton component sizes
    double compDensity[MAX_NUM_COMPS]; // hold non-singleton component densities
    int compPtr = 0;                   // pointer for these two arrays

    if (BFS_WG2 && !RETAINSYMMETRIC)
        fatal("Need to retain symmetric edges for .wg2 output");

    // srand(time(NULL)); // assign random value between 1 - 1000
    randomVal = rand() % 1000 + 1;

    if (((visited = new int[numVertices]) == NULL) || ((component = new int[numVertices]) == NULL))
        fatal("memory not allocated"); // allocate memory
    if (((queue = new int[numVertices]) == NULL) || ((clusterNum = new int[numVertices]) == NULL))
        fatal("memory not allocated"); // allocate memory
    if (((globalClusterNum = new int[numVertices]) == NULL))
        fatal("memory not allocated"); // allocate memory

    for (int i = 0; i < numVertices; i++)
        visited[i] = 0; // initialize values to not visited

    for (int i = 0; i < numVertices; i++)
    {
        int degree = getDegree(i);

        if (degree == 0)
        {
            numSingle++;
            clusterNum[i] = -1; // designate node as singleton
            globalClusterNum[i] = -1;
            continue;
        }

        if (visited[i] == 0)
        { // unvisited vertex
            visited[i] = 1;
            numNonSingle++; // another nonsingleton cluster
            ptr = 0;        // start a new component
            startPtr = 0;   // start a new queue
            queue[0] = i;   // first vertex in queue
            endPtr = 1;

            while (startPtr != endPtr)
            {                                 // while queue is not empty
                int node = queue[startPtr++]; // remove first vertex on queue
                component[ptr++] = node;      // this vertex is in component
                visited[node] = 1;            // mark as visited

                Edge *edgePtr;                         // pointer to move through linked list of edges
                edgePtr = &(vertices[node].firstEdge); // point to first edge

                while (edgePtr->next != NULL)
                {                                // follow until find last edge
                    edgePtr = edgePtr->next;     // pointer points at next edge in list
                    int endpt = edgePtr->target; // find adjacent vertex

                    if (!visited[endpt])
                    {
                        visited[endpt] = 1;
                        queue[endPtr++] = endpt; // add to queue
                    }
                }
            }

            long int numberEdges = 0; // count edges in component
            for (int j = 0; j < ptr; j++)
                numberEdges += getDegree(component[j]);
            numberEdges /= 2; // count (i,j) and (j,i) only once
            // totalEdges += numberEdges; // total edges in graph

            long double complete = (double)numberEdges; // compute completeness of cluster
            complete /= (long double)ptr * (long double)(ptr - 1);
            complete *= (double)(2.0);
            if (complete < 0 - TOL)
                fatal("Negative density of component computed");

            for (int j = 0; j < ptr; j++) {// assign cluster number
                clusterNum[component[j]] = k;
                globalClusterNum[component[j]] = k;
            }
            k++; // move to next cluster number

            if (compPtr == MAX_NUM_COMPS)
                fatal("MAX_NUM_COMPS exceeded.  Change to larger value in network.h");

            compSize[compPtr] = ptr;           // record component size
            compDensity[compPtr++] = complete; // record density

            double sum_degree = 0.0f;
            for (int j = 0; j < ptr; j++)
                sum_degree += getOriginalDegree(component[j]);

            float first = (float)numberEdges / originalNumEdges;
            float second = pow((float)sum_degree / (2 * originalNumEdges), 2);

            if (0)
            {
                std ::cout << "--------------" << numberEdges << "," << originalNumEdges << "," << ptr << "," << originalNumVertices << std ::endl;
                for (int j = 0; j < ptr; j++) // loop through the components in cluster
                    std ::cout << component[j] << ",";
                std ::cout << "-=-" << sum_degree;
                std ::cout << "\nci = " << first << " , ri = " << second << std ::endl;
                std ::cout << std ::endl;
            }

            modularity_value += first - second;

            if (1)
            {
                std::cout << numNonSingle << ": " << ptr << " vertices, " << numberEdges << " edges (";
                std::cout << complete << " complete), modularity_value : " << first - second << std::endl;
            }

            if (ptr == 2) { // only 2 vertices in component
                numTwo++;
            }
            else
            {
                numComp3orMore++; // one more with 3 or more vertices

                if (complete < 0.99999) // not a clique
                    numNotCliques++;

                for (int k = 0; k < 10000; k++)
                    if (k == ptr)
                    {
                        num[k]++;
                        break;
                    }

                // At this point we have the vertices and edges in a community.
                // Calcluating the q_value her

                if (0)
                    std::cout << ptr << "  " << complete << std::endl;

                if (DESCRIPTIVE_OUTPUT)
                {
                    fprintf(output, "%d vertices, %d edges (%.4f complete)\n", ptr, numberEdges, complete);
                    for (int j = 0; j < ptr; j++)
                    {
                        std::cout << component[j] << " ";
                        fprintf(output, "%d ", component[j]);
                        if ((j + 1) % 10 == 0) // 10 vertices per line
                            fprintf(output, "\n");
                    }

                    std::cout << std::endl;
                    fprintf(output, "\n");
                }
            }

            if (ptr > maxCompAll)
                maxCompAll = ptr; // keep track of largest component overall

            if ((ptr > 2) && (BFS_GML || BFS_WG2)) // print out nontrivial components
                if (complete < MIN_COMPLETE + TOL)
                {                               // print only if not too complete
                    const char base[] = "comp"; // make up new filename for current comp.
                    const char suffix1[] = ".gml";
                    const char suffix2[] = ".nn";
                    const char suffix3[] = ".wg2";
                    char filename[50]; // .gml file
                    char filenn[50];   // node numbering file
                    char fileWG2[50];  //.wg2 file
                    printNum++;        // increase to next cluster printed
                    //sprintf(filename, "%s%d_%d%s", base, randomVal, printNum, suffix1);
                    //sprintf(filenn, "%s%d_%d%s", base, randomVal, printNum, suffix2);
                    //sprintf(fileWG2, "%s%d_%d%s", base, randomVal, printNum, suffix3);

                    sprintf(filename, "%s%d%s", base, printNum, suffix1);
                    sprintf(filenn, "%s%d%s", base, printNum, suffix2);
                    sprintf(fileWG2, "%s%d%s", base, printNum, suffix3);

                    if (BFS_GML)
                        if ((outputGML = fopen(filename, "w")) == NULL)
                            fatal("Component output file could not be opened");

                    if ((outputCompnn = fopen(filenn, "w")) == NULL)
                        fatal("Component node numbering file could not be opened");

                    if (BFS_WG2)
                        if ((outputWG2 = fopen(fileWG2, "w")) == NULL)
                            fatal("Component output file could not be opened");

                    fprintf(outputCompnn, "%d nodes\n", ptr);

                    for (int j = 0; j < ptr; j++) // output node numbering
                        fprintf(outputCompnn, "%d ", component[j]);
                    fprintf(outputCompnn, "\n");
                    fclose(outputCompnn);

                    if (BFS_GML)
                    {
                        fprintf(outputGML, "Graph with %d nodes and %d edges.\ngraph\n[", ptr, numberEdges);

                        for (int j = 0; j < ptr; j++)
                            fprintf(outputGML, "\n  node\n  [\n    id %d\n  ]", j + 1);
                        fclose(outputGML);

                        if ((outputGML = fopen(filename, "a")) == NULL)
                            fatal("Unable to open component output file");
                    }

                    if (ptr > maxSizeComp)
                        maxSizeComp = ptr; // keep track of largest component printed

                    int *nn; // array to hold node numbering
                    if ((nn = new int[numVertices]) == NULL)
                        fatal("memory not allocated"); // allocate memory

                    for (int j = 0; j < ptr; j++)
                        nn[j] = -1; // initialize
                    for (int j = 0; j < ptr; j++)
                        nn[component[j]] = j; // insert new node numbers

                    Edge *edgePtr; // pointer to move through linked list of edges

                    for (int j = 0; j < ptr; j++)
                    { // find edges in component
                        int node = component[j];
                        edgePtr = &(vertices[node].firstEdge); // point first edge

                        if (nn[node] < 0)
                            fatal("Node numbering array error");
                        if (BFS_WG2)
                            fprintf(outputWG2, "%d %d, X %d [ ", j, j, vertices[node].degree);

                        while (edgePtr->next != NULL)
                        {                            // follow until find last edge
                            edgePtr = edgePtr->next; // pointer points at next edge in list
                            if (nn[edgePtr->target] < 0)
                                fatal("node numbering array error");

                            if (BFS_GML)
                                if (!RETAINSYMMETRIC || (node < edgePtr->target)) // print once
                                    fprintf(outputGML, "\n\tedge\n\t[\n\tsource %d\n\ttarget %d\n\tweight %f\n\t]", nn[node] + 1, nn[edgePtr->target] + 1, edgePtr->weight);

                            if (BFS_WG2)
                                fprintf(outputWG2, "(%d, 1.0) ", nn[edgePtr->target]);
                        }

                        if (BFS_WG2)
                            fprintf(outputWG2, "] \n");
                    }

                    //for (int j = 0; j < ptr; j++)
                    //vertices[component[j]].printEdges(component[j], filename);

                    if (BFS_GML)
                    {
                        fprintf(outputGML, "\n]\n");
                        fclose(outputGML);
                    }

                    if (BFS_WG2)
                        fclose(outputWG2);
                }
        }
    }

    //for (int j = 0; j < numVertices; j++)
    //std::cout << visited[j] << std::endl;

    //if (totalEdges != numEdges) fatal("Error in number of edges found in bfs");

    if (!DESCRIPTIVE_OUTPUT)
    {
        fprintf(output, "%d nodes %d clusters %d edges\n", numVertices, k, numEdges);
        for (int i = 0; i < numVertices; i++)
            fprintf(output, "%d ", clusterNum[i]);
        fprintf(output, "\n\n");
    }

    std::cout << std::endl;
    std::cout << numSingle << " vertices with degree zero" << std::endl;
    std::cout << numTwo << " components with only 2 vertices" << std::endl;
    std::cout << numComp3orMore << " components with 3 or more vertices" << std::endl;

    std::cout << "Total Modularity value for the given network : " << modularity_value << std::endl;
    std::cout << "Total communities in the network : " << k << std::endl;

    fprintf(output, "%d singletons, %d components with only 2 vertices,\n%d components with 3 or more vertices\n", numSingle, numTwo, numComp3orMore);

    // compute statistics about nonsingleton component sizes and densities
    float avgdata = 0.0; // first compute for component sizes
    float sdata = 0.0;
    float confUpData, confLowData;
    float tempData;
    float median;
    int size = compPtr;

    if (size > 1)
    {
        for (int i = 0; i < size; i++) // sort data
            for (int j = size - 1; j > i; j--)
                if (compSize[j] < compSize[j - 1])
                {
                    tempData = compSize[j];
                    compSize[j] = compSize[j - 1];
                    compSize[j - 1] = (int)tempData;
                }

        if ((size % 2) == 0) // an even number of data
            median = 0.5 * (compSize[size / 2] + compSize[size / 2 - 1]);

        else // odd number of data
            median = compSize[size / 2];

        for (int i = 0; i < size; i++) // find mean of data
            avgdata += compSize[i];
        avgdata /= (float)size;

        for (int i = 0; i < size; i++) // find 95% confidence interval
            sdata += (compSize[i] - avgdata) * (compSize[i] - avgdata);
        sdata /= (float)((size - 1) * size);
        sdata = sqrt(sdata);
    }

    else
    {
        median = avgdata = compSize[0];
        sdata = 0;
    }

    confUpData = avgdata + (2.0 * sdata);
    confLowData = avgdata - (2.0 * sdata);

    std::cout << "\nSizes of non-singleton components:" << std::endl;
    std::cout << "\tmax = " << compSize[size - 1];
    std::cout << ", min = " << compSize[0];
    std::cout << ", median = " << median;
    std::cout << ", average = " << avgdata << std::endl;
    std::cout << "\t95% confidence interval: " << confUpData << ", " << confLowData << std::endl;

    float meanSize = avgdata; // save in case PRINTCOMPSIZES

    avgdata = 0.0; // repeat for component densities
    sdata = 0.0;
    if (size > 1)
    {
        for (int i = 0; i < size; i++) // sort data
            for (int j = size - 1; j > i; j--)
                if (compDensity[j] < compDensity[j - 1])
                {
                    tempData = compDensity[j];
                    compDensity[j] = compDensity[j - 1];
                    compDensity[j - 1] = tempData;
                }

        if ((size % 2) == 0) // an even number of data
            median = 0.5 * (compDensity[size / 2] + compDensity[size / 2 - 1]);

        else // odd number of data
            median = compDensity[size / 2];

        for (int i = 0; i < size; i++) // find mean of data
            avgdata += compDensity[i];
        avgdata /= (float)size;

        for (int i = 0; i < size; i++) // find 95% confidence interval
            sdata += (compDensity[i] - avgdata) * (compDensity[i] - avgdata);
        sdata /= (float)((size - 1) * size);
        sdata = sqrt(sdata);
    }

    else
    {
        median = avgdata = compDensity[0];
        sdata = 0;
    }

    confUpData = avgdata + (2.0 * sdata);
    confLowData = avgdata - (2.0 * sdata);

    std::cout << "\nDensities of non-singleton components:" << std::endl;
    std::cout << "\tmax = " << compDensity[size - 1];
    std::cout << ", min = " << compDensity[0];
    std::cout << ", median = " << median;
    std::cout << ", average = " << avgdata << std::endl;
    std::cout << "\t95% confidence interval: " << confUpData << ", " << confLowData << std::endl;
    std::cout << std::endl;

    if (1)
    {
        num[1] = numSingle;
        num[2] = numTwo;
        for (int i = 1; i < 10000; i++)
            if (num[i] > 0)
                //std::cout << num[i] << std::endl;
                std::cout << "\t" << num[i] << " components with " << i << " vertices" << std::endl;
        //std::cout << num3 << " " << num4 << " " << num5 << " " << num6 << " " << num7;
        //std::cout << " " << num8 << " " << num9 << " " << num10 << std::endl;

        if (0)
            for (int i = 1; i < 64; i++)
                std::cout << num[i] << std::endl;
    }

    if (BFS_GML || BFS_WG2)
    {
        fprintf(output, "%d vertices in largest component printed to separate file\n", maxSizeComp);
        std::cout << "\n"
                  << maxSizeComp << " vertices in largest component printed to separate file" << std::endl;
    }

    fprintf(output, "%d vertices in largest component overall\n", maxCompAll);
    std::cout << maxCompAll << " vertices in largest component overall" << std::endl;

    fprintf(output, "%d components are not cliques\n", numNotCliques);
    std::cout << numNotCliques << " components are not cliques" << std::endl;

    fclose(output);

    if (PRINTCOMPSIZES)
    { // print out summary of component sizes to "saveCompSizes.txt"
        FILE *compsize;

        // check if new file, if so print out header
        if ((compsize = fopen("saveCompSizes.txt", "r")) == NULL)
        {
            if ((compsize = fopen("saveCompSizes.txt", "w")) == NULL)
                fatal("'saveCompSizes.txt' file could not be opened.\n");

            fprintf(compsize, "NumSingle NumDouble Num3orMore NumNotCliques MaxSize AvgSize\n");
            fclose(compsize);
        }

        if ((compsize = fopen("saveCompSizes.txt", "a")) == NULL)
            fatal("'saveCompSizes.txt' file could not be opened.\n");

        fprintf(compsize, "%d %d %d %d %d %.2f\n", numSingle, numTwo, numComp3orMore, numNotCliques, maxCompAll, meanSize);
        fclose(compsize);
    }

    delete[] visited;
    delete[] component;
    delete[] queue;
    delete[] clusterNum;

    return modularity_value;
}