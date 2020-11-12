# ERBGA
-----------------



Project Structure:

```
bfs-orig|->
	|README
	|README.md
	|calc.cpp - Implementation of the calculations
	|datasets.txt - Datasets to be run on Lewis
	|erbga_driver.cpp - Read and Populates Network and Runs the GA
	|erbga_driver.h - Control flags for Debugging, Output
	|erbga.cpp - Genetic Algorithm Core
	|erbga.h - Data Objects Definition, Parameters, Flags for GA
	|network.cpp - Network core class
	|network.h - Vertex, Edge, Network class Definitions
	|run_ga.sh - Script to run each dataset read from dataset.txt a pecified number of epochs
	|run_multiple_islands.sh - dependency file for run_ga.sh
	|run_sbatch_ga.sh - Core sbatch file used to run the program
	|timer.h - Timer Core file
	|testFile.gml - a Test GML file with 15 nodes and 15 edges
	|Makefile
		-|clean - cleans and object and erbga binary
		-|cleanlog - cleans *.log, log files
```

Compiling:

```
$ make
```

Running the project

```
$ ./erbga input.gml out.bfs | tee input-screen
```

Help:

```
$ ./erbga
Usage : ./erbga -n processes -k loops -m timetosleep
n = number of process to fork
k = number of times to loop on fprintf
m = time to sleep in millis
```

Clean the project *.obj and erbga binary:

```
$ make clean
```

Clean the project only *.log:

```
$ make cleanlog
```

## Network Library - Sharlee Climer

October/November 2009.

Breadth-First Search (BFS) is a simple exploration of a network, where each
connected component is identified.  The nodes in each component are output.

--

July 2017 - March 2018

Modifications added to the Network Library by Aditya Karnam:

* Add New Edge to the front of the list
* Moved ID and invID from BFS to Network
* assignID and getID functions added
* Added remove edge
* Fixed garbage collection
* Added Qs - Fitness Calculation    
* Added Modularity - Standard Fitness
* Added compatibility data objects to interface with ERBGA

--

### Psuedocode:

----
```
initialize an array visited[n] to all 0’s

for each node v
  if (visited[v] = 0)
    put v in a queue 
    visited[v] = 1
    while (queue is not empty)
      -remove node w from queue 
      -add w to component
      -add all unvisited nodes that are incident to w to queue and mark as visited
    record component
```
----

* The input file must be in .gml format.  This common format is described in the 'formatSummary' file that can be found in this library.  It is a list of the nodes and edges of a network.  

* The main output file has a suffix '.bfs' and a custom format, as described in 'formatSummary'.

* Files in .gml format can also be output for each component.  Components with densities that are less than a threshold are output.  The density is equal to the number of edges divided by the number of edges possible.  Adjust MIN\_COMPLETE in 'network.h' to the desired minimum density.  The component will be printed out only if it has a density that is less than MIN\_COMPLETE.  (We used 0.5 as a default value.)

* These component output files are named 'compx.gml' with 'x' equal to the number of component as it is written out.  These numbers start with zero and are consecutive.  Each run will produce new compx.gml files, so the program aborts with a message if there is already a file named 'comp1.gml' in the directory.

* The node numbers in each 'compx.gml' file are consecutively numbered and start with '1'.  A companion file is output with each 'compx.gml'.  It is named 'compx.nn' and contains the original node numbers for the component. 'compx.nn' files are used for input to the annotation program 'annotate'.

* There is also an option to output a .wg2 file for each component.  Setting BFS_WG2 to '1' in network.h will utilize this option.  .wg2 files can be used as input to Mark Newman's modularity program.  Their format is defined in 'formatSummary.txt'.

___

### Network Library:

* The Network library is defined in ***network.h***.  network.cpp and network.h are object-oriented and are intended to be used by application programs.  

* The data is protected and can only be accessed through the Network object, as shown in network.h.  

* When writing an application program using the Network objects, first create a network, then call the functions via that Network.  See bfsNet.cpp for an example application program that uses these objects.

* Note that this library has been optimized for sparse networks and might not be very efficient for dense networks.


***ERBGA*** : Developed with :heart: by Aditya Karnam. 
