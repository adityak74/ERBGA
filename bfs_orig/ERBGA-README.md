# ERBGA
-----------------

Project Structure:

```
|->
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
	|run_multiple_islands - dependency file for run_ga.sh
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

## 

***ERBGA*** : Developed with :heart: by Aditya Karnam. 