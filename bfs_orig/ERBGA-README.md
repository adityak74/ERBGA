# ERBGA
-----------------

Project Structure:

```
/
	|README
	|README.md
	|calc.cpp - Implementation of the calculations
	|
	|forkos_4.c - Prob#4
	|forkos_5.c - Prob#5
	|forkos_6.c - Prob#6
	|forkos_7.c - Prob#7
	|Makefile
```

Compiling:

```
$ make all
```

Running the project (Prob #5)

```
$ ./forkos -n3 -k4 -m10
```

Help:

```
$ ./forkos -h
Usage : ./forkos -n processes -k loops -m timetosleep
n = number of process to fork
k = number of times to loop on fprintf
m = time to sleep in millis
```

Clean the project:

```
$ make clean
```

Clean the project only *.o:

```
$ make cleanobj
```

***Process Chains*** : Project #1 as a part of CS4760. 