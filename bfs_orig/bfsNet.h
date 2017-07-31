// -------------------------------------------------------------------------
// bfsNet.h -   Header file for breadth-first search of sparse networks
//
// written by Sharlee Climer, October 2007
//
// update : added DEBUG for Debugging mode in the BFSNET Driver file.
//
// ------------------------------------------------------------------------

#ifndef _BFSNET_H
#define _BFSNET_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include "timer.h"

const int QUIET = 1;  // set to one to eliminate output to screen
const int VERBOSE = 0;  // set to one to display maximum output to screen
const int DIRECTED = 0; // set to one for directed graph
const int DEBUG = 0; // set to one for debug mode

inline void warning(const char* p) { fprintf(stderr,"Warning: %s \n",p); }
inline void fatal(const char* string) {fprintf(stderr,"Fatal: %s\n",string);
                                 exit(1); }

#endif
