// -------------------------------------------------------------------------
// cd_spnet.h -   Header file for breadth-first search of sparse networks
//
// written by Aditya Karnam, July 2017
//
// ------------------------------------------------------------------------

#ifndef _CDSPNET_H
#define _CDSPNET_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include "timer.h"

const int QUIET = 1;  // set to one to eliminate output to screen
const int VERBOSE = 0;  // set to one to display maximum output to screen


inline void warning(const char* p) { fprintf(stderr,"Warning: %s \n",p); }
inline void fatal(const char* string) {fprintf(stderr,"Fatal: %s\n",string);
                                 exit(1); }

#endif
