// -------------------------------------------------------------------------
// erbga_driver.h -   Header file for breadth-first search of sparse networks
//
// written by Sharlee Climer, October 2007
//
//  Modified : added DEBUG for Debugging mode in the BFSNET Driver file.
//  Modified by Aditya Karnam, July,2017 - March, 2018
//
// ------------------------------------------------------------------------

#ifndef _BFSNET_H
#define _BFSNET_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <iostream>
#include <iomanip>
#include <getopt.h>
#include "timer.h"

const int QUIET = 1;  // set to one to eliminate output to screen
const int VERBOSE = 0;  // set to one to display maximum output to screen
const int DIRECTED = 0; // set to one for directed graph
const int DEBUG = 0; // set to one for debug mode

inline void warning(const char* p) { fprintf(stderr,"Warning: %s \n",p); }
inline void fatal(const char* string) {fprintf(stderr,"Fatal: %s\n",string);
                                 exit(1); }
/* Flag set by ‘--verbose’. */
static int verbose_flag;
// popsize, gen, input, output are mandatory params
// rest params are optional, picked from header file if not specified
static struct option long_options[] =
{
    /* These options set a flag. */
    {"verbose", no_argument, &verbose_flag, 1},
    
    /* These options don’t set a flag.
        We distinguish them by their indices. */
    {"popsize", required_argument, 0, 'a'},
    {"gen", required_argument, 0, 'b'},
    {"in", required_argument, 0, 'c'},
    {"out", required_argument, 0, 'd'},
    {"rpop", no_argument, 0, 'e'},
    {"tsize", no_argument, 0, 'f'},
    {"elirate", no_argument, 0, 'g'},
    {"mutrate", no_argument, 0, 'h'},
    {"grepsize", no_argument, 0, 'i'},
    {"greprate", no_argument, 0, 'j'},
    {"crosstype", no_argument, 0, 'k'},
    {"popcheckpt", no_argument, 0, 'l'},
    {"help", no_argument, 0, 'm'},
    {0, 0, 0, 0}
};

#endif
