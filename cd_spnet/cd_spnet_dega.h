// -------------------------------------------------------------------------
// cd_spnet_dega.h -   Header file for Community Detection Using Differential Evolution
//
// Friday, June 23 2017
// Aditya Karnam
//
// ------------------------------------------------------------------------

#ifndef _CDSPNETDEGA_H
#define _CDSPNETDEGA_H


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>


#include "timer.h"

const int QUIET = 0;  // set to 1 to eliminate most output to screen (Boolean)
const int VERBOSE = 0;  // set to 1 to display maximum output to screen (Boolean)

// the following can be adjusted if needed
const int MAX_STRNG = 1000; // maximum string length
const double TOL = 0.00001; // tolerance


// inline functions
inline void warning(const char* p) { fprintf(stderr,"\nWarning: %s \n",p); }
inline void fatal(const char* string) {fprintf(stderr,"\nFatal: %s\n\n",string); exit(1); }

#endif
