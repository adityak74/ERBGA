/****************************************************************************
*
*	calc.cpp:	Implementation of calculations.
*
*                       Sharlee Climer
*                       January 2003
*
****************************************************************************/

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

float calcAvg(float * values, int size)  // calculate average of values
{
  float avg = 0;
  for (int i = 0; i < size; i++)
    avg += values[i];
  return (avg / (float)size);
}

void sortArray(float *data, int size)
{
  float tempData;

  if (size > 1) {
    for (int i = 0; i < size; i++)   // sort data
      for (int j = size-1; j > i; j--)
	if (data[j] < data[j-1]) {
	  tempData = data[j];
	  data[j] = data[j-1];
	  data[j-1] = tempData;
	}
  }
}

void printStatistics(float *data, int size, char *outputFile)
{
  float avgdata = 0.0;
  float sdata = 0.0;
  float confUpData, confLowData;
  float tempData;
  float median;
 
  if (size > 1) {
    for (int i = 0; i < size; i++)   // sort data
      for (int j = size-1; j > i; j--)
	if (data[j] < data[j-1]) {
	  tempData = data[j];
	  data[j] = data[j-1];
	  data[j-1] = tempData;
	}
    
    if ((size % 2) == 0)  // an even number of data
      median = 0.5 * (data[size/2] + data[size/2 - 1]);
    
    else        // odd number of data
      median = data[size/2];

    for (int i = 0; i < size; i++)     // find mean of data
      avgdata += data[i];
    avgdata /= (float)size;

    for (int i = 0; i < size; i++)  // find 95% confidence interval
      sdata += (data[i] - avgdata) * (data[i] - avgdata);
    sdata /= (float)((size - 1) * size);
    sdata = sqrt(sdata);
  }

  else {
    median = avgdata = data[0];
    sdata = 0;
  }

  confUpData = avgdata + (2.0 * sdata);
  confLowData = avgdata - (2.0 * sdata);
 
  FILE *f1;
  f1 = fopen(outputFile, "a");                    // output the results
  fprintf(f1, "max = %0.5f, min = %.05f \nmedian = %0.5f \naverage = %0.5f \n95confup = %0.5f, 95confdown = %0.5f \n", data[size-1], data[0], median, avgdata,
           confUpData, confLowData);
  fclose(f1);
}
