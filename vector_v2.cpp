#include <iostream>
#include <cstdlib>
#include <vector>
#include <ctime>
#include "timer.h"
using namespace std; 

#define VERBOSE 0

// Random number reference : https://www.daniweb.com/programming/software-development/threads/1769/c-random-numbers
 
int get_random_number(int min, int max){
   int random_integer; 
   int lowest=min, highest=max; 
   int range=(highest-lowest)+1; 
   random_integer = lowest+int((long)range*rand()/(RAND_MAX + 1.0)); 
   return random_integer;
}

int main() {
   
   long long int MAX_VECTOR_LENGTH = 5000000;
   long long int MAX_ITERATIONS = MAX_VECTOR_LENGTH;
   long long int MAX_PUSH_POPS = 10000;
   // start timer
   timer t;
   t.start("\nTimer started.");

   // 
   srand((unsigned)time(0));

   // create a vector to store int
   vector<int> vec; 
   long long int i;

   // display the original size of vectors
   cout << "vector size = " << vec.size() << endl;

   // push 5 values into the vector
   for(i = 0; i < MAX_VECTOR_LENGTH; i++){
      vec.push_back(i);
   }

   // display extended size of vector
   cout << "vector size = " << vec.size() << endl;

   // 5000000 (5 million) iterations
   for (int i = 0; i < MAX_ITERATIONS; ++i)
   {
      if (i%2==0)
      {
         // use push operations
         for (int j = 0; j < get_random_number(0, MAX_PUSH_POPS); ++j)
         {
            // v2.0
            // push value = (i+j) into the vector. this value gets pushed at the end of the vector
            // vec.push_back(i+j);

            // v2.1
            // modify the valie to i+j at random locations
            long long int loc = get_random_number(0, MAX_VECTOR_LENGTH - 1);
            //cout << "Cur loc : " << loc << endl;
            vec.at( loc ) = (i+j);
         }
         // end of push

      }else{
         // use pop operations
         for (int j = 0; j < get_random_number(0, MAX_PUSH_POPS); ++j)
         {
            // v2.0
            // vec.pop_back();
            // erase or pop an element at a random position in the vector
            // vec.erase( vec.begin() + get_random_number(0, vec.size() - 1) );

            // v2.1
            // modify value to -1 at random locations in the vector
            long long int loc = get_random_number(0, MAX_VECTOR_LENGTH - 1);
            //cout << "Cur loc : " << loc << endl;
            vec.at( loc ) = -1;

         }
         // end of push
      }
   }

   // use iterator to access the values
   if(VERBOSE){
   vector<int>::iterator v = vec.begin();
      while( v != vec.end()) {
         cout << " value of v = " << *v << endl;
         v++;
      }
   }

   t.stop();
   cout << "Time Elapsed : " << t << " seconds." << endl;

   return 0;
}