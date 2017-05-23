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
   for(i = 0; i < 5000000; i++){
      vec.push_back(i);
   }

   // display extended size of vector
   cout << "extended vector size = " << vec.size() << endl;

   // 5000000 (5 million) iterations
   for (int i = 0; i < 5000000; ++i)
   {
      if (i%2==0)
      {
         // use push operations
         for (int j = 0; j < get_random_number(0, 10000); ++j)
         {
            // push value = (i+j) into the vector. this value gets pushed at the end of the vector
            vec.push_back(i+j);
         }
         // end of push

      }else{
         // use pop operations
         for (int j = 0; j < get_random_number(0, 10000); ++j)
         {
            //vec.pop_back();
            // erase or pop an element at a random position in the vector
            vec.erase( vec.begin() + get_random_number(0, vec.size() - 1) );
         }
         // end of push
      }
   }

   // use iterator to access the values
   if(VERBOSE){
   vector<int>::iterator v = vec.begin();
      while( v != vec.end()) {
         cout << "value of v = " << *v << endl;
         v++;
      }
   }

   t.stop();
   cout << "Time Elapsed : " << t << " seconds." << endl;

   return 0;
}