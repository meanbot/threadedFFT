// Threaded two-dimensional Discrete FFT transform
// Parav Nagarsheth
// ECE6122 Project 2


#include <iostream>
#include <string>
#include <math.h>

#include "Complex.h"
#include "InputImage.h"

// You will likely need global variables indicating how
// many threads there are, and a Complex* that points to the
// 2d image being transformed.

struct barrier_struct {
  pthread_mutex_t mutex;
  int threads_left;
  int n_threads;
  int g_sense;
};

using namespace std;

// Function to reverse bits in an unsigned integer
// This assumes there is a global variable N that is the
// number of points in the 1D transform.
unsigned ReverseBits(unsigned v)
{ //  Provided to students
  unsigned n = N; // Size of array (which is even 2 power k value)
  unsigned r = 0; // Return value
   
  for (--n; n > 0; n >>= 1)
    {
      r <<= 1;        // Shift return value
      r |= (v & 0x1); // Merge in next bit
      v >>= 1;        // Shift reversal value
    }
  return r;
}

// GRAD Students implement the following 2 functions.
// Undergrads can use the built-in barriers in pthreads.

// Call MyBarrier_Init once in main
int MyBarrier_Init(struct barrier_struct *barrier, n_threads)// you will likely need some parameters)
{
  int ret;
  ret = pthread_mutex_init(&(barrier->mutex), NULL);
  if (ret < 0)
    return ret;
  barrier->threads_left = n_threads;
  barrier->n_threads = n_threads;
  barrier->g_sense = 0;
  return 0;
}

// Each thread calls MyBarrier after completing the row-wise DFT
void MyBarrier(struct barrier_struct *barrier, int l_sense) // Again likely need parameters
{
  l_sense = !l_sense;
  pthread_mutex_lock(&(barrier->mutex));
  barrier->threads_left--;
  pthread_mutex_unlock(&(barrier->mutex));  
  if (!barrier->threads_left)
  {
    barrier->threads_left = barrier->n_threads;
    barrier->g_sense = l_sense;
  }
  else
    while(barrier->g_sense != l_sense);
}
                    
void Transform1D(Complex* h, int N)
{
  // Implement the efficient Danielson-Lanczos DFT here.
  // "h" is an input/output parameter
  // "N" is the size of the array (assume even power of 2)
}

void* Transform2DTHread(void* v)
{ // This is the thread startign point.  "v" is the thread number
  // Calculate 1d DFT for assigned rows
  // wait for all to complete
  // Calculate 1d DFT for assigned columns
  // Decrement active count and signal main if all complete
  return 0;
}

void Transform2D(const char* inputFN) 
{ // Do the 2D transform here.
  InputImage image(inputFN);  // Create the helper object for reading the image
  // Create the global pointer to the image array data
  // Create 16 threads
  // Wait for all threads complete
  // Write the transformed data
}

int main(int argc, char** argv)
{
  string fn("Tower.txt"); // default file name
  if (argc > 1) fn = string(argv[1]);  // if name specified on cmd line
  // MPI initialization here
  Transform2D(fn.c_str()); // Perform the transform.
}  
