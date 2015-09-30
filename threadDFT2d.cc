// Threaded two-dimensional Discrete FFT transform
// Parav Nagarsheth
// ECE6122 Project 2

#include <iostream>
#include <string>
#include <math.h>
#include <cstring>

#include "Complex.h"
#include "InputImage.h"

// You will likely need global variables indicating how
// many threads there are, and a Complex* that points to the
// 2d image being transformed.
int nRows, nCols, nRowsPerThread;
InputImage *image;
Complex *weights, *invweights;

struct barrier_struct {
  pthread_mutex_t mutex;
  int threads_left;
  int n_threads;
  int g_sense;
};

Complex* W_array=NULL;
struct barrier_struct mybarrier;

using namespace std;

// Function to reverse bits in an unsigned integer
// This assumes there is a global variable N that is the
// number of points in the 1D transform.
unsigned ReverseBits(unsigned v, unsigned N)
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
int MyBarrier_Init(struct barrier_struct *barrier, int n_threads)// you will likely need some parameters)
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
void MyBarrier(struct barrier_struct *barrier, int* l_sense) // Again likely need parameters
{
  *l_sense = !(*l_sense);
  pthread_mutex_lock(&(barrier->mutex));
  barrier->threads_left--;
  pthread_mutex_unlock(&(barrier->mutex));
  if (!barrier->threads_left)
  {
    barrier->threads_left = barrier->n_threads;
    barrier->g_sense = (*l_sense);
  }
  else
    while(barrier->g_sense != *l_sense);
}

/*
 * Quick Transpose of square matrix
 */
void SqMatrixTranspose(Complex *mat, int w, int h)
{
  int i, j;
  Complex *temp = new Complex[w*h];
  memcpy(temp, mat, w*h*sizeof(Complex));
  for (i = 0; i < w; i++) {
      for (j = 0; j < h; j++) {
        mat[j*h+i] = temp[i*w+j];
      }
   }
   delete[] temp;
}

/*
 * This precomputes N/2 weights for DFT of length N
 * 
 */
void precompute_weights(Complex* weights, int N, int sign = -1)
{
  for (int i = 0; i < N/2; i++)
  {
    weights[i].real = cos(2*M_PI*i/N);
    weights[i].imag = sign*sin(2*M_PI*i/N);
  }
  cout<<endl;
}

Complex W(int i, int j, unsigned N)
{
  return(W_array[i*(N/j)]);
}

/*
 * In-place recursive implementation of the Danielson-Lanczos FFT
 */
 
void Transform1D(Complex* h, Complex* weights, int direction, int N)
{
  // Implement the efficient Danielson-Lanczos DFT here.
  // "h" is an input/output parameter
  // "N" is the size of the array (assume even power of 2)
  int i = 1;
  Complex tmp;
  int count;
  Complex *H = new Complex [N];
  for(int f=0; f<N; f++)
  {
    H[f] = h[ReverseBits(f, N)];
  }
  memcpy(h, H, N*sizeof(Complex));
  delete[] H;
  while(i != N)
  {
    i *= 2;
    count = 0;
    for (int k = 0; k < N; k+=i)
    {
      count++;
      for (int j = 0; j < i/2; j++)
      {
        tmp = h[j+k];
        h[j+k] = tmp + weights[j*N/i]*h[j+k+i/2];
        h[j+k+i/2] = tmp - weights[j*N/i]*h[j+k+i/2];
      }
    }
  }
  if (direction == -1)
  {
    for (int f = 0; f < N; f++)
    {
      h[f].real = h[f].real/N;
      h[f].imag = h[f].imag/N;
    }
  }
}

void* Transform2DTHread(void* v)
{ // This is the thread startign point.  "v" is the thread number
  // Calculate 1d DFT for assigned rows
  // wait for all to complete
  // Calculate 1d DFT for assigned columns
  // Decrement active count and signal main if all complete
  Complex* image_ptr = image->GetImageData();
  long threadNumber  = (long)v;
  Complex* thread_image_ptr = image_ptr + nRowsPerThread*threadNumber*nCols;
  int localsense = mybarrier.g_sense;
  printf("Inside thread %lu\n", threadNumber);
  for (int i = 0; i < nRowsPerThread; i++)
  { 
    Transform1D(thread_image_ptr+i*nCols, weights, 1, nCols);
  }
  MyBarrier(&mybarrier, &localsense);
  if (threadNumber == 0)
  {
      image->SaveImageData("MyAfter1d.txt", image_ptr, nCols, nRows);
      SqMatrixTranspose(image_ptr, nCols, nRows);
  }
  MyBarrier(&mybarrier, &localsense);
  for (int i = 0; i < nRowsPerThread; i++)
  { 
    Transform1D(thread_image_ptr+i*nCols, weights, 1, nCols);
  }
  MyBarrier(&mybarrier, &localsense);
  if (threadNumber == 0)
  {
    SqMatrixTranspose(image_ptr, nCols, nRows);
    image->SaveImageData("Tower-DFT2D.txt", image_ptr, nCols, nRows);
    printf("Thread %lu done saving\n", threadNumber);
  }
  MyBarrier(&mybarrier, &localsense);
  printf("Thread %lu going for 1st inverse\n", threadNumber);

  // Inverse transform
  for (int i = 0; i < nRowsPerThread; i++)
  { 
    Transform1D(thread_image_ptr+i*nCols, invweights, -1, nCols);
  }
  MyBarrier(&mybarrier, &localsense);
  if (threadNumber == 0)
  {
      image->SaveImageData("MyAfterInverse1d.txt", image_ptr, nCols, nRows);
      SqMatrixTranspose(image_ptr, nCols, nRows);
      printf("Thread %lu done with saving\n", threadNumber);
  }
  MyBarrier(&mybarrier, &localsense);
  printf("Thread %lu going for 2nd inverse transform\n", threadNumber);
  for (int i = 0; i < nRowsPerThread; i++)
  {
    Transform1D(thread_image_ptr+i*nCols, invweights, -1, nCols);
  }
  MyBarrier(&mybarrier, &localsense);
  if (threadNumber == 0)
  {
    SqMatrixTranspose(image_ptr, nCols, nRows);
    image->SaveImageData("TowerInverse.txt", image_ptr, nCols, nRows);
  }

  printf("Exiting thread %lu\n", threadNumber);
  pthread_exit(NULL);
}

void Transform2D(const char* inputFN, unsigned numThreads) 
{ // Do the 2D transform here.
  image = new InputImage(inputFN);
  // Create the global pointer to the image array data
  // Create 16 threads
  // Wait for all threads complete
  // Write the transformed data
  nRows = image->GetHeight();
  nCols = image->GetWidth();
  nRowsPerThread = nRows/numThreads;
  weights = new Complex[nCols/2];
  invweights = new Complex[nCols/2];

  precompute_weights(weights, nCols, -1);
    // Inverse transform
  precompute_weights(invweights, nCols, 1);

  pthread_t *threads = new pthread_t[numThreads];
  MyBarrier_Init(&mybarrier, numThreads);

  for (unsigned int i = 0; i < numThreads; i++)
    pthread_create(&threads[i], NULL, Transform2DTHread, (void*)i);

  for (unsigned int i = 0; i < numThreads; i++)
    pthread_join(threads[i], NULL);
  delete [] threads;
}

int main(int argc, char** argv)
{
  string fn("Tower.txt"); // default file name
  if (argc > 1) fn = string(argv[1]);  // if name specified on cmd line
  int nThreads = 16;
  // MPI initialization here
  Transform2D(fn.c_str(), nThreads); // Perform the transform.
}  
