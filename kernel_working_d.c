#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "immintrin.h"

#define MAX_FREQ 3.4
#define BASE_FREQ 2.4

//timing routine for reading the time stamp counter
static __inline__ unsigned long long rdtsc(void) {
  unsigned hi, lo;
  __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
  return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}

void kernel
(
 int               m,
 int               n,
 int               k,
 double*     restrict a,
 double*     restrict b,
 double*     restrict c
 );

void naive
(
 int               m,
 int               n,
 int               k,
 double*     restrict a,
 double*     restrict b,
 double*     restrict c
 ){

  for (int p = 0; p != k; ++p)
    {
      for (int i = 0; i != m; ++i)
	{
	  for (int j = 0; j !=n; ++j)
	    {
	      c[i*n+j] += a[p*m+i] * b[p*n+j];
	    }
	}

    }

}

int main(){

  double *a;
  double *b;
  double *c, *c_check;

  unsigned long long t0, t1;

  //TODO: Change this according to your calculations for the size of the kernel
  int m = 8;  //m is the number of rows of C
  int n = 9;  //n is the number of columns of C

  /*
    Assume the following
        - A is stored in column major order
        - B is stored in row major order
        - C must be written out in row major order
  */
  
  //double *A;
  
  //posix_memalign((void**) &A, 64, m * n * sizeof(double));
  
  double A[8][9] = 
    {
        {4.51249166e+00, 1.57616189e-01, -2.21716994e+00, -1.09471621e+00, -3.82371892e-02, 5.37878417e-01, -2.03524844e+00, -7.10889077e-02, 1.00000000e+00},
        {4.55803664e+00, -2.34418662e-01, 1.95273973e+00, -1.77873502e+00, 9.14798888e-02, -7.62040068e-01, 2.33417520e+00, -1.20046035e-01, 1.00000000e+00},
        {1.77002133e-01, 6.52226503e-03, -3.41981902e-01, -7.60576940e-02, -2.80261277e-03,  1.46949387e-01, -5.17577485e-01, -1.90719596e-02, 1.00000000e+00},
        {4.64699012e-01, -2.03392431e-01, 5.67007552e-01, -5.09929101e-01, 2.23189025e-01, -6.22195536e-01, 8.19564062e-01, -3.58712032e-01, 1.00000000e+00},
        {5.04420718e-01, -3.66671837e-01,  5.32046419e-01, -9.48466587e-01, 6.89456188e-01, -1.00041143e+00,  9.48076522e-01, -6.89172643e-01, 1.00000000e+00},
        {3.19842188e-02,  4.01381137e-02, -1.76711092e-01, 3.32296863e-02, 4.17010943e-02, -1.83592233e-01, -1.80997233e-01, -2.27139752e-01, 1.00000000e+00},
        {4.07584825e+00, -1.33173690e+00, -1.90251974e+00, -3.07245014e+00, 1.00388802e+00,  1.43415473e+00, -2.14234216e+00,  6.99985852e-01, 1.00000000e+00},
        {4.53285432e-02,  1.66627095e-01, -2.05315656e-01,  1.61924150e-01, 5.95230924e-01, -7.33435504e-01, -2.20774899e-01, -8.11565462e-01, 1.00000000e+00}
   };

  for (int k = 32; k <= 32; k+= 32){
 
    //create memory aligned buffers
    posix_memalign((void**) &a, 64, m * k * sizeof(double));
    posix_memalign((void**) &b, 64, n * k * sizeof(double));
    posix_memalign((void**) &c, 64, m * n * sizeof(double));
    posix_memalign((void**) &c_check, 64, m * n * sizeof(double));

    //initialize A
    for (int i = 0; i != k * m; ++i){
      a[i] = ((double) rand())/ ((double) RAND_MAX);
    }
    //initialize B
    for (int i = 0; i != k * n; ++i){
      b[i] = ((double) rand())/ ((double) RAND_MAX);
    }
    //initialize C
    for (int i = 0; i != m * n; ++i){
      c[i] = 0.0;
      c_check[i] = 0.0;
    }

    naive(m, n, k, a, b, c_check);

    t0 = rdtsc();

    kernel(m, n, k, a, b, c);

    t1 = rdtsc();


    int correct = 1;
    for (int i = 0; i != m * n; ++i) {
      correct &= (fabs(c[i] - c_check[i]) < 1e-13);
    }

    //printf("%d\t %d\t %d\t %lf %d\n", m, n, k, (2.0*m*n*k)/((double)(t1-t0)*MAX_FREQ/BASE_FREQ), correct);

    free(a);
    free(b);
    free(c);
    free(c_check);
  }

  return 0;
}
