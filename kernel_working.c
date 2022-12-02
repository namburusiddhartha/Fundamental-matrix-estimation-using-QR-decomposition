//#include <algorithm>
//using namespace std;
//
//
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "linearalgebra.h"


#define FMADD(dest, src1, src2) \
  __asm__ __volatile__ (      \
  "vfmadd231pd %[rsrc1], %[rsrc2], %[rdest]\n"  \
    : [rdest] "+x"(dest)     \
    : [rsrc1] "x" (src1) , [rsrc2] "x"(src2));
    
#define MATMUL(m, n, k, a, b, c)


//void householder (double **a, double **v, int m, int n) {
//    int i, j;
//    double vnorm, vTa, vpartdot;
//
//    for(i = 0; i < n; i++) {
//        /* set v[i] equal to subvector a[i][i : m] */
//        partialvec_copy(a[i], v[i], m - i, i);
//
//        /* vpartdot = ||v[i]||^2 - v[i][0] * v[i][0]; since vpartdot 
//           is unaffected by the change in v[i][0], storing this value 
//           prevents the need to recalculate the entire norm of v[i] 
//           after updating v[i][0] in the following step              */
//        vpartdot = partialdot_product(v[i], v[i], m - i, 1);
//
//        /* set v[i][0] = v[i][0] + sign(v[i][0]) * ||v[i]|| */
//        if(v[i][0] < 0) {
//            v[i][0] -= sqrt(v[i][0] * v[i][0] + vpartdot);
//        }
//        else {
//            v[i][0] += sqrt(v[i][0] * v[i][0] + vpartdot);
//        }
//
//        /* normalize v[i] */
//        vnorm = sqrt(v[i][0] * v[i][0] + vpartdot);
//        scalar_div(v[i], vnorm, m - i, v[i]);
//    
//        for(j = i; j < n; j++) {
//            /* set a[j][i:m] = a[j][i:m] - 2 * (v[i]^T a[j][i:m]) * v[i] */
//            vTa = subdot_product(a[j], v[i], m - i, i);
//            vTa *= 2;
//            partialscalar_sub(v[i], vTa, m - i, i, a[j]);
//        }
//    }
//}




void kernel
(
  int              qq,
  int              qqq,
  int              qqqq,
  double* restrict qqqqq,
  double* restrict qqqqqq,
  double* restrict qqqqqqq
){
 
 //TODO: Implement a small mxk * kxn matrix multiplication here
 
   //int rn = n % 4;
   //n = (n/4)*4;
   int num_points = 1; // k * 4 * 4;
   
   //printf("%d\n", num_points);
   int size = 0;

   int m = 9, n = 9;
   const int p = (m < n) ? m : n;
   const int q = (m > n) ? m : n;
   double X1[16] = { 335., 1763., 831., 1268., 1310., 941., 300., 928., 587., 571., 604., 493., 385., 536., 839., 345.} ;
 
   double X2[16] = { 189., 1501., 779., 1065., 1054., 831., 288., 822., 732., 323., 609., 367., 248., 505., 1014., 332.} ;
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

   int i, j, k; // m, n;
    double x;
    //const int m = 8, n = 9;
    double * a;
    double * ac;
    double AT[n][m];
    //double * v;
    double * b;
    double * P;
    double * Qc;
    double * Q;
    double u[m], v[m];
    posix_memalign((void**) &a, 64, n * n * sizeof(double));
    posix_memalign((void**) &b, 64, n * n * sizeof(double));
    //posix_memalign((void**) &v, 64, n * n * sizeof(double));
    posix_memalign((void**) &P, 64, n * n * sizeof(double));
    posix_memalign((void**) &Q, 64, n * n * sizeof(double));
    posix_memalign((void**) &Qc, 64, n * n * sizeof(double));
    posix_memalign((void**) &ac, 64, n * n * sizeof(double));

    
   for (int r = 0; r < n; r++){
     for (int c = 0; c < m; c++){
        AT[r][c] = A[c][r];
        }
    }

   for(int r = 0; r < n; r++){
     for(int c = 0; c < n; c++){
       a[c + n*r] = 0;
       }
   }

   double a3[3][3] = 
    {
      {-1, -1, 1},
      {1, 3, 3},
      {-1, -1, 5}	

    };
    
   m = 8;
   for (int r = 0; r < n; r++){
     for (int c = 0; c < n; c++){
       for (int k = 0; k < m; k++){
         a[c + n*r] += AT[r][k] * A[k][c];
       }
     }
   }


   m = 9;

    printf("A = \n");
    double alpha = 0;
    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {

            printf("%9.6g ", a[j + n*i]);
        }
        printf("\n");
    }
    printf("\n");
    
    for(i = 0; i < n; i++) {
        for(j = 0; j < m; j++) {
        
            if(i != j){
              b[j + n*i] = 0;
              P[j + n*i] = 0;
              Q[j + n*i] = 0;}
            
            else{
              b[j + n*i] = 1;
              P[j + n*i] = 1;
              Q[j + n*i] = 1;}
    
           }
        }
    
    
    for(i = 0; i < n; i ++){
      
      for(int k = 0; k < m; k++){
        u[k] = 0;
        v[k] = 0;
      }
      
      for(j = i; j < m; j ++){
        u[j] = a[n*j + i];
      }
      
      double ulen = 0;
      
      for(j = 0; j < m; j ++){
        ulen += u[j]*u[j];
      }
      ulen = sqrt(ulen);
      if(u[i] < 0){
        alpha = ulen;
      }
      else{
        alpha = -ulen;
      }
      
      for(j = 0; j < m; j ++){
        if(j == i){
          v[j] = u[j] + alpha;
        }
        else{
          v[j] = u[j];
        }
        
      }
      double vnorm = 0;
      for(j = 0; j < m; j ++){
        vnorm += v[j]*v[j];
      }
      vnorm = sqrt(vnorm);
      //printf("vnorm %f \n", vnorm);
      if(vnorm < 1e-18){
        i = i + 1;
        if(i >= n){
          break;
        }
      }
      
      
      for(j = 0; j < m; j ++){
        v[j] = v[j]/vnorm;
      }

      for(j = 0; j < m; j ++){
        for(int k = 0; k < m; k++){
          
          ac[j*m + k] = a[j*m + k];
          Qc[j*m + k] = Q[j*m + k];
        
        }
      }
      
      for(j = 0; j < m; j ++){
        for(int k = 0; k < m; k++){
        
          
          P[j*m + k] = b[j*m + k] - 2 * v[j] * v[k];;
        
        }
      }
      

      for(int l = 0; l < m; l ++){
        for(int j = 0; j < n; j++){
        
          a[l*m + j] = 0;
          Q[l*m + j] = 0;
        
          for(int k = 0; k < m ; k++){
        
            a[l*m + j] += P[m*l + k] * ac[k*m + j];
            Q[l*m + j] += P[m*l + k] * Qc[k*m + j];
          
          }
               
        }

      }
      
    }
    
    printf("R = \n");
    for(i = 0; i < m; i++) {
        for(j = 0; j < n; j++) {
            printf("%9.6g ", a[j + n*i]);
        }
        printf("\n");
    }
    printf("\n");
    
    printf("Q = \n");
    for(i = 0; i < m; i++) {
        for(j = 0; j < m; j++) {
            printf("%9.6g ", Q[j + m*i]);
        }
        printf("\n");
    }
    printf("\n");
 
}
