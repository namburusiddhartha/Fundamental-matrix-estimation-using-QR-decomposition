//#include <algorithm>
//using namespace std;
//
//
#include <math.h>
#include <stdlib.h>
#include <stdio.h>



#define FMADD(dest, src1, src2) \
  __asm__ __volatile__ (      \
  "vfmadd231pd %[rsrc1], %[rsrc2], %[rdest]\n"  \
    : [rdest] "+x"(dest)     \
    : [rsrc1] "x" (src1) , [rsrc2] "x"(src2));
    



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




void naive
(
  int              m,
  int              n,
  double* restrict A,      
  double* restrict F,
  double* restrict Q
){
 
    int i, j, k;
    double temp1, temp2, ulen, vnorm;

    double * a;
    double * ac;
    double * Qc;
    posix_memalign((void**) &a, 64, n * n * sizeof(double));
    posix_memalign((void**) &Qc, 64, n * n * sizeof(double));
    posix_memalign((void**) &ac, 64, n * n * sizeof(double));


    // Find A^T A

    
   for (int r = 0; r < n; r++){
     for (int c = 0; c < n; c++){
       a[c + n*r] = 0;
       for (int k = 0; k < m; k++){
         a[c + n*r] += A[k*n + r] * A[k*n + c];
       }
     }
   }


   m = n;
    
    // Loop across i starts
    for(i = 0; i < n; i ++){
      
      ulen = 0;
      for(j = i; j < m; j ++){
        ulen += a[n*j + i]*a[n*j + i];
      }
      ulen = sqrt(ulen);
      
      // Initial vnorm
      if(a[n*i + i] < 0)
        ulen = ulen;
      else
        ulen = -ulen;

      // Modified vnorm
      vnorm = 0;
      for(j = i; j < m; j ++){

        if(j == i){
          vnorm += (a[n*j + i] + ulen)*(a[n*j + i] + ulen);
        }
        else{
          vnorm += (a[n*j + i])*(a[n*j + i]);
        }
      }
      vnorm = sqrt(vnorm);
      if(vnorm < 1e-18){
        i = i + 1;
        if(i >= n){
          break;
        }
      }
      
      // Copy value for matrix multiplication
      for(j = 0; j < m; j ++){
        for(int k = 0; k < m; k++){
          
          ac[j*m + k] = a[j*m + k];
          
          if(i == 0){
            if(j == k)
              Qc[j*m + k] = 1;
            else
              Qc[j*m + k] = 0;
          }
          else
          Qc[j*m + k] = Q[j*m + k];
        
        }
      }
      
      // Q = Q * vvT (ignoring zero multiplications)
      for(int l = i; l < m; l ++){
        for(int j = 0; j < n; j++){
        
          Q[l*m + j] = 0;
          
        
          for(int k = i; k < m ; k++){
          
            if(l == i)
            temp1 = ac[n*l + i] + ulen;
            else
            temp1 = ac[n*l + i];
            if(k == i)
            temp2 = ac[n*k + i] + ulen;
            else
            temp2 = ac[n*k + i];
          
            if(l == k){
              temp1 = 1 - 2 * (temp1/vnorm) * (temp2/vnorm);}
            else{
              temp1 = 0 - 2 * (temp1/vnorm) * (temp2/vnorm);}
        
            Q[l*m + j] += temp1 * Qc[k*m + j];
          
          }
          
               
        }

      }
      
      
      
      // R = R * vvT (ignoring zero multiplications)
      for(int l = i; l < m; l ++){
        for(int j = i; j < n; j++){
        
          a[l*m + j] = 0;
          
        
          for(int k = i; k < m ; k++){
          
            if(l == i)
            temp1 = ac[n*l + i] + ulen;
            else
            temp1 = ac[n*l + i];
            if(k == i)
            temp2 = ac[n*k + i] + ulen;
            else
            temp2 = ac[n*k + i];
          
            if(l == k){
              temp1 = 1 - 2 * (temp1/vnorm) * (temp2/vnorm);}
            else{
              temp1 = 0 - 2 * (temp1/vnorm) * (temp2/vnorm);}
        
            a[l*m + j] += temp1 * ac[k*m + j];
          
          }
          
               
        }

      }
      
    }
    
    
    
    // Print out the results
    
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
    
    for(int i = 0; i < n; i++){
      F[i] = Q[n*(n-1) + i];
    }
    
    
 
}


void kernel
(
  int              m,
  int              n,
  double* restrict A,      
  double* restrict F,
  double* restrict Q
){

//   __m256d Asimd1, Asimd2, Asimd3, Asimd4, Asimd5, Asimd6, Asimd7, Asimd8, Asimd9, Asimd10, Asimd11, Asimd12, Asimd13, Asimd14, Asimd15, Asimd16;
//   __m256d Csimd1, Csimd2, Csimd3, Csimd4, Csimd5, Csimd6, Csimd7, Csimd8, Csimd9, Csimd10, Csimd11, Csimd12, Csimd13, Csimd14, Csimd15, Csimd16; 
//   
//   
//   Asimd1 = _mm256_load_pd(A + n*0);
//   Asimd2 = _mm256_load_pd(A + n*1);
//   Asimd3 = _mm256_load_pd(A + n*2);
//   Asimd4 = _mm256_load_pd(A + n*3);
//   Asimd5 = _mm256_load_pd(A + n*4);
//   Asimd6 = _mm256_load_pd(A + n*5);
//   Asimd7 = _mm256_load_pd(A + n*6);
//   Asimd8 = _mm256_load_pd(A + n*7);
//   
//   
//   __m256d b1, b2, b3, b4, b5, b6, b7, b8; 
//   
//   for(int i = 0; i < n; i++){
//   
//   broadcastan = _mm256_broadcast_sd(n*0 + 0);
//   FMADD(Csimd1, broadcastan, Asimd1);
//   
//   broadcastan = _mm256_broadcast_sd(n*1 + 0);
//   FMADD(Csimd1, broadcastan, Asimd2);
//   
//   broadcastan = _mm256_broadcast_sd(n*2 + 0);
//   FMADD(Csimd1, broadcastan, Asimd3);
//   
//   broadcastan = _mm256_broadcast_sd(n*3 + 0);
//   FMADD(Csimd1, broadcastan, Asimd4);
//   
//   broadcastan = _mm256_broadcast_sd(n*4 + 0);
//   FMADD(Csimd1, broadcastan, Asimd5);
//   
//   broadcastan = _mm256_broadcast_sd(n*5 + 0);
//   FMADD(Csimd1, broadcastan, Asimd6);
//   
//   broadcastan = _mm256_broadcast_sd(n*6 + 0);
//   FMADD(Csimd1, broadcastan, Asimd7);
//   
//   broadcastan = _mm256_broadcast_sd(n*7 + 0);
//   FMADD(Csimd1, broadcastan, Asimd8);
//   
//   
//   b2 = _mm256_broadcast_sd(a*i + 1);
//   b3 = _mm256_broadcast_sd(a*i + 2);
//   b4 = _mm256_broadcast_sd(a*i + 3);

//   __m256d csimd1, csimd2, csimd3, csimd4, csimd5, csimd6, csimd7, csimd8, csimd9, csimd10, csimd11, csimd;
//   __m256d bsimd1, bsimd2, bsimd3;
//   __m256d broadcastan;
//
//  csimd  = _mm256_load_pd(c);
//  csimd1  = _mm256_load_pd(c + 4);
//  csimd2  = _mm256_load_pd(c + 8);
//  csimd3  = _mm256_load_pd(c + 12);
//  csimd4  = _mm256_load_pd(c + 16);
//  csimd5  = _mm256_load_pd(c + 20);
//  csimd6  = _mm256_load_pd(c + 24);
//  csimd7  = _mm256_load_pd(c + 28);
//  csimd8  = _mm256_load_pd(c + 32);
//  
//for (int x = 0; x < k; x ++){
//     
//     
//     bsimd1 = _mm256_load_pd(a + 9*x);
//     bsimd2 = _mm256_load_pd(a + 9*x + 4);
//     
//     
//     broadcastan = _mm256_broadcast_sd(a + 0 + 4*x);
//     FMADD(csimd, broadcastan, bsimd1);
//     FMADD(csimd1, broadcastan, bsimd2);
//     
//     broadcastan = _mm256_broadcast_sd(a + 1 + 4*x);
//     FMADD(csimd3, broadcastan, bsimd1);
//     FMADD(csimd4, broadcastan, bsimd2);
//     FMADD(csimd5, broadcastan, bsimd3);
//     
//     broadcastan = _mm256_broadcast_sd(a + 2 + 4*x);
//     FMADD(csimd6, broadcastan, bsimd1);
//     FMADD(csimd7, broadcastan, bsimd2);
//     FMADD(csimd8, broadcastan, bsimd3);
//     
//     broadcastan = _mm256_broadcast_sd(a + 3 + 4*x);
//     FMADD(csimd9, broadcastan, bsimd1);
//     FMADD(csimd10, broadcastan, bsimd2);
//     FMADD(csimd11, broadcastan, bsimd3);
//   
//    }
////    
//    _mm256_store_pd(c + 0, csimd);
//    _mm256_store_pd(c + 4, csimd1);
//    _mm256_store_pd(c + 8, csimd2);
//    _mm256_store_pd(c + 12, csimd3);
//    _mm256_store_pd(c + 16, csimd4);
//    _mm256_store_pd(c + 20, csimd5);
//    _mm256_store_pd(c + 24, csimd6);
//    _mm256_store_pd(c + 28, csimd7);
//    _mm256_store_pd(c + 32, csimd8);
//    _mm256_store_pd(c + 36, csimd9);
//    _mm256_store_pd(c + 40, csimd10);
//    _mm256_store_pd(c + 44, csimd11);
   
   
   
   
   
   
   //}

   
   






}
