//#include <algorithm>
//using namespace std;
//
//
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <immintrin.h>



#define FMADD(dest, src1, src2) \
  __asm__ __volatile__ (      \
  "vfmadd231pd %[rsrc1], %[rsrc2], %[rdest]\n"  \
    : [rdest] "+x"(dest)     \
    : [rsrc1] "x" (src1) , [rsrc2] "x"(src2));
    



//void householder (double **a, double **v, int m, int n) {
//    int i, j;
//    double asimd12, vTa, vpartdot;
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
//        asimd12 = sqrt(v[i][0] * v[i][0] + vpartdot);
//        scalar_div(v[i], asimd12, m - i, v[i]);
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
    double temp11, temp21, temp12, temp22, temp13, temp23, temp14, temp24, ulen1, ulen2, ulen3, ulen4, vnorm1, vnorm2, vnorm3, vnorm4;

    double * a;
    double * ac;
    double * Qc;
    posix_memalign((void**) &a, 64, 4 * n * n * sizeof(double));
    posix_memalign((void**) &Qc, 64, 4 * n * n * sizeof(double));
    posix_memalign((void**) &ac, 64, 4 * n * n * sizeof(double));

    
    // Find A^T A
//      for (int j = 0; j < n; j++){
//            
//            for (int k = 0; k < 4; k++){
//            //printf("\n   j  = %d", j);
//            //printf("\n   j = %d", j);      
//            //printf("\n   k = %d \n", k); 
//            for (int i = 0; i < m; i++){
//                    
//            printf("%lf \n",  A[4*n*i + 4*j + k]);
//            
//          }
//          printf("\n");
//        }
//        printf("\n");
//      }
      
      for (int j = 0; j < n; j++){
                
            
          for (int p = 0; p < n; p++){
            
            for (int k = 0; k < 4; k++){
//            printf("\n   j  = %d", j);
//            printf("\n   p = %d", p);      
//            printf("\n   k = %d \n", k); 
//                 
              a[4*n*j + 4*p + k] = 0;
              for (int i = 0; i < m; i++){
                    
            a[4*n*j + 4*p + k] += A[4*n*i + 4*j + k] * A[4*n*i + 4*p + k]; 
            //printf("%lf_ %lf \n", A[4*n*i + 4*j + k], A[4*n*i + 4*p + k]);
            }
          
          }
          //printf("\n");
        }
        //printf("\n");
      }
    int test = 0;
    for (int j = 0; j < n; j++){
          for (int p = 0; p < n; p++){
            if(abs(a[4*n*j + 4*p + 0] - a[4*n*j + 4*p + 1]) > 1e-20 || abs(a[4*n*j + 4*p + 2] - a[4*n*j + 4*p + 3]) > 1e-20 || abs(a[4*n*j + 4*p + 1] - a[4*n*j + 4*p + 2]) > 1e-20){
            printf("j = %d p = %d", j, p);
            test = 1;}
            //for (int k = 0; k < 4; k++){
                  
          //printf("%lf ", a[4*n*j + 4*p + k]);
          //}
          //printf("\n");
        }
        //printf("\n");
      }
    if( test == 1){
    printf("\n\n wrong input \n\n");}
    
//   for (int r = 0; r < n; r++){
//     for (int c = 0; c < n; c++){
//       a[c + n*r] = 0;
//       for (int k = 0; k < m; k++){
//         a[c + n*r] += A[k*n + 4*r] * A[k*n + 4*c];
//       }
//     }
//   }


   m = n;
    
    // Loop across i starts
    for(i = 0; i < n; i ++){
      
      ulen1 = 0;
      ulen2 = 0;
      ulen3 = 0;
      ulen4 = 0;
      
      for(j = i; j < m; j ++){
        
        ulen1 += a[4*n*j + 4*i + 0]*a[4*n*j + 4*i + 0];
        ulen2 += a[4*n*j + 4*i + 1]*a[4*n*j + 4*i + 1];
        ulen3 += a[4*n*j + 4*i + 2]*a[4*n*j + 4*i + 2];
        ulen4 += a[4*n*j + 4*i + 3]*a[4*n*j + 4*i + 3];
      }
      ulen1 = sqrt(ulen1);
      ulen2 = sqrt(ulen2);
      ulen3 = sqrt(ulen3);
      ulen4 = sqrt(ulen4);
      
      // Initial vnorm
      if(a[4*n*i + 4*i + 0] < 0)
        ulen1 = ulen1;
      else
        ulen1 = -ulen1;
        
      if(a[4*n*i + 4*i + 1] < 0)
        ulen2 = ulen2;
      else
        ulen2 = -ulen2;
        
      if(a[4*n*i + 4*i + 2] < 0)
        ulen3 = ulen3;
      else
        ulen3 = -ulen3;
        
      if(a[4*n*i + 4*i + 3] < 0)
        ulen4 = ulen4;
      else
        ulen4 = -ulen4;

      
      // Modified vnorm
      vnorm1 = 0;
      vnorm2 = 0;
      vnorm3 = 0;
      vnorm4 = 0;
      
      for(j = i; j < m; j ++){

        if(j == i){
          vnorm1 += (a[4*n*j + 4*i + 0] + ulen1)*(a[4*n*j + 4*i+ 0] + ulen1);
          vnorm2 += (a[4*n*j + 4*i + 1] + ulen2)*(a[4*n*j + 4*i+ 1] + ulen2);
          vnorm3 += (a[4*n*j + 4*i + 2] + ulen3)*(a[4*n*j + 4*i+ 2] + ulen3);
          vnorm4 += (a[4*n*j + 4*i + 3] + ulen4)*(a[4*n*j + 4*i+ 3] + ulen4);
        
        }
        else{
          vnorm1 += (a[4*n*j + 4*i + 0])*(a[4*n*j + 4*i + 0]);
          vnorm2 += (a[4*n*j + 4*i + 1])*(a[4*n*j + 4*i + 1]);
          vnorm3 += (a[4*n*j + 4*i + 2])*(a[4*n*j + 4*i + 2]);
          vnorm4 += (a[4*n*j + 4*i + 3])*(a[4*n*j + 4*i + 3]);
        }
      }
      
      
      
      vnorm1 = sqrt(vnorm1);
      vnorm2 = sqrt(vnorm2);
      vnorm3 = sqrt(vnorm3);
      vnorm4 = sqrt(vnorm4);
      
      
      
      
      
      
      
      if(vnorm1 < 1e-18){
        i = i + 1;
        if(i >= n){
          break;
        }
      }
      
      // Copy value for matrix multiplication
      for(j = 0; j < m; j ++){
        for(int k = 0; k < m; k++){
          
          ac[4*j*m + 4*k + 0] = a[4*j*m + 4*k + 0];
          ac[4*j*m + 4*k + 1] = a[4*j*m + 4*k + 1];
          ac[4*j*m + 4*k + 2] = a[4*j*m + 4*k + 2];
          ac[4*j*m + 4*k + 3] = a[4*j*m + 4*k + 3];
          
          if(i == 0){
            if(j == k){
              Qc[4*j*m + 4*k + 0] = 1;
              Qc[4*j*m + 4*k + 1] = 1;
              Qc[4*j*m + 4*k + 2] = 1;
              Qc[4*j*m + 4*k + 3] = 1;}
            else{
              Qc[4*j*m + 4*k + 0] = 0;
              Qc[4*j*m + 4*k + 1] = 0;
              Qc[4*j*m + 4*k + 2] = 0;
              Qc[4*j*m + 4*k + 3] = 0;}
          }
          else{
          Qc[4*j*m + 4*k + 0] = Q[4*j*m + 4*k + 0];
          Qc[4*j*m + 4*k + 1] = Q[4*j*m + 4*k + 1];
          Qc[4*j*m + 4*k + 2] = Q[4*j*m + 4*k + 2];
          Qc[4*j*m + 4*k + 3] = Q[4*j*m + 4*k + 3];}
        
        }
      }
      
      // Q = Q * vvT (ignoring zero multiplications)
      for(int l = i; l < m; l ++){
        for(int j = 0; j < n; j++){
        
          Q[4*l*m + 4*j + 0] = 0;
          Q[4*l*m + 4*j + 1] = 0;
          Q[4*l*m + 4*j + 2] = 0;
          Q[4*l*m + 4*j + 3] = 0;
        
          for(int k = i; k < m ; k++){
          
            if(l == i){
            temp11 = ac[4*n*l + 4*i + 0] + ulen1;
            temp12 = ac[4*n*l + 4*i + 1] + ulen2;
            temp13 = ac[4*n*l + 4*i + 2] + ulen3;
            temp14 = ac[4*n*l + 4*i + 3] + ulen4;}
            
            else{
            temp11 = ac[4*n*l + 4*i + 0];
            temp12 = ac[4*n*l + 4*i + 1];
            temp13 = ac[4*n*l + 4*i + 2];
            temp14 = ac[4*n*l + 4*i + 3];}
            
            
            if(k == i){
            temp21 = ac[4*n*k + 4*i + 0] + ulen1; 
            temp22 = ac[4*n*k + 4*i + 1] + ulen2;
            temp23 = ac[4*n*k + 4*i + 2] + ulen3;
            temp24 = ac[4*n*k + 4*i + 3] + ulen4;}
            
            else{
            temp21 = ac[4*n*k + 4*i + 0];
            temp22 = ac[4*n*k + 4*i + 1];
            temp23 = ac[4*n*k + 4*i + 2];
            temp24 = ac[4*n*k + 4*i + 3];}
          
            if(l == k){
              temp11 = 1 - 2 * (temp11/vnorm1) * (temp21/vnorm1);
              temp12 = 1 - 2 * (temp12/vnorm2) * (temp22/vnorm2);
              temp13 = 1 - 2 * (temp13/vnorm3) * (temp23/vnorm3);
              temp14 = 1 - 2 * (temp14/vnorm4) * (temp24/vnorm4);}
                      
            else{
              temp11 = 0 - 2 * (temp11/vnorm1) * (temp21/vnorm1);
              temp12 = 0 - 2 * (temp12/vnorm2) * (temp22/vnorm2);
              temp13 = 0 - 2 * (temp13/vnorm3) * (temp23/vnorm3);
              temp14 = 0 - 2 * (temp14/vnorm4) * (temp24/vnorm4);}
            
            Q[4*l*m + 4*j + 0] += temp11 * Qc[4*k*m + 4*j + 0];
            Q[4*l*m + 4*j + 1] += temp12 * Qc[4*k*m + 4*j + 1];
            Q[4*l*m + 4*j + 2] += temp13 * Qc[4*k*m + 4*j + 2];
            Q[4*l*m + 4*j + 3] += temp14 * Qc[4*k*m + 4*j + 3];          
          }
          
               
        }

      }
      
      // print
//      if(i == 0){
//      printf("\n For i = %d \n Q = \n", i);
//      for(int ind = 0; ind < m; ind++) {
//        for(j = 0; j < m; j++) {
//            printf("%9.6g ", Q[4*j + 4*n*ind + 1]);
//        }
//        printf("\n");
//      }
//      printf("\n");
//      }
      
      
      // R = R * vvT (ignoring zero multiplications)
      for(int l = i; l < m; l ++){
        for(int j = i; j < n; j++){
        
          a[4*l*m + 4*j + 0] = 0;
          a[4*l*m + 4*j + 1] = 0;
          a[4*l*m + 4*j + 2] = 0;
          a[4*l*m + 4*j + 3] = 0;
          
        
          for(int k = i; k < m ; k++){
          
            if(l == i){
            temp11 = ac[4*n*l + 4*i + 0] + ulen1;
            temp12 = ac[4*n*l + 4*i + 1] + ulen2;
            temp13 = ac[4*n*l + 4*i + 2] + ulen3;
            temp14 = ac[4*n*l + 4*i + 3] + ulen4;}
            
            else{
            temp11 = ac[4*n*l + 4*i + 0];
            temp12 = ac[4*n*l + 4*i + 1];
            temp13 = ac[4*n*l + 4*i + 2];
            temp14 = ac[4*n*l + 4*i + 3];}
            
            if(k == i){
            temp21 = ac[4*n*k + 4*i + 0] + ulen1;
            temp22 = ac[4*n*k + 4*i + 1] + ulen2;
            temp23 = ac[4*n*k + 4*i + 2] + ulen3;
            temp24 = ac[4*n*k + 4*i + 3] + ulen4;}
            
            else{
            temp21 = ac[4*n*k + 4*i + 0];
            temp22 = ac[4*n*k + 4*i + 1];
            temp23 = ac[4*n*k + 4*i + 2];
            temp24 = ac[4*n*k + 4*i + 3];}
                      
            if(l == k){
              temp11 = 1 - 2 * (temp11/vnorm1) * (temp21/vnorm1);
              temp12 = 1 - 2 * (temp12/vnorm2) * (temp22/vnorm2);
              temp13 = 1 - 2 * (temp13/vnorm3) * (temp23/vnorm3);
              temp14 = 1 - 2 * (temp14/vnorm4) * (temp24/vnorm4);}
                        
            else{
            
              temp11 = 0 - 2 * (temp11/vnorm1) * (temp21/vnorm1);
              temp12 = 0 - 2 * (temp12/vnorm2) * (temp22/vnorm2);
              temp13 = 0 - 2 * (temp13/vnorm3) * (temp23/vnorm3);
              temp14 = 0 - 2 * (temp14/vnorm4) * (temp24/vnorm4);}
        
            a[4*l*m + 4*j + 0] += temp11 * ac[4*k*m + 4*j + 0];
            a[4*l*m + 4*j + 1] += temp12 * ac[4*k*m + 4*j + 1];
            a[4*l*m + 4*j + 2] += temp13 * ac[4*k*m + 4*j + 2];
            a[4*l*m + 4*j + 3] += temp14 * ac[4*k*m + 4*j + 3];
                    
          }
          
               
        }

      }
      
    }
    
    
    
    // Print out the results
    
    printf("R = \n");
    for(i = 0; i < m; i++) {
        for(j = 0; j < n; j++) {
            printf("%9.6g ", a[4*j + 4*n*i + 1]);
        }
        printf("\n");
    }
    printf("\n");
    
    printf("Q = \n");
    for(i = 0; i < m; i++) {
        for(j = 0; j < m; j++) {
            printf("%9.6g ", Q[4*j + 4*n*i + 3]);
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
  double* restrict Q,
  double* restrict Qc,
  double* restrict a,
  double* restrict ac
  
){
    
   
   

    __m256d asimd1, asimd2, asimd3, asimd4, asimd5, asimd6, asimd7, asimd8, asimd9, asimd10, asimd11, asimd12, asimd13, asimd14, asimd15, asimd16;
    //__m256d Asimd1, Asimd2, Asimd3;
    //__m256d asimd11, asimd12;

//    double * a;
//    double * ac;
//    double * Qc;
//    posix_memalign((void**) &a, 64, 4 * n * n * sizeof(double));
//    posix_memalign((void**) &Qc, 64, 4 * n * n * sizeof(double));
//    posix_memalign((void**) &ac, 64, 4 * n * n * sizeof(double));

    
    // Find A^T A

      for (int j = 0; j < n; j++){
          for (int p = 0; p < n; p++){
            
            asimd1  = _mm256_load_pd(a + 4*n*j + 4*p + 0);
            asimd2 = _mm256_load_pd(A + 4*n*0 + 4*j + 0);
            asimd3 = _mm256_load_pd(A + 4*n*0 + 4*p + 0);
            asimd4 = _mm256_load_pd(A + 4*n*1 + 4*j + 0);
            asimd5 = _mm256_load_pd(A + 4*n*1 + 4*p + 0);
            asimd6 = _mm256_load_pd(A + 4*n*2 + 4*j + 0);
            asimd7 = _mm256_load_pd(A + 4*n*2 + 4*p + 0);
            FMADD(asimd1, asimd2, asimd3);
            FMADD(asimd1, asimd4, asimd5);
            FMADD(asimd1, asimd6, asimd7);
            
            asimd8 = _mm256_load_pd(A + 4*n*3 + 4*j + 0);
            asimd9 = _mm256_load_pd(A + 4*n*3 + 4*p + 0);
            asimd10 = _mm256_load_pd(A + 4*n*4 + 4*j + 0);
            asimd11 = _mm256_load_pd(A + 4*n*4 + 4*p + 0);
            asimd12 = _mm256_load_pd(A + 4*n*5 + 4*j + 0);
            asimd13 = _mm256_load_pd(A + 4*n*5 + 4*p + 0);
            asimd14 = _mm256_load_pd(A + 4*n*6 + 4*j + 0);
            asimd15 = _mm256_load_pd(A + 4*n*6 + 4*p + 0);

            FMADD(asimd1, asimd8, asimd9);
            FMADD(asimd1, asimd10, asimd11);
            FMADD(asimd1, asimd12, asimd13);
            
            asimd2 = _mm256_load_pd(A + 4*n*7 + 4*j + 0);
            asimd3 = _mm256_load_pd(A + 4*n*7 + 4*p + 0);
            FMADD(asimd1, asimd14, asimd15);
            FMADD(asimd1, asimd2, asimd3);

            
            _mm256_store_pd(a + 4*n*j + 4*p + 0, asimd1);
            

        }
      }
      
      //matrix_multiply_kernel()
//    int test = 0;
//    for (int j = 0; j < n; j++){
//          for (int p = 0; p < n; p++){
//            if(abs(a[4*n*j + 4*p + 0] - a[4*n*j + 4*p + 1]) > 1e-20 || abs(a[4*n*j + 4*p + 2] - a[4*n*j + 4*p + 3]) > 1e-20 || abs(a[4*n*j + 4*p + 1] - a[4*n*j + 4*p + 2]) > 1e-20){
//            printf("j = %d p = %d", j, p);
//            test = 1;}
//        }
//      }
//    if( test == 1){
//    printf("\n\n wrong input \n\n");}
//    

   m = n;
    
    // Loop across i starts
    for(int i = 0; i < n; i ++){
      printf("iteration i, %d\n", i);
      
      asimd11 = _mm256_setzero_pd();
      
      //printf("asimd11 zero simd %lf, %lf, %lf, %lf\n", asimd11[0], asimd11[1], asimd11[2], asimd11[3]);
      
      //asimd112 = _mm256_setzero_pd();
      //asimd113 = _mm256_setzero_pd();
      //asimd114 = _mm256_setzero_pd();
      
      asimd1  = _mm256_load_pd(a + 4*n*0 + 4*i + 0);
      asimd2  = _mm256_load_pd(a + 4*n*1 + 4*i + 0);
      asimd3  = _mm256_load_pd(a + 4*n*2 + 4*i + 0);
      asimd4  = _mm256_load_pd(a + 4*n*3 + 4*i + 0);
      
      if(i<1)
      FMADD(asimd11, asimd1, asimd1);
      if(i<2)
      FMADD(asimd11, asimd2, asimd2);
      if(i<3)
      FMADD(asimd11, asimd3, asimd3);
      if(i<4)
      FMADD(asimd11, asimd4, asimd4);
      
      asimd5  = _mm256_load_pd(a + 4*n*4 + 4*i + 0);
      asimd6  = _mm256_load_pd(a + 4*n*5 + 4*i + 0);
      asimd7  = _mm256_load_pd(a + 4*n*6 + 4*i + 0);
      asimd8  = _mm256_load_pd(a + 4*n*7 + 4*i + 0);
      asimd9  = _mm256_load_pd(a + 4*n*8 + 4*i + 0);
      
      if(i<5)
      FMADD(asimd11, asimd5, asimd5);
      if(i<6)
      FMADD(asimd11, asimd6, asimd6);
      if(i<7)
      FMADD(asimd11, asimd7, asimd7);
      if(i<8)
      FMADD(asimd11, asimd8, asimd8);
      if(i<9)
      FMADD(asimd11, asimd9, asimd9);


      asimd11 = _mm256_sqrt_pd(asimd11);
      
      
      
      
      asimd10 = _mm256_load_pd(a + 4*n*i + 4*i + 0);

      
      asimd12 = _mm256_mul_pd ( _mm256_and_pd(asimd10, _mm256_set1_pd(-0.0)),  _mm256_set1_pd(-1.0));
      asimd11 = _mm256_xor_pd(asimd11, asimd12);
      


      asimd12 = _mm256_setzero_pd();
      
      if(i == 0){
      asimd1 = _mm256_add_pd(asimd1, asimd11);
      FMADD(asimd12, asimd1, asimd1);
      FMADD(asimd12, asimd2, asimd2);
      FMADD(asimd12, asimd3, asimd3);
      FMADD(asimd12, asimd4, asimd4);     
      FMADD(asimd12, asimd5, asimd5);
      FMADD(asimd12, asimd6, asimd6);
      FMADD(asimd12, asimd7, asimd7);
      FMADD(asimd12, asimd8, asimd8);
      FMADD(asimd12, asimd9, asimd9);}
      else if(i == 1){
      asimd2 = _mm256_add_pd(asimd2, asimd11);
      FMADD(asimd12, asimd2, asimd2);
      FMADD(asimd12, asimd3, asimd3);
      FMADD(asimd12, asimd4, asimd4);     
      FMADD(asimd12, asimd5, asimd5);
      FMADD(asimd12, asimd6, asimd6);
      FMADD(asimd12, asimd7, asimd7);
      FMADD(asimd12, asimd8, asimd8);
      FMADD(asimd12, asimd9, asimd9);}
      else if(i == 2){
      asimd3 = _mm256_add_pd(asimd3, asimd11);
      FMADD(asimd12, asimd3, asimd3);
      FMADD(asimd12, asimd4, asimd4);     
      FMADD(asimd12, asimd5, asimd5);
      FMADD(asimd12, asimd6, asimd6);
      FMADD(asimd12, asimd7, asimd7);
      FMADD(asimd12, asimd8, asimd8);
      FMADD(asimd12, asimd9, asimd9);}
      else if(i == 3){
      asimd4 = _mm256_add_pd(asimd4, asimd11);
      FMADD(asimd12, asimd4, asimd4);     
      FMADD(asimd12, asimd5, asimd5);
      FMADD(asimd12, asimd6, asimd6);
      FMADD(asimd12, asimd7, asimd7);
      FMADD(asimd12, asimd8, asimd8);
      FMADD(asimd12, asimd9, asimd9);}
      else if(i == 4){
      asimd5 = _mm256_add_pd(asimd5, asimd11);    
      FMADD(asimd12, asimd5, asimd5);
      FMADD(asimd12, asimd6, asimd6);
      FMADD(asimd12, asimd7, asimd7);
      FMADD(asimd12, asimd8, asimd8);
      FMADD(asimd12, asimd9, asimd9);}
      else if(i == 5){
      asimd6 = _mm256_add_pd(asimd6, asimd11);
      FMADD(asimd12, asimd6, asimd6);
      FMADD(asimd12, asimd7, asimd7);
      FMADD(asimd12, asimd8, asimd8);
      FMADD(asimd12, asimd9, asimd9);}
      else if(i == 6){
      asimd7 = _mm256_add_pd(asimd7, asimd11);
      FMADD(asimd12, asimd7, asimd7);
      FMADD(asimd12, asimd8, asimd8);
      FMADD(asimd12, asimd9, asimd9);}
      else if(i == 7){
      asimd8 = _mm256_add_pd(asimd8, asimd11);
      FMADD(asimd12, asimd8, asimd8);
      FMADD(asimd12, asimd9, asimd9);}
      else if(i == 8){
      asimd9 = _mm256_add_pd(asimd9, asimd11);
      FMADD(asimd12, asimd9, asimd9);}

      
      asimd12 = _mm256_add_pd(asimd12, _mm256_set1_pd(1.0*1e-24));
      
      
      // Copy value for matrix multiplication
      for(int j = 0; j < m; j ++){
        for(int k = 0; k < m; k++){
          
          ac[4*j*m + 4*k + 0] = a[4*j*m + 4*k + 0];
          ac[4*j*m + 4*k + 1] = a[4*j*m + 4*k + 1];
          ac[4*j*m + 4*k + 2] = a[4*j*m + 4*k + 2];
          ac[4*j*m + 4*k + 3] = a[4*j*m + 4*k + 3];
          
          Qc[4*j*m + 4*k + 0] = Q[4*j*m + 4*k + 0];
          Qc[4*j*m + 4*k + 1] = Q[4*j*m + 4*k + 1];
          Qc[4*j*m + 4*k + 2] = Q[4*j*m + 4*k + 2];
          Qc[4*j*m + 4*k + 3] = Q[4*j*m + 4*k + 3];
        
        }
      }
      
          // Q = Q * vvT (ignoring zero multiplications)
          for(int l = i; l < m; l ++){
            for(int j = 0; j < n; j++){

              
              asimd16 = _mm256_load_pd(a + 4*n*l + 4*i + 0);
              if(l == i){
                asimd16 = _mm256_add_pd(asimd16, asimd11);
              }
              
              
            
              if(i == 0){
              
              asimd15 = _mm256_load_pd(Qc + 4*0*m + 4*j + 0);    
              asimd13 = _mm256_mul_pd(asimd16, asimd1);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 0)), asimd13);
              asimd15 = _mm256_mul_pd(asimd13, asimd15);
              
              
              asimd14 = _mm256_load_pd(Qc + 4*1*m + 4*j + 0); 
              asimd13 = _mm256_mul_pd(asimd16, asimd2);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 1)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              
              
              asimd14 = _mm256_load_pd(Qc + 4*2*m + 4*j + 0);   
              asimd13 = _mm256_mul_pd(asimd16, asimd3);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 2)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(Qc + 4*3*m + 4*j + 0);   
              asimd13 = _mm256_mul_pd(asimd16, asimd4);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 3)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(Qc + 4*4*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd5);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 4)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(Qc + 4*5*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd6);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 5)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(Qc + 4*6*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd7);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 6)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(Qc + 4*7*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd8);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 7)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(Qc + 4*8*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd9);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 8)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              }
              
              else if(i == 1){
              asimd15 = _mm256_load_pd(Qc + 4*1*m + 4*j + 0);  
                   
              asimd13 = _mm256_mul_pd(asimd16, asimd2);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 1)), asimd13);
              
              asimd15 = _mm256_mul_pd(asimd13, asimd15);

              asimd14 = _mm256_load_pd(Qc + 4*2*m + 4*j + 0);   
              asimd13 = _mm256_mul_pd(asimd16, asimd3);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 2)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(Qc + 4*3*m + 4*j + 0);   
              asimd13 = _mm256_mul_pd(asimd16, asimd4);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 3)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(Qc + 4*4*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd5);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 4)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(Qc + 4*5*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd6);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 5)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(Qc + 4*6*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd7);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 6)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(Qc + 4*7*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd8);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 7)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(Qc + 4*8*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd9);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 8)), asimd13);
              FMADD(asimd15, asimd13, asimd14);}
              
              
              else if(i == 2){
              asimd15 = _mm256_load_pd(Qc + 4*2*m + 4*j + 0);    
              asimd13 = _mm256_mul_pd(asimd16, asimd3);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 2)), asimd13);
              asimd15 = _mm256_mul_pd(asimd13, asimd15);
              
              asimd14 = _mm256_load_pd(Qc + 4*3*m + 4*j + 0);   
              asimd13 = _mm256_mul_pd(asimd16, asimd4);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 3)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(Qc + 4*4*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd5);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 4)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(Qc + 4*5*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd6);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 5)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(Qc + 4*6*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd7);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 6)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(Qc + 4*7*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd8);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 7)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(Qc + 4*8*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd9);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 8)), asimd13);
              FMADD(asimd15, asimd13, asimd14);}
              
              
              else if(i == 3){
              asimd15 = _mm256_load_pd(Qc + 4*3*m + 4*j + 0);    
              asimd13 = _mm256_mul_pd(asimd16, asimd4);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 3)), asimd13);
              asimd15 = _mm256_mul_pd(asimd13, asimd15);
              
              asimd14 = _mm256_load_pd(Qc + 4*4*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd5);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 4)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(Qc + 4*5*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd6);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 5)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(Qc + 4*6*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd7);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 6)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(Qc + 4*7*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd8);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 7)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(Qc + 4*8*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd9);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 8)), asimd13);
              FMADD(asimd15, asimd13, asimd14);}
              
              
              else if(i == 4){
              asimd15 = _mm256_load_pd(Qc + 4*4*m + 4*j + 0);    
              asimd13 = _mm256_mul_pd(asimd16, asimd5);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 4)), asimd13);
              asimd15 = _mm256_mul_pd(asimd13, asimd15);
              
              asimd14 = _mm256_load_pd(Qc + 4*5*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd6);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 5)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(Qc + 4*6*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd7);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 6)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(Qc + 4*7*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd8);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 7)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(Qc + 4*8*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd9);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 8)), asimd13);
              FMADD(asimd15, asimd13, asimd14);}
              
              else if(i == 5){
              asimd15 = _mm256_load_pd(Qc + 4*5*m + 4*j + 0);                
              asimd13 = _mm256_mul_pd(asimd16, asimd6);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 5)), asimd13);
              asimd15 = _mm256_mul_pd(asimd13, asimd15);
              
              asimd14 = _mm256_load_pd(Qc + 4*6*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd7);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 6)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(Qc + 4*7*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd8);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 7)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(Qc + 4*8*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd9);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 8)), asimd13);
              FMADD(asimd15, asimd13, asimd14);}
              
              else if(i == 6){
              asimd15 = _mm256_load_pd(Qc + 4*6*m + 4*j + 0);    
              asimd13 = _mm256_mul_pd(asimd16, asimd7);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 6)), asimd13);
              asimd15 = _mm256_mul_pd(asimd13, asimd15);
              
              asimd14 = _mm256_load_pd(Qc + 4*7*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd8);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 7)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(Qc + 4*8*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd9);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 8)), asimd13);
              FMADD(asimd15, asimd13, asimd14);}
              
              else if(i == 7){
              asimd15 = _mm256_load_pd(Qc + 4*7*m + 4*j + 0);    
              asimd13 = _mm256_mul_pd(asimd16, asimd8);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 7)), asimd13);
              asimd15 = _mm256_mul_pd(asimd13, asimd15);
              
              asimd14 = _mm256_load_pd(Qc + 4*8*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd9);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 8)), asimd13);
              FMADD(asimd15, asimd13, asimd14);}
              
              else if(i == 8){
              asimd15 = _mm256_load_pd(Qc + 4*8*m + 4*j + 0);    
              asimd13 = _mm256_mul_pd(asimd16, asimd9);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 8)), asimd13);
              asimd15 = _mm256_mul_pd(asimd13, asimd15);}
              
          
          _mm256_store_pd(Q + 4*l*m + 4*j + 0, asimd15);    
        }

      }
      
    
    for(int l = i; l < m; l ++){
            for(int j = i; j < n; j++){

              
              asimd16 = _mm256_load_pd(ac + 4*n*l + 4*i + 0);
              if(l == i){
                asimd16 = _mm256_add_pd(asimd16, asimd11);
              }
              
              
            
              if(i == 0){
              
              asimd15 = _mm256_load_pd(ac + 4*0*m + 4*j + 0);    
              asimd13 = _mm256_mul_pd(asimd16, asimd1);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 0)), asimd13);
              asimd15 = _mm256_mul_pd(asimd13, asimd15);
              
              
              asimd14 = _mm256_load_pd(ac + 4*1*m + 4*j + 0); 
              asimd13 = _mm256_mul_pd(asimd16, asimd2);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 1)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              
              
              asimd14 = _mm256_load_pd(ac + 4*2*m + 4*j + 0);   
              asimd13 = _mm256_mul_pd(asimd16, asimd3);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 2)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + 4*3*m + 4*j + 0);   
              asimd13 = _mm256_mul_pd(asimd16, asimd4);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 3)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + 4*4*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd5);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 4)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + 4*5*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd6);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 5)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + 4*6*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd7);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 6)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + 4*7*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd8);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 7)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + 4*8*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd9);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 8)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              }
              
              else if(i == 1){
              asimd15 = _mm256_load_pd(ac + 4*1*m + 4*j + 0);  
                   
              asimd13 = _mm256_mul_pd(asimd16, asimd2);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 1)), asimd13);
              
              asimd15 = _mm256_mul_pd(asimd13, asimd15);

              asimd14 = _mm256_load_pd(ac + 4*2*m + 4*j + 0);   
              asimd13 = _mm256_mul_pd(asimd16, asimd3);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 2)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + 4*3*m + 4*j + 0);   
              asimd13 = _mm256_mul_pd(asimd16, asimd4);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 3)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + 4*4*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd5);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 4)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + 4*5*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd6);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 5)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + 4*6*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd7);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 6)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + 4*7*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd8);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 7)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + 4*8*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd9);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 8)), asimd13);
              FMADD(asimd15, asimd13, asimd14);}
              
              
              else if(i == 2){
              asimd15 = _mm256_load_pd(ac + 4*2*m + 4*j + 0);    
              asimd13 = _mm256_mul_pd(asimd16, asimd3);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 2)), asimd13);
              asimd15 = _mm256_mul_pd(asimd13, asimd15);
              
              asimd14 = _mm256_load_pd(ac + 4*3*m + 4*j + 0);   
              asimd13 = _mm256_mul_pd(asimd16, asimd4);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 3)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + 4*4*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd5);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 4)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + 4*5*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd6);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 5)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + 4*6*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd7);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 6)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + 4*7*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd8);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 7)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + 4*8*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd9);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 8)), asimd13);
              FMADD(asimd15, asimd13, asimd14);}
              
              
              else if(i == 3){
              asimd15 = _mm256_load_pd(ac + 4*3*m + 4*j + 0);    
              asimd13 = _mm256_mul_pd(asimd16, asimd4);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 3)), asimd13);
              asimd15 = _mm256_mul_pd(asimd13, asimd15);
              
              asimd14 = _mm256_load_pd(ac + 4*4*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd5);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 4)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + 4*5*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd6);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 5)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + 4*6*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd7);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 6)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + 4*7*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd8);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 7)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + 4*8*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd9);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 8)), asimd13);
              FMADD(asimd15, asimd13, asimd14);}
              
              
              else if(i == 4){
              asimd15 = _mm256_load_pd(ac + 4*4*m + 4*j + 0);    
              asimd13 = _mm256_mul_pd(asimd16, asimd5);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 4)), asimd13);
              asimd15 = _mm256_mul_pd(asimd13, asimd15);
              
              asimd14 = _mm256_load_pd(ac + 4*5*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd6);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 5)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + 4*6*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd7);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 6)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + 4*7*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd8);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 7)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + 4*8*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd9);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 8)), asimd13);
              FMADD(asimd15, asimd13, asimd14);}
              
              else if(i == 5){
              asimd15 = _mm256_load_pd(ac + 4*5*m + 4*j + 0);                
              asimd13 = _mm256_mul_pd(asimd16, asimd6);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 5)), asimd13);
              asimd15 = _mm256_mul_pd(asimd13, asimd15);
              
              asimd14 = _mm256_load_pd(ac + 4*6*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd7);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 6)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + 4*7*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd8);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 7)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + 4*8*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd9);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 8)), asimd13);
              FMADD(asimd15, asimd13, asimd14);}
              
              else if(i == 6){
              asimd15 = _mm256_load_pd(ac + 4*6*m + 4*j + 0);    
              asimd13 = _mm256_mul_pd(asimd16, asimd7);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 6)), asimd13);
              asimd15 = _mm256_mul_pd(asimd13, asimd15);
              
              asimd14 = _mm256_load_pd(ac + 4*7*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd8);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 7)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + 4*8*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd9);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 8)), asimd13);
              FMADD(asimd15, asimd13, asimd14);}
              
              else if(i == 7){
              asimd15 = _mm256_load_pd(ac + 4*7*m + 4*j + 0);    
              asimd13 = _mm256_mul_pd(asimd16, asimd8);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 7)), asimd13);
              asimd15 = _mm256_mul_pd(asimd13, asimd15);
              
              asimd14 = _mm256_load_pd(ac + 4*8*m + 4*j + 0);  
              asimd13 = _mm256_mul_pd(asimd16, asimd9);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 8)), asimd13);
              FMADD(asimd15, asimd13, asimd14);}
              
              else if(i == 8){
              asimd15 = _mm256_load_pd(ac + 4*8*m + 4*j + 0);    
              asimd13 = _mm256_mul_pd(asimd16, asimd9);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 8)), asimd13);
              asimd15 = _mm256_mul_pd(asimd13, asimd15);}
              
          
          _mm256_store_pd(a + 4*l*m + 4*j + 0, asimd15);    
        }

      }
}
      
      
    
    
    
    // Print out the results
    
    printf("R = \n");
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++) {
            printf("%9.6g ", a[4*j + 4*n*i + 1]);
        }
        printf("\n");
    }
    printf("\n");
    
    printf("Q = \n");
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < m; j++) {
            printf("%9.6g ", Q[4*j + 4*n*i + 3]);
        }
        printf("\n");
    }
    printf("\n");
    
    for(int i = 0; i < n; i++){
      F[i] = Q[n*(n-1) + i];
    }
    
    
   
   






}
