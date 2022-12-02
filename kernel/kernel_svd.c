#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <immintrin.h>
#include <math.h>
//#include <algorithm>
//using namespace std;
#define FMADD(dest, src1, src2) \
  __asm__ __volatile__ (      \
  "vfmadd231pd %[rsrc1], %[rsrc2], %[rdest]\n"  \
    : [rdest] "+x"(dest)     \
    : [rsrc1] "x" (src1) , [rsrc2] "x"(src2));
    
#define MATMUL(m, n, k, a, b, c) \

//sad

double determinant(double[25][25],double);
void cofactor(double[25][25],double);
void transpose(double[25][25],double[25][25],double);



double determinant(double matrix[25][25],double size)
{
    double s=1,det=0,m_minor[25][25];
    int i,j,m,n,c;
    if (size==1)
    {
        return (matrix[0][0]);
    }
    else
    {
        det=0;
        for (c=0;c<size;c++)
        {
            m=0;
            n=0;
            for (i=0;i<size;i++)
            {
                for (j=0;j<size;j++)
                {
                    m_minor[i][j]=0;
                    if (i != 0 && j != c)
                    {
                       m_minor[m][n]=matrix[i][j];
                       if (n<(size-2))
                          n++;
                       else
                       {
                           n=0;
                           m++;
                       }
                    }
                }
            }
            det=det + s * (matrix[0][c] * determinant(m_minor,size-1));
            s=-1 * s;
        }
    }
 
    return (det);
}
 
 /*calculate cofactor of matrix*/
void cofactor(double matrix[25][25],double size)
{
     double m_cofactor[25][25],matrix_cofactor[25][25];
     int p,q,m,n,i,j;
     for (q=0;q<size;q++)
     {
         for (p=0;p<size;p++)
         {
             m=0;
             n=0;
             for (i=0;i<size;i++)
             {
                 for (j=0;j<size;j++)
                 {
                     if (i != q && j != p)
                     {
                        m_cofactor[m][n]=matrix[i][j];
                        if (n<(size-2))
                           n++;
                        else
                        {
                            n=0;
                            m++;
                        }
                     }
                 }
             }
             matrix_cofactor[q][p]=pow(-1,q + p) * determinant(m_cofactor,size-1);
         }
     }
     transpose(matrix,matrix_cofactor,size);
}

/*Finding transpose of cofactor of matrix*/ 
void transpose(double matrix[25][25],double matrix_cofactor[25][25],double size)
{
     int i,j;
     double m_transpose[25][25],m_inverse[25][25],d;
 
     for (i=0;i<size;i++)
     {
         for (j=0;j<size;j++)
         {
             m_transpose[i][j]=matrix_cofactor[j][i];
         }
     }
     d=determinant(matrix,size);
     for (i=0;i<size;i++)
     {
         for (j=0;j<size;j++)
         {
             m_inverse[i][j]=m_transpose[i][j] / d;
         }
     }
}
  

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
   const int m = 8, n = 9;
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

   //double A[8][3] = 
   // {
   //     {2., 5., 3.},
   //     {1., 2., 1.},
   //     {4., 1., 1.},
   //     {3., 5., 2.},
   //     {5., 3., 1.},
   //     {4., 5., 5.},
   //     {2., 4., 2.},
   //     {2., 2., 5.},
   // };

   double matrixFor1D[m][n];
   memcpy ( &matrixFor1D, &A, sizeof(A));
   double past_U[q];  
   double past_singularval = 0.0;
   double pastV[p];
   double C[n][n];
   double AT[n][m];
   double currentV[p];
   double norm = 0.0;
   int iterations;
   double lastv[p];
   double current_U[m];
   for(int iter = 0; iter < n; ++iter){
       printf("\n\n %d \n\n", iter);
       if(iter != 0){
       for(int r = 0; r < m; r++){
        for(int c = 0; c < n; c++){
          A[r][c] -= past_singularval *past_U[r]*pastV[c]  ;
         }}}
         printf("\n\n A = \n");
        // for(int r = 0; r < m; r++){
        //for(int c = 0; c < n; c++){
        //   printf("%f ", A[r][c]);
        //   }
        //   printf("\n");
        // }
        // printf("\n");

         // svd_1d(double A[m][n])
         
             norm = 0.0;
             for(int i = 0; i < n; i++){
              
              // Change this to uniform random distribution
              currentV[i] = 0.0;
              currentV[i] = ((double) rand() / (RAND_MAX)) * 2 - 1 ; //(double)0.5; //(rand()%1000000)/1000000 ;
              norm += currentV[i]*currentV[i];
               }

	      printf("Random values \n");
	      for(int i = 0; i < n; i++){
		printf("%f \t", currentV[i]);
	      }
	      printf("\n");
              norm = sqrt(norm);
              for(int i = 0; i < n; i++){
              currentV[i] = currentV[i]/norm;
               }    
             for (int r = 0; r < n; r++){
               for (int c = 0; c < m; c++){
                AT[r][c] = A[c][r];
                 }
               }

	     for(int r = 0; r < n; r++){
                for(int c = 0; c < n; c++){
                   C[r][c] = 0;
                }
                }

             for (int r = 0; r < n; r++){
               for (int c = 0; c < n; c++){
                 for (int k = 0; k < m; k++){
                   C[r][c] += AT[r][k] * A[k][c];
                   }
                 }
               }


      
         iterations = 0;
         while (1)
         {
         
             iterations ++;
             memcpy ( &lastv, &currentV, sizeof(currentV));
             //auto lastv = currentV;
             // 3x3 @ 3x1
             norm = 0.0;
             double sum = 0.0;
             for (int r = 0; r < n; r++){
                 for (int k = 0; k < n; k++){
                       currentV[r] += C[r][k] * lastv[k];
                       }
                     norm += currentV[r]*currentV[r]; 
                      }
     
             norm = sqrt(norm);
             for (int r = 0; r < n; r++){
                 currentV[r] = currentV[r]/norm;
                 // Calculating sum
                 sum += currentV[r]*lastv[r];
                 }
                 //break;   
             printf("iter = %d  sum = %f \n ", iterations, sum);
             if (iterations > 20) break;
          }
          //printf(" currentV[0] = %f \t ", currentV[0]);
          //printf(" currentV[1] %f \t ", currentV[1]);
          //printf(" currentV[2] %f \t \n", currentV[2]);

	  for(int i = 0; i < n; i++){
		printf("V = %f \t", currentV[i]);
		
	  }
          printf("End V \n");

          // We have Updated currentV
          double singularval = 0.0;
          double unorm = 0.0;
          for (int r = 0; r < m; r++){
             for (int c = 0; c < n; c++){
             
             current_U[r] += matrixFor1D[r][c]*currentV[c];
             }
             unorm += current_U[r]*current_U[r];
          }
          unorm = sqrt(unorm);
          singularval = unorm; 
          for (int r = 0; r < m; r++){
          current_U[r] = current_U[r]/unorm;
          }
   printf("eigen val = %f", singularval);
   memcpy ( &past_U, &current_U, sizeof(current_U));
   memcpy ( &pastV, &currentV, sizeof(currentV));
   memcpy ( &past_singularval, &singularval, sizeof(singularval));
   
   }
 
}
