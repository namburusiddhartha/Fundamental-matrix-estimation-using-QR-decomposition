#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <immintrin.h>
#include <math.h>

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
  int              m,
  int              n,
  int              k,
  double* restrict a,
  double* restrict b,
  double* restrict c
){
 
 //TODO: Implement a small mxk * kxn matrix multiplication here
 
   //int rn = n % 4;
   //n = (n/4)*4;
   int num_points = k * 4 * 4;
   
   //printf("%d\n", num_points);
   int size = 0;
   
   double X1[16];
   double X2[16];
   for(int i=0;i<16;i++){
      X1[i]=rand()%1000;
      X2[i]=rand()%1000;
   }
   
   for(int i = 0; i < num_points; i ++){
       
       double A[8][9];
       double det = 0;
       
       for (int j = 0; j < 8; j++) {
       A[j][0] = X1[j]*X2[j];
       A[j][1] = X1[j]*X2[8 + j];
       A[j][2] = X1[8 + j]*X2[j];
       A[j][3] = X1[8 + j]*X2[8 + j];
       A[j][4] = X1[j]*X1[j];
       A[j][5] = X2[j]*X2[j];
       A[j][6] = X1[8 + j]*X1[8 + j];
       A[j][7] = X2[8 + j]*X2[8 + j];
       A[j][8] = 1;   
       
       }
       
       double B[9][8];
       
       for (int k = 0; k < 8; k ++){
         for (int l = 0; l < 9; l++){
          B[l][k] = A[k][l];         
         }
       }
       
       double C[9][9];
       
       for (int k = 0; k < 9; k ++){
         for (int l = 0; l < 9; l++){
           for (int m = 0; m < 8; m++){
             C[l][k] += A[m][k] * B[l][m];
           }         
         }
       }
       
       int n = 9;
       

      size=n;
      double matrix[25][25],size,d;
      int i,j;

    for (i=0;i<size;i++)
    {
        for (j=0;j<size;j++)
        {
             
             matrix[i][j] = C[i][j];
        }

    }
    d=determinant(matrix,size);
    if (d==0)
    {
        d = 3490;
    }
    else
        cofactor(matrix,size);
        
        
        
        
    double v[9];
    for(int j=0;j<9;j++){
      v[j]=rand()%1000;
   }
   double V[9];
   
   for (int j = 0; j < 9; j ++){
     for (int k = 0; k < 9; k ++) {
       V[j] += (matrix[j][k] * v[k]) / v[k];
     }
   }
   
   double eigenvalue = V[0];
   
   for (int j = 0; j < 1000; j ++){
     double V[9];
   
     for (int j = 0; j < 9; j ++){
       for (int k = 0; k < 9; k ++) {
         V[j] += (matrix[j][k] * v[k]) / v[k];
       }
     }
     eigenvalue = V[0];
     
     double norm;
     
     for (int k = 0; k < 9; k ++){
         norm += v[k]*v[k];
     }
     
     for (int j = 0; j < 9; j ++){
       for (int k = 0; k < 9; k ++) {
         v[j] += (matrix[j][k] * v[k]) / norm;
       }
     }
     
     
   
   }
   
   double finaleigvalue = V[0];
   
   double error;
   double brr1, brr2, brr3;
   
   for(int k = 0; k < 8; k++){
     brr1 =  v[0] * X1[k] + v[1] * X1[j + 8] + v[2];
     brr2 =  v[3] * X1[k] + v[4] * X1[j + 8] + v[5];
     brr3 =  v[6] * X1[k] + v[7] * X1[j + 8] + v[8];
     
     error += (brr1 - X2[j])*(brr1 - X2[j]) + (brr2 - X2[j + 8])*(brr2 - X2[j + 8]) + (brr2 - 1)*(brr2 - 1);
   
   }
   
   
    


   }
   
   
   
     
   
   // Scalar multiply when n is not a multiple of 4
//   for (int i = 0; i != m * (n + rn); ++i) {
//      printf("%lf \n", c[i]);
//    }
//   
//   for (int p = 0; p != k; ++p)
//    {
//      for (int i = 0; i != m; ++i)
//	{
//	  for (int j = 0; j < n + rn; ++j)
//	    {
//	      c[m*(n) + j] += a[p*m + i] * b[p*n + j];
//
//	    }
//	}

//    }
   

   
 
}
