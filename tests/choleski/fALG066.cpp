/*
*   CHOLESKI'S ALGORITHM 6.6
*
*   To factor the positive definite n by n matrix A into LL**T,
*   where L is lower triangular.
*
*   INPUT:   the dimension n; entries A(I,J), 1<=I, J<=n of A.
*
*   OUTPUT:  the entries L(I,J), 1<=J<=I, 1<=I<=n of L.
*
*   the entries of U = L**T are U(I,J) = L(J,I), I<=J<=n, 1<=I<=n
*/

#include<stdio.h>
#include<math.h>

#include <Random.h> 
#include <iostream>

#include <chrono>
using namespace std::chrono;

#define true 1
#define false 0

void INPUT(FILE *, double [][10]);

int main()
{
  int N,I,J,K,NN,JJ,KK;
  FILE *INP;
  char NAME[30];
  printf("Input the file name in the form - drive:name.ext\n");
  printf("for example:   A:DATA.DTA\n");
  scanf("%s", NAME);
  INP = fopen(NAME, "r");

  N = 8;
 
  for(int i = 0; i < 10; i ++){
     high_resolution_clock::time_point t_start = high_resolution_clock::now();
     
     double A_float[10][10];
     
     INPUT(INP, A_float);
     
     float A[N][N];
     float S;
     
     for(int i = 0; i < N; i ++){
       for(int j = 0; j < N; j ++){
	 A[i][j] = A_float[i][j];
       }
     }
     
     A[0][0] = sqrt(A[0][0]);
      /* STEP 2 */
      for (J=2; J<=N; J++)
	A[J-1][0] = A[J-1][0] / A[0][0];
      /* STEP 3 */
      NN = N - 1;
      for (I=2; I<=NN; I++) {
         /* STEP 4 */
         KK = I - 1;
         S = 0.0;
         for (K=1; K<=KK; K++){
	   S = S - A[I-1][K-1] * A[I-1][K-1];
	 }
	 float temp;
	 temp = A[I-1][I-1] + S;
         A[I-1][I-1] = sqrt(temp);
         /* STEP 5 */
         JJ = I + 1;
         for (J=JJ; J<=N; J++) {
            S = 0.0;
            KK = I - 1;
            for (K=1; K<=KK; K++){
	      S = S - A[J-1][K-1] * A[I-1][K-1];
	    }
	    float temp2;
	    temp2 = A[J-1][I-1] + S;
            A[J-1][I-1] = temp2 / A[I-1][I-1]; 
         }  
      }
      /* STEP 6 */
      S = 0.0;
      for (K=1; K<=NN; K++){
	S = S - A[N-1][K-1] * A[N-1][K-1];
      }
      float temp3;
      temp3 = A[N-1][N-1] + S;
      A[N-1][N-1] = sqrt(temp3);

      high_resolution_clock::time_point t_end = high_resolution_clock::now();
      auto duration = duration_cast<microseconds>( t_end - t_start).count();
      std::cout << "duration " << duration << std::endl; 
   }
   fclose(INP);
   return 0;
}

void INPUT(FILE *INP, double A[][10])
{
   int I, J;
   for (I=1; I<=8; I++) {
     for (J=1; J<=8; J++) fscanf(INP, "%lf", &A[I-1][J-1]);
     fscanf(INP, "\n");
   }
}
