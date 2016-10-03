/*
*   LDL-t ALGORITHM 6.5
*
*   To factor the positive definite n by n matrix A into LDL**T,
*   where L is a lower triangular matrix with ones along the diagonal
*   and D is a diagonal matrix with positive entries on the
*   diagonal.
*
*   INPUT:   the dimension n; entries A(I,J), 1<=I, J<=n of A.
*
*   OUTPUT:  the entries L(I,J), 1<=J<I, 1<=I<=N of L and D(I),
*            1<=I<=n of D.
*/

#include<stdio.h>
#include<math.h>

#include <Random.h>
#include <iostream>

#define true 1
#define false 0

#include <chrono>
using namespace std::chrono;

void INPUT(FILE *, double [][10]);

int main()
{
  int N,I,J,K,OK;

  FILE *INP;
  char NAME[30];
  printf("Input the file name in the form - drive:name.ext\n");
  printf("for example:   A:DATA.DTA\n");
  scanf("%s", NAME);
  INP = fopen(NAME, "r");

  N = 8;
  OK = 1;
  for(int k = 0; k < 10; k ++){


    high_resolution_clock::time_point t_start = high_resolution_clock::now();
    double A_float[10][10];
    INPUT(INP, A_float);
    
    float A[N][N];
    
    for(int i = 0; i < N; i ++){
      for(int j = 0; j < N; j ++){
	A[i][j] = A_float[i][j];
      }
    }


    float D[N],V[N];
  
    if (OK) {
      /* STEP 1 */
      for (I=1; I<=N; I++) { 
	/* STEP 2 */
	for (J=1; J<=I-1; J++)
	  V[J-1] = A[I-1][J-1] * D[J-1];
	/* STEP 3 */
	D[I-1] = A[I-1][I-1];
	for (J=1; J<=I-1; J++){

	  D[I-1] = D[I-1] - A[I-1][J-1] * V[J-1];


	} 
	/* STEP 4 */
	for (J=I+1; J<=N; J++) {

	  for (K=1; K<=I-1; K++){	    

	    A[J-1][I-1] = A[J-1][I-1] - A[J-1][K-1] * V[K-1];

	  }
	  

	  A[J-1][I-1] = A[J-1][I-1] / D[I-1];
	}  
      }

      high_resolution_clock::time_point t_end = high_resolution_clock::now();
      auto duration = duration_cast<microseconds>( t_end - t_start).count();
      std::cout << "duration " << duration << std::endl; 



    }
  }
  
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
