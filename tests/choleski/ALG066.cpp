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

#include <Random.h> 
#include <InstableFloat.h>
#include <InstableEnv.h>
#include <InstableMath.h>

#include <chrono>
using namespace std::chrono;

#define true 1
#define false 0
#define dimention 15

void INPUT(FILE *, double [][dimention]);

int main()
{
  int N,I,J,K,NN,JJ,KK;
  FILE *INP;
  INP = fopen("choleski.txt", "r");

  N = dimention;
  
   for(int i = 0; i < 10; i ++){
     high_resolution_clock::time_point t_start = high_resolution_clock::now();
     
     double A_float[dimention][dimention];
     
     INPUT(INP, A_float);
     
     InstableFloat A[dimention][dimention];
     InstableFloat S;
     
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
	     
	   InstableEnv::open_label(1);
	   S = S - A[I-1][K-1] * A[I-1][K-1];
	   InstableEnv::close_label();
	 }
	 
	 InstableFloat temp;
	 temp = A[I-1][I-1] + S;
         A[I-1][I-1] = sqrt(temp);
         /* STEP 5 */
         JJ = I + 1;
         for (J=JJ; J<=N; J++) {
            S = 0.0;
            KK = I - 1;
            for (K=1; K<=KK; K++){
	      InstableEnv::open_label(2);
	      S = S - A[J-1][K-1] * A[I-1][K-1];
	      InstableEnv::close_label();

	    }

	    InstableFloat temp2;
	    temp2 = A[J-1][I-1] + S;
            A[J-1][I-1] = temp2 / A[I-1][I-1];
         }  
      }

      /* STEP 6 */
      S = 0.0;
      for (K=1; K<=NN; K++){
	InstableEnv::open_label(3);
	S = S - A[N-1][K-1] * A[N-1][K-1];
	InstableEnv::close_label();
      }
      InstableFloat temp3;
      temp3 = A[N-1][N-1] + S;

      A[N-1][N-1] = sqrt(temp3);

      high_resolution_clock::time_point t_end = high_resolution_clock::now();
      auto duration = duration_cast<microseconds>( t_end - t_start).count();
      std::cout << "duration " << duration << std::endl; 

      float bound_width = 0;
      InstableFloat r = A[0][0];
      int row, column;
      for (I=1; I<=N; I++) {
	for (J=1; J<=I; J++){
	  FP_Interval v_bound = A[I-1][J-1].getVolatileBound();
	  if(v_bound.upper() - v_bound.lower() > bound_width){
	    row = I-1;
	    column = J-1;
	    r = A[I-1][J-1];
	    bound_width = v_bound.upper() - v_bound.lower();
	  }
	}
      }
      std::cout << row << column << std::endl;
      std::cout << r << std::endl;
   }
   fclose(INP);
   return 0;
}

void INPUT(FILE *INP, double A[][dimention])
{
   int I, J;
   for (I=1; I<=dimention; I++) {
     for (J=1; J<=dimention; J++) fscanf(INP, "%lf", &A[I-1][J-1]);
     fscanf(INP, "\n");
   }
}
