
/*
*   CROUT REDUCTION FOR TRIDIAGONAL LINEAR SYSTEMS ALGORITHM 6.7
*
*   To solve the n x n linear system
*
*   E1:  A[1,1] X[1] + A[1,2] X[2]                  = A[1,n+1]
*   E2:  A[2,1] X[1] + A[2,2] X[2] + A[2,3] X[3]    = A[2,n+1]
*   :
*   .
*   E(n):          A[n,n-1] X[n-1] + A[n,n] X[n]    = A[n,n+1]
*
*   INPUT:   the dimension n; the entries of A.
*
*   OUTPUT:  the solution X(1), ..., X(N).
*/
#include <InstableFloat.h>
#include <InstableEnv.h>
#include <Random.h>

#include<stdio.h>

#include <chrono>
using namespace std::chrono;

#define true 1
#define false 0

int main()
{
   InstableFloat A[10],B[10],C[10],BB[10],Z[10],X[10];
   int N,I,NN,II;

   N = 10;
   int RANDOM_SEED[10] = { 100001, 101010,  
			  111010,  110001, 110011,
			   101001, 101101,110010, 110110,110010};
   
   for(int i = 0; i < 1; i ++){
     Random R = new_Random_seed(RANDOM_SEED[i]);
     double *A_float = RandomVector(N, R);
     double *B_float = RandomVector(N, R);
     double *C_float = RandomVector(N, R);
     double *BB_float = RandomVector(N, R);

     high_resolution_clock::time_point t_start = high_resolution_clock::now();

     for (I=1; I<=N; I++)
       A[I-1] = A_float[I-1];

     /* the lower sub-diagonal A(I,I-1) is stored
	in B(I), 2 <= I <= n */ 
     for (I=2; I<=N; I++)
       B[I-1] = B_float[I-1];
     
     /* the upper sub-diagonal A(I,I+1) is stored
	in C(I), 1 <= I <= n-1 */ 
     NN = N - 1;
     for (I=1; I<=NN; I++)
       C[I-1] = C_float[I-1];
     
     for (I=1; I<=N; I++)
       BB[I-1] = BB_float[I - 1];
     


      /* Steps 1-3 set up and solve LZ = B   */
      /* STEP 1 */
      /* the entries of U overwrite C and
            the entries of L overwrite A */
      C[0] = C[0] / A[0]; 
      Z[0] = BB[0] / A[0];
      /* STEP 2 */
      for (I=2; I<=NN; I++) {
	InstableEnv::open_label(1);
	A[I-1] = A[I-1] - B[I-1] * C[I-2];
	InstableEnv::close_label();
	
	C[I-1] = C[I-1] / A[I-1];
	InstableFloat t;
	InstableEnv::open_label(2);
	t = BB[I-1]-B[I-1]*Z[I-2];
	InstableEnv::close_label();
	Z[I-1] = t/A[I-1];
      }
      /* STEP 3 */
      InstableEnv::open_label(3);
      A[N-1] = A[N-1] - B[N-1] * C[N-2];
      InstableEnv::close_label();
      InstableFloat t;
      InstableEnv::open_label(4);
      t = BB[N-1]-B[N-1]*Z[N-2];
      InstableEnv::close_label();
      Z[N-1] = t/A[N-1];
      /* STEP 4 */
      /* STEPS 4, 5 solve UX = Z */
      X[N-1] = Z[N-1];
      /* STEP 5 */
      for (II=1; II<=NN; II++) {
	I = NN - II + 1;
	InstableEnv::open_label(5);
	X[I-1] = Z[I-1] - C[I-1] * X[I];
	InstableEnv::close_label();
      }

      high_resolution_clock::time_point t_end = high_resolution_clock::now();
      auto duration = duration_cast<microseconds>( t_end - t_start).count();
      std::cout << "duration " << duration << std::endl;
      std::cout << X[0] << std::endl;
      /* STEP 6 */
      float bound_width = 0;
      InstableFloat r=X[0];
      int index;
      for(int i = 0; i < N; i ++){
	  FP_Interval v_bound = X[i].getVolatileBound();
	  if(v_bound.upper() - v_bound.lower() > bound_width){
	    index = i;
	    r = X[i];
	    bound_width = v_bound.upper() - v_bound.lower();
	  }
      }

      std::cout << index << std::endl;
      std::cout << r << std::endl;

      Random_delete(R);
      free(A_float);
      free(B_float);
      free(C_float);
      free(BB_float);
   }
   return 0;
}



