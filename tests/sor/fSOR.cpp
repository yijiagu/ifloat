#include <fstream>

#include <stdio.h>
#include <Random.h> 

#include <iostream>

#include <chrono>
using namespace std::chrono;

double SOR_num_flops(int M, int N, int num_iterations)
{
  double Md = (double) M;
  double Nd = (double) N;
  double num_iterD = (double) num_iterations;

  return (Md-1)*(Nd-1)*num_iterD*6.0;
}

void SOR_execute(int M, int N, float omega, float G[][100], int 
		 num_iterations)
{

  float omega_over_four = omega * 0.25;
  float one_minus_omega = 1.0 - omega;

  /* update interior points */

  int Mm1 = M-1;
  int Nm1 = N-1; 
  int p;
  int i;
  int j;
  float *Gi;
  float *Gim1;
  float *Gip1;

  for (p=0; p<num_iterations; p++){
    for (i=1; i<Mm1; i++){
      Gi = G[i];
      Gim1 = G[i-1];
      Gip1 = G[i+1];
      for (j=1; j<Nm1; j++){
	float temp;

	temp  = Gim1[j] + Gip1[j] + Gi[j-1] + Gi[j+1];


	Gi[j] = omega_over_four * temp + one_minus_omega * Gi[j];

      }
    }		
  }
}
     
int main(){
  int RANDOM_SEED[10] = {101010, 110110, 110010, 111010, 110010, 110001, 110011,
                         100001, 101001, 101101};
  int N = 100;
  
  high_resolution_clock::time_point t_start = high_resolution_clock::now();
  for(int i = 0; i < 10; i ++){
    Random R = new_Random_seed(RANDOM_SEED[i]);
    double **G = RandomMatrix(N, N, R);
    //double result = 0.0;
    float InstableG[100][100];  
    for(int i = 0; i < N; i ++){
      for(int j = 0; j < N; j ++){
	InstableG[i][j] = G[i][j];
      }
    }


    int cycles=1;


    SOR_execute(N, N, 1.25, InstableG, cycles);

    Random_delete(R);
    free(G);
  }
  
  high_resolution_clock::time_point t_end = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>( t_end - t_start).count();
  std::cout << "duration " << duration <<   std::endl; 
}       
