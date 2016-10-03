#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <iostream>

#include <Random.h>

#include <chrono>
using namespace std::chrono;

#define PI  3.1415926535897932

/*-----------------------------------------------------------------------*/

static int int_log2(int n);

double FFT_num_flops(int N)
{

     double Nd = (double) N;
     double logN = (double) int_log2(N);

     return (5.0*Nd-2)*logN + 2*(Nd+1);
}

static int int_log2 (int n)
{
    int k = 1;
    int log = 0;
    for(/*k=1*/; k < n; k *= 2, log++);
    if (n != (1 << log))
    {
      printf("FFT: Data length is not a power of 2!: %d ",n);
      exit(1);
    }
    return log; 
}


void FFT_bitreverse(int N, float *data) {
  /* This is the Goldrader bit-reversal algorithm */
  int n=N/2;
  int nm1 = n-1;
  int i=0; 
  int j=0;
  for (; i < nm1; i++) {

    /*int ii = 2*i; */
    int ii = i << 1;

    /*int jj = 2*j; */
    int jj = j << 1;

    /* int k = n / 2 ; */
    int k = n >> 1;

    if (i < j) {
      float tmp_real, tmp_imag;
      tmp_real = data[ii];
      tmp_imag = data[ii+1];
      data[ii]   = data[jj];
      data[ii+1] = data[jj+1];
      data[jj]   = tmp_real;
      data[jj+1] = tmp_imag; }

    while (k <= j) 
      {
        /*j = j - k ; */
        j -= k;

        /*k = k / 2 ;  */
        k >>= 1 ; 
      }
    j += k ;
  }
}

static void FFT_transform_internal (int N, float *data, int direction) {
  int n = N/2;
  int bit = 0;
  int logn;
  int dual = 1;

  if (n == 1) return;         /* Identity operation! */
  logn = int_log2(n);


  if (N == 0) return;    

  /* bit reverse the input data for decimation in time algorithm */
  FFT_bitreverse(N, data) ;

  /* apply fft recursion */
  /* this loop executed int_log2(N) times */
  for (bit = 0; bit < logn; bit++, dual *= 2) {
    float w_real = 1.0;
    float w_imag = 0.0;
    int a;
    int b;

    double theta = 2.0 * direction * PI / (2.0 * (double) dual);
    double s = sin(theta);
    double t = sin(theta / 2.0);
    double s2 = 2.0 * t * t;

    for (a=0, b = 0; b < n; b += 2 * dual) {
      int i = 2*b ;
      int j = 2*(b + dual);

      float wd_real, wd_imag;
      wd_real = data[j] ;
      wd_imag = data[j+1] ;
          
      data[j]   = data[i]   - wd_real;
      data[j+1] = data[i+1] - wd_imag;
      data[i]  = data[i] + wd_real;
      data[i+1] = data[i+1] + wd_imag;
    }
      

    /* a = 1 .. (dual-1) */
    for (a = 1; a < dual; a++) {
      /* trignometric recurrence for w-> exp(i theta) w */
      {

	float tmp_real, tmp_imag;

	tmp_real = w_real - s * w_imag - s2 * w_real;

	  

	tmp_imag = w_imag + s * w_real - s2 * w_imag;


	w_real = tmp_real;
	w_imag = tmp_imag;
      }
      for (b = 0; b < n; b += 2 * dual) {
	int i = 2*(b + a);
	int j = 2*(b + a + dual);

	float z1_real, z1_imag;
	z1_real = data[j];
	z1_imag = data[j+1];
              
	float wd_real, wd_imag;

	wd_real = w_real * z1_real - w_imag * z1_imag;

	  

	wd_imag = w_real * z1_imag + w_imag * z1_real;

	  
	data[j]   = data[i]   - wd_real;
	data[j+1] = data[i+1] - wd_imag;
	data[i]  = data[i] + wd_real;
	data[i+1] = data[i+1] + wd_imag;
      }
    }
  }
}

void FFT_transform(int N, float *data)
{
    FFT_transform_internal(N, data, -1);
}


void FFT_inverse(int N, float *data)
{
    int n = N/2;
    double norm = 0.0;
    int i=0;
    FFT_transform_internal(N, data, +1);

    /* Normalize */

    norm=1/((double) n);
    
    for(i=0; i<N; i++)
      data[i] = data[i] * norm;
  
}

int main() {  
  /* initialize FFT data as complex (N real/img pairs) */
  int RANDOM_SEED[10] = {101010, 110110, 110010, 111010, 110010, 110001, 110011,
                         100001, 101001, 101101};
  int N = 16;
  for(int j = 0; j < 10; j ++){
    high_resolution_clock::time_point t_start = high_resolution_clock::now();
    Random R = new_Random_seed(RANDOM_SEED[j]);

    int twoN = 2*N;
    double *x = RandomVector(twoN, R);
    long cycles = 1;



    float InstableX[twoN];
    for(int i = 0; i < twoN; i ++){
      InstableX[i] = x[i];
    }


    int i=0;

    for (i = 0; i < cycles; i ++){
      FFT_transform(twoN, InstableX);     /* forward transform */
      FFT_inverse(twoN, InstableX);       /* backward transform */
    }
    high_resolution_clock::time_point t_end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>( t_end - t_start).count();
    std::cout << "duration " << duration <<   std::endl;   

    Random_delete(R);
    free(x);
  }
  
  return 0;
}


