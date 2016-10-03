#include <iostream>

#include <InstableFloat.h>
#include <InstableEnv.h>
#include <InstableMath.h>

//Dot Product 3-Vectors
InstableFloat dot3(InstableFloat *a, InstableFloat *b){     
  InstableFloat r;
  r = a[0] * b[0] + a[1] * b[1] 
  + a[2] * b[2];
  return r;
}

int raySphere(InstableFloat *r, InstableFloat *s, InstableFloat radiusSq) {
  InstableFloat A, B , C, D;
  InstableFloat temp0, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8;
  radiusSq = radiusSq;
  InstableEnv::open_label(1);
  A = dot3(r,r);
  InstableEnv::close_label();
  std::cout << A << std::endl;

  InstableEnv::open_label(2);
  B = -2.0 * dot3(s,r);
  InstableEnv::close_label();
  std::cout << B << std::endl;

  InstableEnv::open_label(3);
  C = dot3(s,s) - radiusSq;
  InstableEnv::close_label();
  std::cout << C << std::endl;
 
  InstableEnv::open_label(4);
  D = B*B - 4 * A * C;
  InstableEnv::close_label();
  
  D = D + 2;
  temp0 = sqrt(D);

  std::cout << temp0 << std::endl;
  return 0;
} 

int main() {

  InstableFloat s[3] = {-10.4194345474,	-15,	-14};
  InstableFloat r[3] = { -10.998046875,	-16,	-15};
  InstableFloat radiusSq (0.015625);
  return raySphere(r, s, radiusSq);
} 
