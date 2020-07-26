#include <iostream>

extern"C"{
    float horner_(float coeff[], int N, float X);
}

int main(int argc, char ** args){
  int I,N;
  const int NMAX = 10;
  float COEF[NMAX];
  float X;
  
  COEF[2] = 3;
  std::cout << horner_(COEF, 4, 2.0) <<std::endl;


  return 0;

}
