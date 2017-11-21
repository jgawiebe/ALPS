/*
gradient.hpp
Jacob Wiebe & James Dolman
Rev1: Nov 2017
*/

#include <iostream>
//#include <armadillo>
//#include "conversion.hpp"
#include "gradient.hpp"

//using namespace std;
//using namespace arma;

//compile with: g++ test.cpp -o test -O2 -larmadillo; ./test
int main(){


  mat A(50, 50, fill::randu);
  mat B, C;
  B.ones(50);
  C.ones(50);

  gradient(A, B, C);


  return 0;
}
