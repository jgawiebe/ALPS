#include <iostream>
#include <armadillo>
#include "Fixed.h"

typedef std::vector<int> ivector;
typedef numeric::Fixed<5, 11> int16;

using namespace std;
using namespace arma;

#define Q 2048 //this is 2^11 (fractional bits)

ivector convertv(vec original){
  ivector output;

  original *= Q;
  output = conv_to<ivector>::from(original);
  return output;
}
