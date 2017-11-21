/*
gradient.hpp
Jacob Wiebe & James Dolman
Rev1: Nov 2017
*/
#include <iostream>
#include <armadillo>
#include "conversion.hpp"

//put these in generic file
//#define Q 2048 //this is 2^11 (fractional bits)
//typedef std::vector<int> ivector; //simplifies writing int vector


using namespace std;
using namespace arma;

//always 3 matrices (unsure of what type yet)
template <class T>
void gradient (T &x_deriv, T &y_deriv, T &I){
  int sigma = 1;
  int limit = 1000;
  int variance = 1;
  int denominator = 2;

  vec numerator;
  vec derivative;
  vec temp;
  ivector itemp;
  ivector ideriv;

  //function assumes variance is constant, otherwise use:
  //variance = sigma*sigma;
  //denominator = 2 * variance;
  //derivative = derivative % (temp / variance);

  //gaussDeriv
  //http://arma.sourceforge.net/docs.html#linspace
  temp = linspace(-limit, limit, (2 * limit + 1));
  numerator = temp % temp; //returns element-wise product

  //returns vector of 0s because exp(negative) is very close to 0
  derivative = exp(-numerator/denominator) / pow((datum::pi * denominator), 0.5);

  derivative = derivative % temp; //can't perform % with std vectors
  //ideriv = convertv(derivative); //convert from vec to v<int>


  //for testing:
  // for(int i = 0; i < derivative.size(); i++) {
  //   cout << derivative[i] << endl;
  // }

  x_deriv = conv2( I, derivative, "same" ) ;
  y_deriv = conv2( I, derivative.t(), "same" ) ;

  return;
}
