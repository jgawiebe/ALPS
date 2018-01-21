/*
gradient.hpp
Jacob Wiebe & James Dolman
Rev1: Nov 2017
*/

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

//M: gaussDeriv
void gradient (mat *x_deriv, mat *y_deriv, mat *I){
  int sigma = 1;
  int limit = 1000;
  int variance = 1;
  int denominator = 2;

  vec numerator, derivative, temp;

//	function assumes variance is constant, otherwise use: variance = sigma * sigma;
	//	denominator = 2 * variance;
	//	derivative = derivative % (temp / variance);



  //M: gaussDeriv
  temp = linspace(-limit, limit, (2 * limit + 1));
  numerator = temp % temp; //returns element-wise product

  //returns vector of 0s because exp(negative) is very close to 0
  derivative = exp(-numerator/denominator) / pow((datum::pi * denominator), 0.5);
  derivative = derivative % temp;

//  *x_deriv = conv2(*I, derivative, "same");
//  *y_deriv = conv2(*I, derivative.t(), "same");
}
