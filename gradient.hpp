/*
gradient.hpp
Jacob Wiebe & James Dolman
Rev1: Nov 2017
 Complete: Jan 21 2018
*/

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

//M: gaussDeriv
tuple<mat, mat> gradient(mat I) {

	printf("Generating gradient of a %u by %u matrix...", I.n_rows, I.n_cols);

	//int sigma = 1;
  int limit = 1000;
	//int variance = 1;
  int denominator = 2;
	double thresh = 1e-6;

  vec numerator, derivative, temp;

//	function assumes variance is constant, otherwise use: variance = sigma * sigma;
	//	denominator = 2 * variance;
	//	derivative = derivative % (temp / variance);


  //M: gaussDeriv
  temp = linspace(-limit, limit, (2 * limit + 1));
  numerator = temp % temp; //returns element-wise product

  //returns vector of 0s because exp(negative) is very close to 0
	derivative = exp(-numerator/denominator) / pow((datum::pi * denominator), 0.5);
	derivative = -derivative % temp;
	derivative = derivative(find(abs(derivative) > thresh));

	//derivative.save("mats/gradient/deriv-c.txt", raw_ascii);

	mat x_deriv, y_deriv;

	x_deriv = conv2(I, derivative.t(), "same");
	y_deriv = conv2(I, derivative, "same");

	cout << " done" << endl;

	return make_tuple(x_deriv, y_deriv);
}
