//gradient.hpp holds the function gradient(). 


#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

//This function produces the spacial gradient of an inputted matrix I. 
//This gradient is the derivative of matrix I in both the x and
//y domain, hence there are two outputs x_deriv and y_deriv. 

tuple<mat, mat> gradient(mat I) {
	
	int limit = 1000;
	int denominator = 2;
	double thresh = 1e-6;

	vec numerator, derivative, temp;

	temp = linspace(-limit, limit, (2 * limit + 1));
	numerator = temp % temp; //returns element-wise product
	derivative = exp(-numerator / denominator) / pow((datum::pi * denominator), 0.5);
	derivative = -derivative % temp;
	derivative = derivative(find(abs(derivative) > thresh));

	mat x_deriv, y_deriv;

	x_deriv = conv2(I, derivative.t(), "same");
	y_deriv = conv2(I, derivative, "same");

	return make_tuple(x_deriv, y_deriv);
}
