/*
energy_calc.hpp
Jacob Wiebe & James Dolman
Rev1: Nov 2017
*/

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

//M: sor
vec successive_overrelaxation (int& failure, mat A, vec x, vec duv, vec b, double omega, int inner_iter, double tolerance){

	double norml = norm(b);
	if (norml == 0) {
		norml = 1;
	}

	vec r = b-A*x;
	error = norm(r)/norml;
	if (error<tolerance) {
		return; //report error
	}



	//non-matrices: gamma, ht, wt
	//vectors: vals, b, cols, ind, pdfaltsumu/v
	//sparse: A

  return duv;
}


void split (mat& M, mat& N, vec&b, mat A, double b, double w) {
	double height = A.n_rows;
	double width = A.n_cols;


}
