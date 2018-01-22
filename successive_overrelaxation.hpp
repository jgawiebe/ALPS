/*
 energy_calc.hpp
 Jacob Wiebe & James Dolman
 Rev1: Nov 2017
 Rev2: Jan 2018
 */

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

tuple<mat, mat, vec> split(mat M, mat N, vec b, mat A, double omega);

//M: sor
//the solution 'x' is the vector 'duv'
//failure flag can also be used outside this function
vec successive_overrelaxation(uword* failure, mat A, vec x, vec b, double omega,
		int inner_iter, vec tolerance) {

	mat M, N; //temp variables for matrix splitting

	double norml = norm(b);
	if (norml == 0) {
		norml = 1;
	}

	vec r(b);
	cout << r;

	vec error(r.n_cols);
	error = (norm(r) / norml);

	if ( all(error < tolerance) ) { //matrix is already within tolerance, done
		mat fail_mat(size(A), fill::zeros);
		return x;
	}

	//M, N are outputs; b is an inout
	tie(M, N, b) = split(M, N, b, A, omega);

	//continue to perform approximations until max iterations or accuracy is below the tolerance level
	for (int i = 0; i < inner_iter; i++) {
		vec x_initial = x;
		mat approx = (N * x) + b;

		x = solve(M, approx);
		error = (norm(x - x_initial) / norm(x));
		if ( all(error < tolerance) ) {
			break; //approximation is within tolerance
		}
	}
	b.save("mats/sor/b-c.txt", raw_ascii);

	//what does this effect? b & r aren't used anywhere after this
	b /= omega;
	r = b - (A * x);

	if ( any(error > tolerance) ) {
		*failure = 1; //convergence not found
	}

	return x; //solution vector (duv)
}

//complete
tuple<mat, mat, vec> split(mat M, mat N, vec b, mat A, double omega) {
	//omega is the relaxation scalar
//	double height = A.n_rows;
//	double width = A.n_cols;
	mat diagA = diagmat(diagmat(A));

	b *= omega;
	M = omega * (trimatl(A, -1) + diagA);
	N = -omega * (trimatu(A, 1) + ((1 - omega) * diagA));

	return make_tuple(M, N, b);
}
