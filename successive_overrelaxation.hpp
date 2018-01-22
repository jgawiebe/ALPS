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

tuple<mat, mat, vec> split(mat M, mat N, vec b, mat A,
		double omega);

//M: sor
//the solution 'x' is the vector 'duv'
//failure flag can also be used outside this function
<<<<<<< Upstream, based on jimmy-remote/james
vec successive_overrelaxation(uword* failure, mat A, vec b, double omega,
		int inner_iter, double tolerance) {
=======
vec successive_overrelaxation(uword* failure, mat A, vec x, vec b, double omega,
		int inner_iter, vec tolerance) {
>>>>>>> bbb520c Tried to grab files from jim's branch

<<<<<<< Upstream, based on jimmy-remote/james
	mat M, N; //temp variables for matrix splittingw
	vec error;
	vec x;
=======
	mat M, N; //temp variables for matrix splitting
>>>>>>> bbb520c Tried to grab files from jim's branch

	double norml = norm(b);
	if (norml == 0) {
		norml = 1;
	}

	vec r(b);
<<<<<<< Upstream, based on jimmy-remote/james

=======
	cout << r;

	vec error(r.n_cols);
>>>>>>> bbb520c Tried to grab files from jim's branch
	error = (norm(r) / norml);

	if (all(error < tolerance)) { //matrix is already within tolerance, done
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
		if (all(error < tolerance)) {
			break; //approximation is within tolerance
		}
	}
<<<<<<< Upstream, based on jimmy-remote/james
	//b.save("mats/sor/b-c.txt", raw_ascii);
=======
	b.save("mats/sor/b-c.txt", raw_ascii);
>>>>>>> bbb520c Tried to grab files from jim's branch

	//what does this effect? b & r aren't used anywhere after this
	b /= omega;
	r = b - (A * x);

	if (any(error > tolerance)) {
		*failure = 1; //convergence not found
	}

	return x; //solution vector (duv)
}


tuple<mat, mat, vec> split(mat M, mat N, vec b, mat A,
		double omega) {
	//omega is the relaxation scalar
<<<<<<< Upstream, based on jimmy-remote/james
    //double height = A.n_rows;
    //double width = A.n_cols;
=======
//	double height = A.n_rows;
//	double width = A.n_cols;
>>>>>>> bbb520c Tried to grab files from jim's branch
	mat diagA = diagmat(diagmat(A));

	b *= omega;
	M = omega * (trimatl(A, -1) + diagA); // -1 parameter
	N = -omega * (trimatu(A, 1) + ((1 - omega) * diagA)); //1 parameter

	return make_tuple(M, N, b);
}
