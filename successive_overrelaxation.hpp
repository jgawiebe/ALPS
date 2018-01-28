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

tuple<mat, mat, vec> split(mat M, mat N, vec b, sp_mat A,
		double omega);

//M: sor
//the solution 'x' is the vector 'duv'
//failure flag can also be used outside this function


tuple<vec, uword> successive_overrelaxation( vec duv, uword failure, sp_mat A,
		vec b, double omega,uword inner_iter, double tolerance) {


/*	  A            REAL matrix
	  duv          REAL initial guess vector
      b            REAL right hand side vector
	  omega        REAL relaxation scalar
      inner_iter   INTEGER maximum number of iterations
	  tolerance    REAL error tolerance */

	cout<<"In successive_overrelaxation"<<endl;

	mat M, N; //temp variables for matrix splittingw
	vec error;
	vec x;
	failure = 0; //init fail to false

	double norml = norm(b);
	if (norml == 0) {
		norml = 1;
	}

	vec r(b);

	error = (norm(r) / norml);

	if (all(error < tolerance)) { //matrix is already within tolerance, done
		mat fail_mat(size(A), fill::zeros);
		duv = x;
		return make_tuple(duv,failure);
	}

	//Need A to be square to be called by split, in sor
	//Had to make an sp_mat for this to work...

	//M, N are outputs; b is an inout
	tie(M, N, b) = split(M, N, b, A, omega);
	cout<<"Out split in successive_overrelaxation"<<endl;

	//continue to perform approximations until max iterations or accuracy is below the tolerance level
	for (uword i = 0; i < inner_iter; i++) {
		vec x_initial = x;
		mat approx = (N * x) + b;

		x = solve(M, approx);
		error = (norm(x - x_initial) / norm(x));
		if (all(error < tolerance)) {
			break; //approximation is within tolerance
		}
	}
	//b.save("mats/sor/b-c.txt", raw_ascii);

	//what does this effect? b & r aren't used anywhere after this
	b /= omega;
	r = b - (A * x);


	if (any(error > tolerance)) {
		failure = 1; //convergence not found set failure to true
	}
	duv = x;
	return make_tuple(duv,failure); //solution vector (duv)
}


tuple<mat, mat, vec> split(mat M, mat N, vec b, sp_mat A,
		double omega) {
			//omega is the relaxation scalar
			//double height = A.n_rows;
			//double width = A.n_cols;
	cout<<"In split"<<endl;
	mat diagA = diagmat(diagmat(A));

	b *= omega;
	M = omega * (trimatl(A, -1) + diagA); // -1 parameter
	N = -omega * (trimatu(A, 1) + ((1 - omega) * diagA)); //1 parameter

	return make_tuple(M, N, b);
}
