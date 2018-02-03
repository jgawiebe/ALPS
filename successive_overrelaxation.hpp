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

tuple<sp_mat, sp_mat, vec> split(sp_mat M, sp_mat N, vec b, sp_mat A,
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

	sp_mat M, N; //temp variables for matrix splittingw
	vec error;
	vec x(b.n_rows , fill::zeros);

	failure = 0; //init fail to false

	double norml = norm(b);
	if (norml == 0) {
		norml = 1;
	}

	vec r(b);
	//note: error is always 1 as norml = norm(b) and r = b.
	error = (norm(r) / norml);
	cout<<"Test line 48"<<endl;
	if (all(error < tolerance)) { //matrix is already within tolerance, done
		mat fail_mat(size(A), fill::zeros);
		duv = x;
		return make_tuple(duv,failure);
	}
	cout<<"Test line 54"<<endl;
	//Need A to be square to be called by split, in sor
	//Had to make an sp_mat for this to work...

	//M, N are outputs; b is an inout
	tie(M, N, b) = split(M, N, b, A, omega);
	cout<<"Out split in successive_overrelaxation"<<endl;
	b.save("mats/test_SOR/Outputs/bPostSplit-c", raw_ascii);
	mat approx;
	mat tmpM;
	//continue to perform approximations until max iterations or accuracy is below the tolerance level
	for (uword i = 0; i < inner_iter; i++) {
		vec x_initial = x;
	    approx = (N * x) + b;
		cout<<"HERE 1 "<<endl;
		tmpM = (mat)M;
		cout<<"HERE 2 "<<endl;
		x = solve(tmpM, approx);
		cout<<"HERE 3 "<<endl;
		error = (norm(x - x_initial) / norm(x));
		if (all(error <= tolerance)) {
			break; //approximation is within tolerance
		}

		 cout<<"Iteration "<<i<<" Error: "<<error<<endl;
	}


	//what does this effect? b & r aren't used anywhere after this
	//just commenting them out for now
    //b /= omega;
	//r = b - (A * x);
	cout<<"Test line 78"<<endl;

	if (any(error > tolerance)) {
		failure = 1; //convergence not found set failure to true
	}
	duv = x;
	return make_tuple(duv,failure); //solution vector (duv)
}


tuple<sp_mat, sp_mat, vec> split(sp_mat M, sp_mat N, vec b, sp_mat A,
		double omega) {
			//omega is the relaxation scalar
			//double height = A.n_rows;
			//double width = A.n_cols;
	cout<<"In split"<<endl;

	sp_mat diagA = diagmat(A);
	sp_mat lwrDiagA = trimatl(A);
	sp_mat uprDiagA = trimatu(A);
	lwrDiagA = lwrDiagA - diagA; //diagonal in lwrdiagA removed
	uprDiagA = uprDiagA - diagA; //diagonal in uprDiagA removed


	b *= omega;
	M = (omega * lwrDiagA) + diagA; // -1 parameter
	N = (-omega * uprDiagA) + ((1 - omega) * diagA); //1 parameter

	return make_tuple(M, N, b);
}
