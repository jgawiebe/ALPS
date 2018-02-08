#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

tuple<sp_mat, sp_mat, vec> split(sp_mat M, sp_mat N, vec b, sp_mat A,
	double omega);

//M: sor
//the solution 'x' is the vector 'duv'
//failure flag can also be used outside this function
tuple<vec, uword> successive_overrelaxation(sp_mat A,
	vec b, double omega, uword inner_iter, double tolerance) {


	/*	  A            REAL matrix
		  duv          REAL initial guess vector
		  b            REAL right hand side vector
		  omega        REAL relaxation scalar
		  inner_iter   INTEGER maximum number of iterations
		  tolerance    REAL error tolerance */

	cout << "Performing successive-overrelaxation >" << endl;

	sp_mat M, N; //temp variables for matrix splittingw
	vec error;
	vec x(b.n_rows, fill::zeros);

	bool failure = 0; //init fail to false

	double norml = norm(b);
	if (norml == 0) {
		norml = 1;
	}

	vec r(b);
	//note: error is always 1 as norml = norm(b) and r = b.
	error = (norm(r) / norml);
	if (all(error < tolerance)) { //matrix is already within tolerance, done
		mat fail_mat(A.n_rows, A.n_cols, fill::zeros);
		return make_tuple(x, failure);
	}
	//Need A to be square to be called by split, in sor
	//Had to make an sp_mat for this to work...

	//M, N are outputs; b is an inout
	tie(M, N, b) = split(M, N, b, A, omega);
	cout << " done" << endl;
	//b.save("mats/test_SOR/Outputs/bPostSplit-c", raw_ascii);
	mat approx;
	//mat tmpM;


	//CAN'T MAKE MAT FROM SPARSE MAT (TOO MUCH MEMORY USAGE)
	
	mat tmpM = (mat)M;
	cout << "FIX" << endl;

	//continue to perform approximations until max iterations or accuracy is below the tolerance level
	for (uword i = 0; i < inner_iter; i++) {
		vec x_initial = x;
		cout << "Making approximation... ";
		approx = (N * x) + b;
		cout << " done" << endl;
		cout << "Solving for x...";
		x = solve(tmpM, approx);
		cout << " done" << endl;
		error = (norm(x - x_initial) / norm(x));
		if (all(error <= tolerance)) {
			break; //approximation is within tolerance
		}

		cout << "Iteration " << i << " of " << inner_iter << " complete. Error is: " << error << endl;
	}


	//what does this effect? b & r aren't used anywhere after this
	//just commenting them out for now
	//b /= omega;
	//r = b - (A * x);

	if (any(error > tolerance)) {
		failure = true; //convergence not found set failure to true
	}
	else {
		cout << "> Solution vector determined" << endl;
	}
	return make_tuple(x, failure); //solution vector (duv)
}


tuple<sp_mat, sp_mat, vec> split(sp_mat M, sp_mat N, vec b, sp_mat A,
	double omega) {
	//omega is the relaxation scalar
	//double height = A.n_rows;
	//double width = A.n_cols;
	cout << "Splitting matrix...";

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