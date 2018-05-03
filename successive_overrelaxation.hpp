//successive_overrelaxation.hpp holds two functions, split() and 
//successive_overrelation().

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;
//This function is a service function for successive_overrelaxation(). 
//It splits up the cooeffient matrix A so that it can be used by the 
//iterative solver (success_overrelaxation()).
tuple<sp_mat, sp_mat, vec> split(sp_mat M, sp_mat N, vec b, sp_mat A,
	double omega);

//This function solves a system of linear equations A(x) = b where 
//the solution is x. The function is an iterative solver where it makes 
//estimations of x until it reaches its max number of tries or the 
//error between estimations of x become so small that is is known to have 
//found the solution. 
tuple<vec, uword> successive_overrelaxation(sp_mat A,
	vec b, double omega, uword inner_iter, double tolerance) {

	sp_mat M, N; //temp variables for matrix splittingw
	double error;
	vec x(b.n_rows, fill::zeros);

	int failure = 0; //init fail to false
	int* fail = &failure;

	double norml = norm(b);
	if (norml == 0) {
		norml = 1;
	}

	vec r(b);
	error = (norm(r) / norml);
	if (error < tolerance) { //matrix is already within tolerance, done
		mat fail_mat(A.n_rows, A.n_cols, fill::zeros);
		return make_tuple(x, *fail);
	}
	
	tie(M, N, b) = split(M, N, b, A, omega);
	
	mat approx;
	mat tmpM = (mat)M;
	vec x_initial = x;
	
	//continue to perform approximations until max iterations or accuracy is below the tolerance level
	for (uword i = 0; i < inner_iter; i++) {
		printf("SOR iteration %d - error is: %1.11f\r", i, error); 
	    x_initial = x;
		approx = (N * x) + b;
		x = solve(tmpM, approx);
		error = (norm(x - x_initial) / norm(x));
		if (error <= tolerance) {
			break; //approximation is within tolerance
		}
	}

	if (error > tolerance) {
		failure = 1; //convergence not found set failure to true
	}
	else {
		cout << "                                               " << endl;
		cout << "Convergence reached" << endl;
		failure = 0;
	}

	return make_tuple(x, *fail); //solution vector (duv)
}


tuple<sp_mat, sp_mat, vec> split(sp_mat M, sp_mat N, vec b, sp_mat A,
	double omega) {

	sp_mat diagA = diagmat(A);
	sp_mat lwrDiagA = trimatl(A);
	sp_mat uprDiagA = trimatu(A);
	lwrDiagA = lwrDiagA - diagA; //diagonal in lwrdiagA removed
	uprDiagA = uprDiagA - diagA; //diagonal in uprDiagA removed

	b *= omega;
	M = (omega * lwrDiagA) + diagA; 
	N = (-omega * uprDiagA) + ((1 - omega) * diagA); 

	return make_tuple(M, N, b);
}
