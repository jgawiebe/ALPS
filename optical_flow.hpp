/*
energy_calc.hpp
Jacob Wiebe & James Dolman
Rev1: Nov 2017
*/

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

//partial derivatives
mat img1_dx, img1_dy;
mat img2_dx, img2_dy;
//change matrices
mat img_z, img_dxz, img dyz;

//M: optic_flow
void compute_derivatives(mat img1, mat img2){

	//produce primary derivatives of img1
	gradient(img1_dx, img1_dy, img1);

	//produce primary derivatives of img2
	gradient(img2_dx, img2_dy, img2);

	//difference between image values
	img_z = img2 - img1;
	//difference between image x partial derivatives
	img_dxz = img2_dx - img1_dx;
	//difference between image y partial derivatives
	img_dyz = img2_dy - img1_dy;
}

//M: resolutionProcess
alpha, gamma, omega, uinit, vinit, outer_iter, inner_iter)
void optical_flow (double alpha, double gamma, double omega, mat u, mat v, int outer_iter, int inner_iter){
	mat img2_ddx


	//non-matrices: alpha, dt, err, flag, gamma, ht, i, inner_iter, it, omega, outer_iter, wt
	//vectors: b, du, dv, duv, tol
	//sparse: A

  return;
}
