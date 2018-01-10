/*
energy_calc.hpp
Jacob Wiebe & James Dolman
Rev1: Nov 2017
*/

#include <iostream>
#include <armadillo>
#include "gradient.hpp"
#include "energy_calc.hpp"

using namespace std;
using namespace arma;

//partial derivatives
mat img1_dx, img1_dy;
mat img2_dx, img2_dy;
//change matrices
mat img_z, img_dxz, img_dyz;

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
void optical_flow (double alpha, double gamma, double omega, mat u, mat v, int outer_iter, int inner_iter){
	int fail_flag = 0;

	mat dxx, dxy;
	mat dyx, dyy;

	//**get ht and wt from size of img_z

	mat du(ht, wt, fill::zeros);
	mat dv(ht, wt, fill::zeros);
	vec duv(ht * wt * 2, fill::zeros);
	vec tolerance(ht * wt * 2);

	//DO MATRIX SIZES NEED TO BE INITIALIZED??
	mat e_data(size(du), fill::zeros);
	mat e_smooth;
	mat e_init;


	//get second derivatives of img2_dx
	gradient(dxx, dxy, img2_dx);

	//get second derivatives of img2_dy
	gradient(dyx, dyy, img2_dy);

	tolerance.fill(1e-8); //fill all values with 1e-8

	for(int i = 0; i < outer_iter; i++){
		mat temp

		//can i call this from another header?
		//3 term matrix function within psi
		e_data = psi_function(
				pow(img_z + (img2_dx % du) + (img2_dy % dv), 2) +
				gamma * pow((img_z + (dxx % du) + (dxy % dv) ), 2) +
				pow(img_dyz + (dxy % du) + (dyy % dv), 2)
		);
		e_smooth = generate_esmooth(u + du, v + dv);

		matrix_builder(mat& A, vec& b, img2_dx, img2_dy, img_z, dxx, dxy, dyy, img_dxz, e_data, (alpha * e_smooth), u, v, gamma);

		successive_overrelaxation(fail_flag, A, duv, b, omega, inner_iter, tolerance);

		//seperate duv into du and dv

		//deal with fail_flag
	}
}
