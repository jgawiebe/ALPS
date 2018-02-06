/*
 energy_calc.hpp
 Jacob Wiebe & James Dolman
 Rev1: Nov 2017
 */

#include <iostream>
#include <armadillo>
#include "gradient.hpp"
#include "energy_calc.hpp"
#include "matrix_builder.hpp"
#include "successive_overrelaxation.hpp"

using namespace std;
using namespace arma;

//partial derivatives
mat img1_dx, img1_dy;
mat img2_dx, img2_dy;
//change matrices
mat img_z, dxz, dyz;

//M: optic_flow
void compute_derivatives(mat img1, mat img2) {
	cout << "Computing derivatives >" << endl;

	//produce primary derivatives of img1
	gradient(img1_dx, img1_dy, img1);

	//produce primary derivatives of img2
	gradient(img2_dx, img2_dy, img2);

	//difference between image values
	img_z = img2 - img1;
	//difference between image x partial derivatives
	dxz = img2_dx - img1_dx;
	//difference between image y partial derivatives
	dyz = img2_dy - img1_dy;

	cout << "> Derivatives saved" << endl;
}

//M: resolutionProcess
tuple<mat, mat, mat, mat> optical_flow(double alpha, double gamma, double omega,
		mat u, mat v, int outer_iter, int inner_iter) {
	cout << "Initializing optical flow procedure >"<< endl;
	bool fail_flag = 0;
	double tolerance = 1e-8;

	mat dxx, dxy;
	mat dyx, dyy;
	mat du(img_z), dv(img_z);
	mat e_smooth, e_data(du.n_rows, du.n_cols, fill::zeros);

	//get ht and wt from size of img_z
	int ht = img_z.n_rows;
	int wt = img_z.n_cols;
	uword side_length = (ht * wt * 2);

	vec duv(side_length, fill::zeros);

	mat A(side_length, side_length);

	//check the size of b
	vec b(side_length, fill::zeros);

	//get second derivatives of img2_dx
	gradient(dxx, dxy, img2_dx);

	//get second derivatives of img2_dy
	gradient(dyx, dyy, img2_dy);

	for (int i = 0; i < outer_iter; i++) {

		//3 term matrix function within psi
		cout << "Generating E_data" << endl;
		e_data = psi_function(
				pow(img_z + (img2_dx % du) + (img2_dy % dv), 2)
						+ gamma * pow((img_z + (dxx % du) + (dxy % dv)), 2)
						+ pow(dyz + (dxy % du) + (dyy % dv), 2));

		e_smooth = generate_esmooth(u + du, v + dv);

		cout << "Building matrix..." << endl;
		//tie(A, b) = build_matrix(A, b, img2_dx, img2_dy, img_z, dxx, dxy, dyy, dxz, dyz, e_data, (alpha * e_smooth), u, v, gamma);

		cout << "Performing successive-overrelaxation..." << endl;
//		duv = successive_overrelaxation(&fail_flag, A, duv, b, omega,
//				inner_iter, tolerance);

		if (fail_flag) {
			cout << "Successive-overrelaxation failed to reach convergence. Retrying." << endl;
			continue; //did not reach convergence, must try again
		}

		for (uword i = 0; i < duv.n_elem; i += 2) {
			du(i) = duv(i); //column major ordering puts values in column by column
		}

		for (uword i = 1; i < duv.n_elem; i += 2) {
			dv(i) = duv(i);
		}
	}
	cout << "> Optical flow complete" << endl;
	return make_tuple(u, v, du, dv);
}
