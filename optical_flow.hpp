/*
 energy_calc.hpp
 Jacob Wiebe & James Dolman
 Rev1: Nov 2017
 */

#include <iostream>
#include <armadillo>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "gradient.hpp"
#include "energy_calc.hpp"
#include "matrix_builder.hpp"
#include "successive_overrelaxation.hpp"
#include "red_black_wrapper.hpp"

using namespace std;
using namespace arma;

//partial derivatives
mat img1_dx, img1_dy;
mat img2_dx, img2_dy;
//change matrices
mat img_z, dxz, dyz;


//extern void sparse_solver(double* A, double* duv);

//M: optic_flow
void compute_derivatives(mat img1, mat img2) {
	cout << "Computing derivatives >" << endl;

	//produce primary derivatives of img1
	tie(img1_dx, img1_dy) = gradient(img1);

	//produce primary derivatives of img2
	tie(img2_dx, img2_dy) = gradient(img2);

	//difference between image values
	img_z = img2 - img1;
	//difference between image x partial derivatives
	dxz = img2_dx - img1_dx;
	//difference between image y partial derivatives
	dyz = img2_dy - img1_dy;
}

//M: resolutionProcess
tuple<mat, mat, mat, mat> optical_flow(double alpha, double gamma, double omega, double tolerance,
	mat u, mat v, int outer_iter, int inner_iter) {
	cout << "Running optical flow procedure >" << endl;
	bool fail_flag = false;
	//double tolerance = 1e-8;
	//double tolerance = 0.5;

	mat dxx, dxy;
	mat dyx, dyy;
	mat du(img_z), dv(img_z);
	mat e_smooth, e_data(du.n_rows, du.n_cols, fill::zeros);

	//get ht and wt from size of img_z
	uword ht = img_z.n_rows;
	uword wt = img_z.n_cols;

	vec duv(ht*wt * 2, fill::zeros);

	sp_mat A;

	cout << "Calculating second derivatives >" << endl;
	//get second derivatives of img2_dx
	tie(dxx, dxy) = gradient(img2_dx);

	//get second derivatives of img2_dy
	tie(dyx, dyy) = gradient(img2_dy);

	for (int i = 0; i < outer_iter + 1; i++) {

		//3 term matrix function within psi
		cout << "Generating E_data >" << endl;
		e_data = psi_function(
			pow(img_z + (img2_dx % du) + (img2_dy % dv), 2)
			+ gamma * pow((img_z + (dxx % du) + (dxy % dv)), 2)
			+ pow(dyz + (dxy % du) + (dyy % dv), 2));
		cout << "e-data done" << endl;

		e_smooth = generate_esmooth(u + du, v + dv);

		tie(A, duv) = build_matrix(img2_dx, img2_dy, img_z, dxx, dxy, dyy, dxz, dyz, e_data, (alpha * e_smooth), u, v, gamma);

		double* duv_raw = duv.memptr();
		
		mat A_dense(A);
		double* A_raw = A_dense.memptr();
		int b_size = duv.n_elem;
		int A_height = A_dense.n_rows;
		int A_width = A_dense.n_cols;

		cv::Mat Acv = to_cvmat(A_dense);

		//cout << "kernel here>>>>>>>>>>>>>" << endl;
		//red_black_sor();
		sparse_solver(Acv, duv_raw, b_size, A_height, A_width);
		cin.get();
		//sparse_solver();
		tie(duv, fail_flag) = successive_overrelaxation(A, duv, omega, inner_iter, tolerance);

		if (fail_flag && i < outer_iter) {
			cout << "> Successive-overrelaxation failed to reach convergence -> Retrying (attempt " << i << " of " << outer_iter - 1 << ")\n\n" << endl;
			continue; //did not reach convergence, must try again
		}

		for (uword i = 0; i < du.n_elem; i += 2) {
			du(i) = duv(i); //column major ordering puts values in column by column
		}

		for (uword i = 1; i < dv.n_elem; i += 2) {
			dv(i) = duv(i);
		}

		if(!fail_flag) {
			duv.save("b-built.txt", raw_ascii);
			return make_tuple(u, v, du, dv);
		}
	}
	cout << "ALGORTIHM COMPLETE: FAILURE\nPress any key to exit" << endl;
	u.save("mats/fin/u-c.txt", raw_ascii);
	v.save("mats/fin/v-c.txt", raw_ascii);
	cin.get();
	exit(EXIT_FAILURE);
}
