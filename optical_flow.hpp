//optical_flow.hpp holds two functions, compute_derivatives() and 
//optical_flow(). 

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

extern double* jacobi_serial(double *val, long long *lrow, long long *lcol, double *b, const unsigned int n, const unsigned int x_size, const unsigned int max_iter);
extern double* parallel_solve_v1(double *val, long long *lrow, long long *lcol, double *b, int *fail, const unsigned int n, const unsigned int x_size, const unsigned int max_iter, double tol);
extern double* parallel_solve_v2(double *val, long long *lrow, long long *lcol, double *b, const unsigned int n, const unsigned int x_size, const unsigned int max_iter, double tol);

//partial derivatives
mat img1_dx, img1_dy;
mat img2_dx, img2_dy;
//change matrices
mat img_z, dxz, dyz;

//This function produces derivative variables for the two 
//inputted matrixes. 
void compute_derivatives(mat img1, mat img2) {
	//produce primary derivatives of img1
	tie(img1_dx, img1_dy) = gradient(img2);
	//produce primary derivatives of img2
	tie(img2_dx, img2_dy) = gradient(img1);
	//difference between image values
	img_z = img2 - img1;
	//difference between image x partial derivatives
	dxz =  img1_dx - img2_dx;
	//difference between image y partial derivatives
	dyz = img1_dy - img2_dy;
}

//This function produces an estimatation of a pixel shift map, expressed 
//as matrixes u and v for a given level in the gaussian pyramid. It has two modes,
//parrallel and serial. 
tuple<mat, mat, mat, mat> optical_flow(double alpha, double gamma, double omega, double tolerance,
	mat u, mat v, int outer_iter, int inner_iter, bool parallel_mode) {
	
	int fail_val = 1;
	int* fail_flag = &fail_val;

	mat dxx, dxy;
	mat dyx, dyy;
	mat du(arma::size(u), fill::zeros);
	mat dv(arma::size(v), fill::zeros);
	mat e_smooth, e_data(du.n_rows, du.n_cols, fill::zeros);

	//get ht and wt from size of img_z
	uword ht = u.n_rows;
	uword wt = u.n_cols;

	vec duv(ht * wt * 2, fill::zeros);
	vec b(ht * wt * 2, fill::zeros);

	// if parallel
	urowvec lrow, lcol;
	vec val;
	int nnz = 0, m = 0;

	//if sequential
	sp_mat A;

	//get second derivatives of img2_dx
	tie(dxx, dxy) = gradient(img1_dx);

	//get second derivatives of img2_dy
	tie(dyx, dyy) = gradient(img1_dy);
	
	for (int i = 0; i < outer_iter + 1; i++) {

		//3 term matrix function within psi
		e_data = psi_function(
			pow(img_z + (img1_dx % du) + (img1_dy % dv), 2)
			+ gamma * (pow((dxz + (dxx % du) + (dxy % dv)), 2)
			+ pow(dyz + (dxy % du) + (dyy % dv), 2)));
		
		e_smooth = generate_esmooth(u + du, v + dv);

		if (parallel_mode) {
			tie(lrow, lcol, val, b) = build_matrix(img1_dx, img1_dy, img_z, dxx, dxy, dyy, dxz, dyz, e_data, alpha*e_smooth, u, v, gamma);
			nnz = (int)val.n_elem;
			m = (int)b.n_elem;
			double* x = jacobi_serial(val.memptr(), (long long*)lrow.memptr(), (long long*)lcol.memptr(), b.memptr(), nnz, m, inner_iter);
			vec x_in(x, m);
			duv = x_in;
		}
		else {
			//serial execution
			tie(A, b) = build_sparse(img2_dx, img2_dy, img_z, dxx, dxy, dyy, dxz, dyz, e_data, alpha*e_smooth, u, v, gamma);
			tie(duv, *fail_flag) = successive_overrelaxation(A, b, omega, inner_iter, tolerance);
		}
		
		uword ik = 0;
		for (uword i = 0; i < du.n_elem; i ++) {	
			du(i) = duv(ik); //column major ordering puts values in column by column
			ik += 2;
		}

		ik = 1;
		for (uword i = 1; i < dv.n_elem; i ++) {
			dv(i) = duv(ik);
			ik += 2;
		}

		return make_tuple(u, v, du, dv);
	}
	
	//only reaches this point if failure occured
	cout << "ALGORTIHM COMPLETE: FAILURE\nPress any key to exit" << endl;
	u.save("u-hori.txt", raw_ascii);
	v.save("v-hori.txt", raw_ascii);
	cout << "Press Enter" << endl;
	cin.get();
	exit(EXIT_FAILURE);
}
