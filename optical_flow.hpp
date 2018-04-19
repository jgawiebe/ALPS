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


extern double* jacobi_serial(double *val, long long *lrow, long long *lcol, double *b, const unsigned int n, const unsigned int x_size, const unsigned int max_iter);
extern double* parallel_solve_v1(double *val, long long *lrow, long long *lcol, double *b, int *fail, const unsigned int n, const unsigned int x_size, const unsigned int max_iter, double tol);
extern double* parallel_solve_v2(double *val, long long *lrow, long long *lcol, double *b, const unsigned int n, const unsigned int x_size, const unsigned int max_iter, double tol);


//partial derivatives
mat img1_dx, img1_dy;
mat img2_dx, img2_dy;
//change matrices
mat img_z, dxz, dyz;


//extern void sparse_solver(double* A, double* duv);

//M: optic_flow
void compute_derivatives(mat img1, mat img2) {
	//cout << "Computing derivatives >" << endl;

	//MISORDER IS DUE TO ERROR IN MATLAB BEING FOLLOWED HERE

	//produce primary derivatives of img1
	tie(img1_dx, img1_dy) = gradient(img2);

	//I KNOW THESE IMAGE ASSIGNMENTS SHOULB BE IM1_DX WITH GRADIENT IMG1
	//BUT THATS NOT HOW THE MATLAB DID IT

	//produce primary derivatives of img2
	tie(img2_dx, img2_dy) = gradient(img1);

	//difference between image values
	img_z = img2 - img1;
	//difference between image x partial derivatives
	dxz =  img1_dx - img2_dx;
	//difference between image y partial derivatives
	dyz = img1_dy - img2_dy;
	//WE ARE NOT DONE HERERE
	//YES WE ARE, PASS IN THIS FORMAT
	
}

//M: resolutionProcess




tuple<mat, mat, mat, mat> optical_flow(double alpha, double gamma, double omega, double tolerance,
	mat u, mat v, int outer_iter, int inner_iter, bool parallel_mode) {
	//cout << "Running optical flow procedure >" << endl;
	
	int fail_val = 1;
	int* fail_flag = &fail_val; //should be init to true
	//bool parallel_mode = true; //temporary!
	//double tolerance = 1e-8;
	//double tolerance = 0.5;

	mat dxx, dxy;
	mat dyx, dyy;
	//mat du, dv;
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

	//cout << "Calculating second derivatives >" << endl;
	//get second derivatives of img2_dx
	tie(dxx, dxy) = gradient(img1_dx);

	//get second derivatives of img2_dy
	tie(dyx, dyy) = gradient(img1_dy);
	
	//SECOND DERIVATIVES PASS THE EXECUTION CHAIN


	for (int i = 0; i < outer_iter + 1; i++) {

		//3 term matrix function within psi
		//cout << "Generating E_data >" << endl;
		e_data = psi_function(
			pow(img_z + (img1_dx % du) + (img1_dy % dv), 2)
			+ gamma * (pow((dxz + (dxx % du) + (dxy % dv)), 2)
			+ pow(dyz + (dxy % du) + (dyy % dv), 2)));
		//cout << "e-data done" << endl;
		//INPUTS INJECTED AFTER RESIZE VIA LOAD
		//PASSED PSI GENERATES E DATA WELL, HIGH ACCURACY
	
		e_smooth = generate_esmooth(u + du, v + dv);

		///////////////////////////////////////////////Testing Matrix Inputs

		//cout << " img1_dx " << endl << img1_dx << endl;
		//cout << " img1_dy " << endl << img1_dy << endl;
		//cout << " img_z " << endl << img_z << endl;
		//cout << " dxx " << endl << dxx << endl;
		//cout << " dxy " << endl << dxy << endl;
		//cout << " dyy " << endl << dyy << endl;
		//cout << " dxz " << endl << dxz << endl;
		//cout << " dyz " << endl << dyz << endl;
		//cout << " e_data " << endl << e_data << endl;
		//cout << " alpha*e_smooth " << endl << alpha*e_smooth << endl;
		//cout << " u " << endl << u << endl;
		//cout << " v " << endl << v << endl;


		//cin.get();


		/////////////////////////////////////////////

		if (parallel_mode) {
			//PARALLEL 
			tie(lrow, lcol, val, b) = build_matrix(img1_dx, img1_dy, img_z, dxx, dxy, dyy, dxz, dyz, e_data, alpha*e_smooth, u, v, gamma);
			//cout << "DONE MATRIX BUILDER" << endl;
			//cin.get();
			
		/*	b.save("b.txt", raw_ascii);
			lrow.save("row.txt", raw_ascii);
			lcol.save("col.txt", raw_ascii);
			val.save("val.txt", raw_ascii);*/

			//for (int ib = 0; ib < 20; ib++) {
			//	printf("(%d, ", lrow[ib]);
			//	printf("%d)\n", lcol[ib]);
			//}
			//cout << lrow << " row " << lcol << " col " endl; // << val << " val" << endl;
			//cout << lrow << "row" << endl;

			//cin.get();

			nnz = (int)val.n_elem;
			m = (int)b.n_elem;

			

			//double* x = parallel_solve_v1(val.memptr(), (long long*)lrow.memptr(), (long long*)lcol.memptr(), b.memptr(), fail_flag, nnz, m, inner_iter, tolerance);
			//cout << "DONE SOLVE" << endl;
			//cin.get();
			
			double* x = jacobi_serial(val.memptr(), (long long*)lrow.memptr(), (long long*)lcol.memptr(), b.memptr(), nnz, m, inner_iter);
			vec x_in(x, m);
			//x_in.save("xin_cuda-c.txt", raw_ascii);
			duv = x_in;
			

			//GETTING GOOD DATA HERE USING PARALLEL SOR IMPLEMENTATION

		}
		else {
			//SERIAL
			tie(A, b) = build_sparse(img2_dx, img2_dy, img_z, dxx, dxy, dyy, dxz, dyz, e_data, alpha*e_smooth, u, v, gamma);
			//cout << "DONE MATRIX BUILDER" << endl;
			//cin.get();
			tie(duv, *fail_flag) = successive_overrelaxation(A, b, omega, inner_iter, tolerance);
			//cout << "DONE SOR" << endl;
			//cin.get();
			//cout << duv << endl;

		}

		//THE PARALLEL IMPLMENTATION (JACOBI) DOES NOT USE THE ERROR FLAGS AT THIS TIME
		//if ( *fail_flag == 1 && i < outer_iter) {
		//	cout << "> Successive-overrelaxation failed to reach convergence -> Retrying (attempt " << i << " of " << outer_iter - 1 << ")\n\n" << endl;
		//	continue; //did not reach convergence, must try again
		//}

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

		//THE PARALLEL IMPLMENTATION (JACOBI) DOES NOT USE THE ERROR FLAGS AT THIS TIME
		//if(*fail_flag == 0) {
		//	//duv.save("FULLTEST-duv.txt", raw_ascii);
			return make_tuple(u, v, du, dv);
		//}
	}
	cout << "ALGORTIHM COMPLETE: FAILURE\nPress any key to exit" << endl;
	u.save("u-hori.txt", raw_ascii);
	v.save("v-hori.txt", raw_ascii);
	cin.get();
	
	
	exit(EXIT_FAILURE);
}
