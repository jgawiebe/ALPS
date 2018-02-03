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

mat bilinear_interpolation(mat img_in, uword height, uword width);
double interpolate(double s, double e, double t);
uword bilinear(uword c00, uword c10, uword c01, uword c11, uword tx, uword ty);

//partial derivatives
mat img1_dx, img1_dy;
mat img2_dx, img2_dy;
//change matrices
mat img_z, dxz, dyz;

//M: optic_flow
void compute_derivatives(mat img1, mat img2) {

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
}

mat bilinear_interpolation(mat img_in, uword height, uword width) {

	//cout << img_in;

	mat img_out(height, width);

	//https://rosettacode.org/wiki/Bilinear_interpolation
//	uword height = img_size.n_rows;
//	uword width = img_size.n_cols;

	double gx, gy;
	//uword y = 0.0;
	double dw = (double) width;
	double dh = (double) height;

//currently ignoring left and bottom edges
	for (uword x = 0, y = 0; y < (height); x++) {

//		cout << "x: " << x;

		gx = (x / dw) * (img_in.n_cols - 1);
		gy = (y / dh) * (img_in.n_rows - 1);

//		cout << "  gx: " << gx;

		double gx_int = floor(gx);
		double gy_int = floor(gy);

		double c00 = img_in(gy_int, gx_int);
		double c10 = img_in(gy_int + 1, gx_int);
		double c01 = img_in(gy_int, gx_int + 1);
		double c11 = img_in(gy_int + 1, gx_int + 1);

		//for (uword i = 0; i < 3; i++) {
		img_out.at(x, y) = interpolate(interpolate(c00, c10, gx - gx_int),
				interpolate(c01, c11, gx - gx_int), gy - gy_int);
		//}

//		cout << "c: " << c00 << " ";
//		cout << img_out.at(x, y) << endl;

		//when x gets to end, reset col and increment row
		if (x == (width - 1)) {
			x = 0;
			y++;
		}
	}
	return img_out;
}

//uword bilinear(uword c00, uword c10, uword c01, uword c11, uword tx, uword ty) {
//	return interpolate(interpolate(c00, c10, tx), interpolate(c01, c11, tx), ty);
//}

double interpolate(double s, double e, double t) {
	return s + (e - s) * t;
}

//mat bilinear_interpolation(mat img_in, mat img_size) {
//	uword srows = img_in.n_rows / img_size.n_rows;
//	uword scols = img_in.n_cols / img_size.n_cols;
//
//	mat row_mat, col_mat;
//
//	for (uword i = 0; i < img_size.n_cols; i++) {
//		col_mat.col(i) = i;
//	}
//
//	for (uword i = 0; i < img_size.n_rows; i++) {
//		row_mat.row(i) = i;
//	}
//
//	col_mat *= scols;
//	row_mat *= srows;
//
//	mat col = floor(col_mat);
//	mat row = floor(row_mat);
//
//	//normalize values in both mats
//	for (uword i = 0; i < col.n_elem; i++) {
//		if (col(i) < 1) {
//			col(i) = 1;
//		} else if (col(i) > img_size.n_cols - 1) {
//			col(i) = img_size.n_cols - 1;
//		}
//		if (row(i) < 1) {
//			row(i) = 1;
//		} else if (row(i) > img_size.n_rows - 1) {
//			row(i) = img_size.n_rows - 1;
//		}
//	}
//
//	mat delta_c = col_mat - col;
//	mat delta_r = row_mat - row;
//
//	//vec ix1 = sub2ind(size(img_in), );
//
//	mat output_img;
//	for (uword i = 0; i < col.n_elem; i++) {
//
//	}
//
//
//	//TEMPORARY
//	return row_mat;
//
//}

//M: resolutionProcess
tuple<mat, mat, mat, mat> optical_flow(double alpha, double gamma, double omega,
		mat u, mat v, int outer_iter, int inner_iter) {
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

		//can i call this from another header?
		//3 term matrix function within psi
		e_data = psi_function(
				pow(img_z + (img2_dx % du) + (img2_dy % dv), 2)
						+ gamma * pow((img_z + (dxx % du) + (dxy % dv)), 2)
						+ pow(dyz + (dxy % du) + (dyy % dv), 2));

		e_smooth = generate_esmooth(u + du, v + dv);

		//tie(A, b) = build_matrix(A, b, img2_dx, img2_dy, img_z, dxx, dxy, dyy, dxz, dyz, e_data, (alpha * e_smooth), u, v, gamma);

//		duv = successive_overrelaxation(&fail_flag, A, duv, b, omega,
//				inner_iter, tolerance);

		if (fail_flag) {
			continue; //did not reach convergence, must try again
		}

		for (uword i = 0; i < duv.n_elem; i += 2) {
			du(i) = duv(i); //column major ordering puts values in column by column
		}

		for (uword i = 1; i < duv.n_elem; i += 2) {
			dv(i) = duv(i);
		}
	}
	return make_tuple(u, v, du, dv);
}
