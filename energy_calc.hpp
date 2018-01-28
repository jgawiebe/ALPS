/*
 energy_calc.hpp
 Jacob Wiebe & James Dolman
 Reu1: Nov 2017
 Rev2: Jan 2018
 */

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

mat psi_function(mat x) {

	//M: psiDerivative
	const double epsilon = 0.001;
	x = 1 / ((2 * (sqrt(x + epsilon))));

	return x;
}

//M: computee_smooth_brox
mat generate_esmooth(mat u, mat v) {
	double height = u.n_rows;
	double width = u.n_cols;

	uword x = 0, y = 0;

	mat e_smooth(2 * height + 1, 2 * width + 1, fill::zeros);
	//mat e_temp(e_smooth);

	mat temp(1, 2), tr_temp(1, 2); //transposition of temp matrix

	temp(0) = 1;
	temp(1) = -1;
	tr_temp = temp.t();

	mat u_dx = conv2(u, temp);
	mat v_dx = conv2(v, temp);

	mat u_dy = conv2(u, tr_temp);
	mat v_dy = conv2(v, tr_temp);

	temp(1) = 1;
	tr_temp = temp.t();

	mat u_dx2 = conv2(u_dx, temp / 2, "same"); //all "valid"
	mat v_dx2 = conv2(v_dx, temp / 2, "same");

	mat u_dy2 = conv2(u_dy, tr_temp / 2, "same");
	mat v_dy2 = conv2(v_dy, tr_temp / 2, "same");

	mat delta_ux = conv2(u_dy2, (temp / 2));      //t
	delta_ux.shed_row(u_dy2.n_rows - 1);
	mat u_pdx = pow(u_dx.submat(0, 0, u_dx.n_rows - 1, delta_ux.n_cols - 1), 2)
			+ pow(delta_ux.submat(0, 0, u_dx.n_rows - 1, delta_ux.n_cols - 1),
					2);

	mat delta_uy = conv2(u_dx2, tr_temp / 2);     //t
	delta_uy.shed_col(u_dx2.n_cols - 1);
	mat u_pdy = pow(u_dy, 2) + pow(delta_uy, 2);

	mat delta_vx = conv2(v_dy2, temp / 2);       //adding row
	delta_vx.shed_row(v_dy2.n_rows - 1);
	mat v_pdx = pow(v_dx.submat(0, 0, v_dx.n_rows - 1, delta_vx.n_cols - 1), 2)
			+ pow(delta_vx.submat(0, 0, v_dx.n_rows - 1, delta_vx.n_cols - 1),
					2);

	mat delta_vy = conv2(v_dx2, tr_temp / 2);     //adding col
	delta_vy.shed_col(v_dx2.n_cols - 1);
	mat v_pdy = pow(v_dy.submat(0, 0, delta_vy.n_rows - 1, v_dy.n_cols - 1), 2)
			+ pow(delta_vy.submat(0, 0, delta_vy.n_rows - 1, v_dy.n_cols - 1),
					2);

	mat temp_x = psi_function(u_pdx + v_pdx);
	mat temp_y = psi_function(u_pdy + v_pdy);

//	e_smooth( 1:2:end, 2:2:end ) = psiDerivative( uypd + vypd ) ;
//	e_smooth( 2:2:end, 1:2:end ) = psiDerivative( uxpd + vxpd ) ;

//	uword y = 0;
//	uword x = 1;
//	for (uword row = 0; row < e_smooth.n_rows; row += 2) {
//		for (uword col = 1; col < e_smooth.n_cols; col += 2) {
//			x = 1;
//			e_smooth.at(row, col) = temp_y.at(x, y);
//			if (x == temp_x.n_cols - 1) {
//				y++;
//				//			} else if (y == temp_y.n_rows - 1) {
//				//				return e_smooth;
//			} else {
//				x++;
//			}
//		}
//	}

//	for (uword row = 0; row < e_smooth.n_rows; row++) {
//		for (uword col = 0; col < e_smooth.n_cols; col++) {
//			if ((col % 2) && !(row % 2)) { //odd col, even row
//				e_smooth.at(row, col) = temp_y.at(x, y);
//			} else if (!(col % 2) && (row % 2)) { //even col, odd row
////				e_smooth.at(row, col) = temp_x.at(x, y);
//			} else {
//				continue; //don't increment x, y if nothing is written
//			}
//			// simple step through matrix
//			if (x == temp_x.n_cols - 1) {
//				x = 0;
//				y++;
////			} else if (y == temp_y.n_rows - 1) {
////				return e_smooth;
//			} else {
//				x++;
//			}
////			if (col == 40) {
////				e_smooth.save("mats/energy/e_smooth-c.txt", raw_ascii);
////			}
//		}
//
//	}
return e_smooth;
}
