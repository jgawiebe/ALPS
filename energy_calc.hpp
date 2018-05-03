//energy_calc.hpp holds two functions; psi_function() and
//generate_esmooth(). 

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

//This function performs a mathematical computation on an inputted 
//matrix x. 
mat psi_function(mat x) {
	const double epsilon = 0.001;
	x = 1 / ((2 * (sqrt(x + epsilon))));
	return x;
}

//This function generates the variable e_smooth. e_smooth is a matrix
//of doubles who's dimension and content is dependant on the inputted 
//u and v matrixes. 
mat generate_esmooth(mat u, mat v) {
	double height = u.n_rows;
	double width = u.n_cols;

	mat e_smooth(2 * height + 1, 2 * width + 1, fill::zeros);
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

	mat u_dx2 = conv2(u_dx, temp / 2, "same"); 
	mat v_dx2 = conv2(v_dx, temp / 2, "same");

	mat u_dy2 = conv2(u_dy, tr_temp / 2, "same");
	mat v_dy2 = conv2(v_dy, tr_temp / 2, "same");

	mat delta_ux = conv2(u_dy2, (temp / 2));      
	delta_ux.shed_row(u_dy2.n_rows - 1);
	mat u_pdx = pow(u_dx.submat(0, 0, u_dx.n_rows - 1, delta_ux.n_cols - 1), 2)
		+ pow(delta_ux.submat(0, 0, u_dx.n_rows - 1, delta_ux.n_cols - 1),
			2);

	mat delta_uy = conv2(u_dx2, tr_temp / 2);     
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

	uword y = 0;
	uword x = 0;

	for (uword row = 1; row < e_smooth.n_rows; row += 2) {
		for (uword col = 1; col < e_smooth.n_cols; col += 2) {

			e_smooth(row - 1, col) = temp_y(x, y); //vert edge case
			e_smooth(row, col - 1) = temp_x(x, y); //right edge case

			if (row == ((e_smooth.n_rows)) - 2) {
				e_smooth(row + 1, col) = temp_y(x, y); //vert edge case
			}

			if (col == ((e_smooth.n_cols)) - 2) {
				e_smooth(row, col + 1) = temp_x(x, y); //right edge case
			}

			y++;
		}
		y = 0;
		x++;
	}
	return e_smooth;
}
