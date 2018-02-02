/*
 energy_calc.hpp
 Jacob Wiebe & James Dolman
 Rev1: Nov 2017
 */

#include <iostream>
#include <armadillo>
#include "optical_flow.hpp"
#include "gaussian_smooth.hpp"

using namespace std;
using namespace arma;

//compile with: g++ <funtion>.cpp -o test -O2 -larmadillo; ./test

//M: optic_flow_brox
int main() {

	//non-matrices: alpha, dt, gamma, ht, i, num_levels, wt
	double alpha = 30.0, gamma = 80.0, omega = 1.8; //check omega value
	int num_levels = 40, outer_iter = 3, inner_iter = 500;
	
	double scale_factor = pow(0.95, num_levels);
	
	mat du, dv;

	//**get images using opencv**
	mat image1, image2;

	//get size of image
	int height = image1.n_rows;
	int width = image1.n_cols;
	
	//perform guassian scaling on images
	mat img1 = g_smooth(image1, scale_factor);
	mat img2 = g_smooth(image2, scale_factor);
	
	//define u and v matrices
	mat u(img1.n_rows, img1.n_cols, fill::zeros);
	mat v(img1.n_rows, img1.n_cols, fill::zeros);

	//check this loop
	for (int i = 0; i < num_levels; i++) {
		
		//resolution increases with each pyramid level
		scale_factor = pow(0.95, i);
		
		//derivatives are stored as optical_flow globals
		compute_derivatives(img1, img2);
		
		//perform optical_flow to get du and dv
		tie(u, v, du, dv) = optical_flow(alpha, gamma, omega, u, v, outer_iter,
				inner_iter);
		
		//add incremental change in x and y domain
		u = u + du;
		v = v + dv;
		
		//scale images to current level of pyramid
		img1 = g_smooth(image1, scale_factor);
		img2 = g_smooth(image2, scale_factor);
		
		//M: u = imresize( u, [size(im1_hr, 1), size(im1_hr, 2)], 'bilinear' );
		
		//resize flow to the current resolution (assuming bilinear)
		u.load("mats/energy/u-m.txt");
		v.load("mats/energy/v-m.txt");

		u.resize(img1.n_rows, img1.n_cols);
		v.resize(img1.n_rows, img1.n_cols);

		u.save("mats/energy/u-c.txt");
		v.save("mats/energy/v-c.txt");
	}

	return 1;
}

