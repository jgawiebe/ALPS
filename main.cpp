/*
 energy_calc.hpp
 Jacob Wiebe & James Dolman
 Rev1: Nov 2017
 */

#include <iostream>
#include <armadillo>
#include "optical_flow.hpp"

using namespace std;
using namespace arma;

//compile with: g++ <funtion>.cpp -o test -O2 -larmadillo; ./test

//M: optic_flow_brox
int main() {

	//non-matrices: alpha, dt, gamma, ht, i, num_levels, wt
	double alpha = 30.0, gamma = 80.0;
	int num_levels = 40, outer_iter = 3, inner_iter = 500;
	
	double scale_factor = pow(0.95, num_levels);
	
	//**get images using opencv**

	//get size of image
	int height = image1.n_rows;
	int width = image1.n_cols;
	
	//perform guassian scaling on images
	mat img1 = gaussian_smooth(image1, scale_factor);
	mat img2 = gaussian_smooth(image2, scale_factor);
	
	//define u and v matrices
	mat u(size(img1), fill::zeros);
	mat v(size(img1), fill::zeros);

	//check this loop
	for (int i = 0; i < num_levels; i++) {
		
		//resolution increases with each pyramid level
		scale_factor = pow(0.95, i);
		
		//derivatives are stored as optical_flow globals
		compute_derivatives(img1, img2);
		
		//perform optical_flow to get du and dv
		optical_flow(alpha, gamma, omega, u, v, outer_iter, inner_iter);
		
		//add incremental change in x and y domain
		u = u + du;
		v = v + dv;
		
		//scale images to current level of pyramid
		img1 = guassian_smooth(image1, scale_factor);
		img2 = guassian_smooth(image2, scale_factor);
		
		//resize flow to the current resolution
		//M: u = imresize( u, [size(im1_hr, 1), size(im1_hr, 2)], 'bilinear' );
		//u.resize()
	}

	return 1;
}

