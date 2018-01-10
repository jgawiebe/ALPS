/*
guassian_smooth.hpp
Jacob Wiebe & James Dolman
Rev1: Nov 2017
*/

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

//M: guassianSmooth
mat g_smooth (mat img, double scale){
	int sigma = 1;
	int mask_size = 100;
	double thresh = 0.001;
	double scale_factor = 1;

	vec limit(13);

	if (1/scale < 0.1) {
		g_smooth = img;
	}

	mat grid = linspace(-mask_size, mask_size, (2 * mask_size) + 1);

	grid = 1 / (sqrt(2*datum::pi)) * exp(pow(-grid, 2)/2));
	limit = find(abs(grid) > abs(thresh));
	grid = grid(limit); //what does this mean?

	grid /= sum(grid);

	//ISSUE HERE
	smooth_img = conv2(grid,grid, img, "same"); //cannot perform this with armadillo

	//non-matrices: endthresh, maxMaskSize, sigma, scale factor
	//vectors: lim, grid
	//img is in 3d in matlab rgb


	//how do I make this bilinear??
	smooth_img.resize()

return smooth_img;
}
