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
	//conv2(u,v,A) first convolves each column of A with the vector u, and then it convolves each row of the result with the vector v.
	for(int i = 0; i < img.n_cols; i++){
		conv2(grid, img.col(i), "same");
	}
	for(int i = 0; i < img.n_rows; i++){
		conv2(grid, img.row(i), "same");
	}
	//smooth_img = conv2(grid, grid, img, "same"); //cannot perform this with armadillo
	smooth_img = imresize(g_smooth, scale_factor, 'bilinear', 0);

return smooth_img;
}
