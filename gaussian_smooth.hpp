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
	//int sigma = 1;
	int mask_size = 100;
	double thresh = 1e-3; //THIS IS CORRECT
	
	double scale_factor = (1 / scale);

	uvec limit;

	mat smooth_img;

	if (1/scale < 0.1) {
		smooth_img  = img;
	}

	mat grid = linspace(-mask_size, mask_size, (2 * mask_size) + 1);

	grid = 1
			/ ((sqrt(2 * datum::pi) * scale_factor)
					* exp(pow(-grid, 2) / (2 * pow(scale_factor, 2))));
	grid.save("mats/g_smooth/grid_init-c.txt", raw_ascii);

	limit = find(abs(grid) > thresh); //limit should only be whole numbers
	limit.save("mats/g_smooth/limit-c.txt", raw_ascii);

	grid = grid(limit); //use limit as indexes for grid elements
	grid.save("mats/g_smooth/grid-c.txt", raw_ascii);

	grid /= accu(grid); //element-wise division of sum of the column

	//rowvec grid_v = conv_to<rowvec>::from(grid);

	//M: smooth_img = conv2(grid, grid, img, "same");

	vec temp;
	rowvec temp_r; // = conv_to<rowvec>::from(temp);
	rowvec temp_row;
	rowvec grid_row(grid.size());
	mat img_build;

	//conv2(u,v,A) first convolves each column of A with the vector u, and then it convolves each row of the result with the vector v.
	//VECS ARE DIFFERENT SIZE. HOW DOES MATLAB DO THIS??
	for (uword i = 0; i < img.n_cols; i++) {
		temp = conv(grid, img.col(i));
		img_build.insert_cols(i, temp);
	}
	for (uword i = 0; i < img.n_rows; i++) {
		temp = conv(grid, img_build.row(i));
		temp_r = conv_to<rowvec>::from(temp);
		smooth_img.insert_rows(i, temp_r);
	}

	smooth_img.save("mats/g_smooth/img_smooth-c.txt", raw_ascii);
//
//	//M: gaussianRescaling
//	//assuming bilinear
//	smooth_img.resize(smooth_img.n_rows * scale_factor,
//			smooth_img.n_cols * scale_factor);



return smooth_img;
}
