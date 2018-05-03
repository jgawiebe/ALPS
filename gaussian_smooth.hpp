//guassian_smooth.hpp holds the g_smooth() function.

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

//g_smooth() outputs a matrix that has been reduced in resolution
//or "smoothed". The inputs are the matrix to be smoothed and the 
//scale which the matrix is to be smoothed down to. The dimensions 
//of the matrix do not change after smoothing - simply the resolution. 

mat g_smooth (mat img, double scale){
	
	int mask_size = 100;
	float thresh = 1e-3; //THIS IS CORRECT
	
	double scale_factor = (1 / scale);

	uvec limit;

	mat smooth_img;

	if (scale_factor < 0.1) {
		smooth_img  = img;
	}

	mat grid = linspace(-mask_size, mask_size, (2 * mask_size) + 1);
	grid = 1 / ((sqrt(2 * datum::pi) * scale_factor)
		* exp(pow(-grid, 2) / (2 * pow(scale_factor, 2))));

	limit = find(abs(grid) > thresh); //use limit as indexes for grid elements
	grid = grid(limit); 
	grid /= accu(grid); //element-wise division of sum of the column

	vec temp;
	rowvec temp_r; // = conv_to<rowvec>::from(temp);
	mat img_build;

	// conv2(u,v,A) first convolves each column of A with the vector u,
	// and then it convolves each row of the result with the vector v.
        for (uword i = 0; i < img.n_cols; i++) {
		temp = conv(grid, img.col(i));
		img_build.insert_cols(i, temp);
	}

	for (uword i = 0; i < img_build.n_rows; i++) {
		temp = conv(grid, img_build.row(i));
		temp_r = conv_to<rowvec>::from(temp);
		smooth_img.insert_rows(i, temp_r);
	}
	
	//This done to trim the matrix down to old size of img.
	//Just the math of the conv in armadillo adds some fat that
	//needs to be trimmed.
	uword grid_length = grid.n_rows;
	smooth_img = smooth_img.submat((grid_length-1)/2,(grid_length-1)/2,(smooth_img.n_rows-1)
		         - ((grid_length)/2),(smooth_img.n_cols-1) - ((grid_length)/2));

        return smooth_img;
}
