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
	double scale_factor = (1 / scale);

	uvec limit;

	mat smooth_img(img);

	if (1/scale < 0.1) {
		smooth_img  = img;
	}

	mat grid = linspace(-mask_size, mask_size, (2 * mask_size) + 1);

	grid = 1 / (sqrt(2*datum::pi) * exp(pow(-grid, 2)/2));

//	mat thresh(grid);
//	thresh.fill(0.001);
	mat ugrid = abs(grid);
	limit = find(ugrid > thresh); //limit should only be whole numbers
	cout << "limit:\n" << limit << endl;
	grid = grid(limit); //use limit as indexes for grid elements


	//double grid_sum = sum(grid);
	grid = (grid / cumsum(grid, 1)); //element-wise division of sum of each row
	cout << "grid:\n" << grid << endl;


	//M: smooth_img = conv2(grid, grid, img, "same");
	//conv2(u,v,A) first convolves each column of A with the vector u, and then it convolves each row of the result with the vector v.
	//TEST THIS
//	for(uword i = 0; i < img.n_cols; i++){
//		conv2(grid, img.col(i), "same");
//	}
//	for(uword i = 0; i < img.n_rows; i++){
//		conv2(grid, img.row(i), "same");
//	}
	
	//M: gaussianRescaling
	//assuming bilinear
	//smooth_img.resize(smooth_img.n_rows * scale_factor, smooth_img.n_cols * scale_factor);

return smooth_img;
}
