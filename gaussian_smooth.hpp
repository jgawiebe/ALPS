/*
guassian_smooth.hpp
Jacob Wiebe & James Dolman
Rev1: Nov 2017
Rev2: Feb 2018 TESTING PASSED 
*/

#include <iostream>
#include <armadillo>
//#include "nvblas.h"

using namespace std;
using namespace arma;

//M: guassianSmooth
mat g_smooth (mat img, double scale){
	//int sigma = 1;
	int mask_size = 100;
	float thresh = 1e-3; //THIS IS CORRECT
	
	double scale_factor = (1 / scale);

	uvec limit;

	mat smooth_img;

	if (scale_factor < 0.1) {
		smooth_img  = img;
	}

	mat grid = linspace(-mask_size, mask_size, (2 * mask_size) + 1);

	grid = 1
			/ ((sqrt(2 * datum::pi) * scale_factor)
					* exp(pow(-grid, 2) / (2 * pow(scale_factor, 2))));
	//assignment is correct dimension and really close to what the matlab has
	//grid.save("mats/test_gaussian_smooth/Outputs/grid-c", raw_ascii);

	limit = find(abs(grid) > thresh); //limit should only be whole numbers
	//consistently 1 under which is incorrect. Should be 1 higher Ah, but
	//it is correct because it is 0 indexed.
	//limit.save("mats/test_gaussian_smooth/Outputs/limit-c", raw_ascii);

	grid = grid(limit); //use limit as indexes for grid elements
	//correct but just slightly off due to grid
	//grid.save("mats/test_gaussian_smooth/Outputs/grid(limit)n-c", raw_ascii);

	grid /= accu(grid); //element-wise division of sum of the column
	//good, just slightly off, barely though
	//grid.save("mats/test_gaussian_smooth/Outputs/grid_end-c", raw_ascii);
	//rowvec grid_v = conv_to<rowvec>::from(grid);

	//M: smooth_img = conv2(grid, grid, img, "same");

	vec temp;
	rowvec temp_r; // = conv_to<rowvec>::from(temp);
	mat img_build;

	// conv2(u,v,A) first convolves each column of A with the vector u,
	// and then it convolves each row of the result with the vector v.
	// VECS ARE DIFFERENT SIZE. HOW DOES MATLAB DO THIS??
	//M:img_smooth = conv2(grid, grid, img, 'same') ;
    for (uword i = 0; i < img.n_cols; i++) {
		temp = conv(grid, img.col(i));
		img_build.insert_cols(i, temp);
		if(i==0){
			/*cout<<"height temp...."<<temp.n_rows<<endl;*/
		}
	}
	/*cout<<"img_build width: "<<img_build.n_cols<<endl;
	cout<<"img_build hight: "<<img_build.n_rows<<endl;*/
	for (uword i = 0; i < img_build.n_rows; i++) {
		temp = conv(grid, img_build.row(i));
		temp_r = conv_to<rowvec>::from(temp);
		smooth_img.insert_rows(i, temp_r);
	}
	/*cout<<"smooth_img width: "<<smooth_img.n_cols<<endl;
	cout<<"smooth_img hight: "<<smooth_img.n_rows<<endl;*/


	//this done to trim the matrix down to old size of img.
	//just the math of the conv in armadillo adds some fat that
	//needs to be trimmed.
	uword grid_length = grid.n_rows;
	smooth_img = smooth_img.submat((grid_length-1)/2,(grid_length-1)/2,(smooth_img.n_rows-1) - ((grid_length)/2),(smooth_img.n_cols-1) - ((grid_length)/2));
	/*cout<<"smooth_img after trim width: "<<smooth_img.n_cols<<endl;
	cout<<"smooth_img after trim hight: "<<smooth_img.n_rows<<endl;*/

return smooth_img;
}
