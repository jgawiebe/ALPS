/*
 * help.cpp
 *
 *  Created on: Jan 18, 2018
 *      Author: s27508
 */

#include <iostream>
#include <fstream>

#include <armadillo>
#include "matrix_builder.hpp"
#include "successive_overrelaxation.hpp"

using namespace std;
using namespace arma;
void init_variables();
	//variables for matrix_builder
    mat A;//inspect matlab code, see if these are edited
	vec b;
	mat img2_dx;
	mat img2_dy;
	mat img_z;
	mat dxx; // was experiencing data loss with umat, not with mat
	mat dxy;
	mat dyy;
	mat dxz;
	mat dyz;
    mat e_data;
	mat e_smooth;
	mat u;
	mat v;
	double gam;

	//variables for successive_overrelaxation (less A and b)
	double omega;
	double max_it;
	double tol;
	vec duv;
	uword* failure;

int main() {

	init_variables();

	tie(A,b) = build_matrix( A, b,img2_dx,img2_dy,img_z,dxx, dxy, dyy, dxz, dyz, e_data,e_smooth, u, v, gam);

	//build_matrix produces good info for A and b
	//last 3 columns of A are 0s and will be trimmed here.
	A.shed_col(3);
	A.shed_col(3);
	A.shed_col(3);
	//0's columns on A trimmed, ready for input to SOR

	duv = successive_overrelaxation(failure,A,b,omega, max_it,tol );

	duv.save("mats/test_SOR/Outputs/duv-c", raw_ascii);

	//smallA = A.submat(0,0,99,99);
	//smallA.save("mats/test_matrix_builder/Outputs/smallA-c", raw_ascii);
	//create print of what A's dimensions are
	return 0;
}


    /* Know size of horizontal matrix is 607x684 which is the correct resolution
	 * for the image.issue here, on inspection of the img_vert_warp matrix it was found that
	 * it was not the 607x684 matrix it should have been - just like it is in
	 * matlab. Instead it was a 2965x1 column vector.
	 *
	 * IF you simply .load the actual bmp file then what is produced is a
	 * 155,702x1 column vector. This is better but the 607x684 is still far
	 * off as that would be a 415,188 row column vector.
	 *
	 * Note the myFile.txt method works for bringing in the correct values and
	 * creating a matrix structure, img_ver_warp has dimensions 607x684 which
	 * is correct. This myFile.txt was created in matlab using the command:
	 * M: c0 = rgb2gray(imread('Horizontal1.bmp'));
	 * M: dlmwrite('myFile.txt',c0).
	 *
	 *
	 * umat img_vert_warp;
	 * img_vert_warp.load("mats/HorizontalBar/myFile.txt");
	 */

void init_variables(){
	//variables loaded in for matrix_builder
	img2_dx.load("mats/test_matrix_builder/Inputs/Ikx.txt");
	img2_dy.load("mats/test_matrix_builder/Inputs/Iky.txt");
	img_z.load("mats/test_matrix_builder/Inputs/Ikz.txt");
	dxx.load("mats/test_matrix_builder/Inputs/Ixx.txt");
	dxy.load("mats/test_matrix_builder/Inputs/Ixy.txt");
	dyy.load("mats/test_matrix_builder/Inputs/Iyy.txt");
	dxz.load("mats/test_matrix_builder/Inputs/Ixz.txt");
    dyz.load("mats/test_matrix_builder/Inputs/Iyz.txt");
	e_data.load("mats/test_matrix_builder/Inputs/E_Data.txt");
	e_smooth.load("mats/test_matrix_builder/Inputs/aE_smooth.txt");
	u.load("mats/test_matrix_builder/Inputs/u.txt");
	v.load("mats/test_matrix_builder/Inputs/v.txt");
    gam = 80; //value stored in gamma.txt

    //variables loaded in for successive_overrelaxation.
    omega = 80;
    max_it = 500;
	duv.load("mats/test_SOR/Inputs/duv.txt");
	tol = 1*(10^(-8));
	*failure = 1;

}

