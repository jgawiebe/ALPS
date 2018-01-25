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

using namespace std;
using namespace arma;
void init_variables();

<<<<<<< HEAD
    sp_mat A;//inspect matlab code, see if these are edited
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


int main() {

	//img_vert_warp.save("mats/HorizontalBar/img_vert_warpmyFile",raw_ascii);
	 /* Currently my plan is to produce these myFiles in matlab at various points
	 * in that program, starting at the end and comparing them to the corresponding
	 * outputs in the matrix_builder program also. If there are discreprencies, then
	 * the program will be cut in half and the first half of the program will be
	 * checked an so on. The notes for this portion of the lab should be very
	 * robust.  */
	init_variables();
	mat smallA;
	//get size of img_vert_warp data structure and print to screen


	tie(A,b) = build_matrix( A, b,img2_dx,img2_dy,img_z,dxx, dxy, dyy, dxz, dyz, e_data,e_smooth, u, v, gam);
	//b produces exact information. 22/jan/18 NO EDITS TO B
	//b.save("mats/test_matrix_builder/Outputs/b-c", raw_ascii);
	//A.save("mats/test_matrix_builder/Outputs/A_sp-c", raw_ascii);

	//uword height = A.n_rows;
	//uword width  = A.n_cols;
	//cout<<"height"<< height << "   width"<< width<<endl;
	//smallA = A.submat(0,0,99,99);
	//smallA.save("mats/test_matrix_builder/Outputs/smallA-c", raw_ascii);
	//create print of what A's dimensions are
=======
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


int main() {

	//img_vert_warp.save("mats/HorizontalBar/img_vert_warpmyFile",raw_ascii);
	 /* Currently my plan is to produce these myFiles in matlab at various points
	 * in that program, starting at the end and comparing them to the corresponding
	 * outputs in the matrix_builder program also. If there are discreprencies, then
	 * the program will be cut in half and the first half of the program will be
	 * checked an so on. The notes for this portion of the lab should be very
	 * robust.  */
	init_variables();

	//get size of img_vert_warp data structure and print to screen


	tie(A,b) = build_matrix( A, b,img2_dx,img2_dy,img_z,dxx, dxy, dyy, dxz, dyz, e_data,e_smooth, u, v, gam);

>>>>>>> refs/heads/jake
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

}

