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

	//variables for successive_overrelaxation (less A and b)
	double omega;
	double max_it;
	double tol;
	vec duv;
	uword failure;


int main() {
	//testing on matrix builder with new matlab done again, all good!
	//Creat new input vectors for test_SOR Then carry on

	init_variables();
	cout << "CHECK IMAGE READINGS" << endl;

	tie(A, duv) = build_matrix(img2_dx, img2_dy, img_z, dxx, dxy, dyy, dxz, dyz, e_data, (10.0 * e_smooth), u, v, gam);
	 cout<<"Out of build_matrix in test main()"<<endl;
	 cin.get();
	//build_matrix produces good info for A and b

	//0's columns on A trimmed, ready for input to SOR
	//note A is a dense matrix at this point of 3 columns which needs to be converted into a
	//square sparse matrix.

	 tie(duv, failure) = successive_overrelaxation(A, duv, omega, max_it, tol);
    cout<<"Out of successive_overrelaxation in test main()"<<endl;
	cin.get();

	duv.save("duv-c.txt", raw_ascii);
	//b.save("mats/test_matrix_builder/OutputsV2/bv2-c", raw_ascii);
	//A.save("mats/test_matrix_builder/OutputsV2/Av2-c", raw_ascii);

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
	img2_dx.load("mats/SmallHori/IkxHoriSmall.txt");
	img2_dy.load("mats/SmallHori/IkyHoriSmall.txt");
	img_z.load("mats/SmallHori/IkzHoriSmall.txt");
	dxx.load("mats/SmallHori/IxxHoriSmall.txt");
	dxy.load("mats/SmallHori/IxyHoriSmall.txt");
	dyy.load("mats/SmallHori/IyyHoriSmall.txt");
	dxz.load("mats/SmallHori/IxzHoriSmall.txt");
    dyz.load("mats/SmallHori/IyzHoriSmall.txt");
	e_data.load("mats/SmallHori/E_DataHoriSmall.txt");
	e_smooth.load("mats/SmallHori/aE_smoothHoriSmall.txt");
	u.load("mats/SmallHori/uHoriSmall.txt");
	v.load("mats/SmallHori/vHoriSmall.txt");
    gam = 80; //value stored in gamma.txt

    //variables loaded in for successive_overrelaxation.
    omega = 1.8;
    //max_it = 500;
	//tol = 1*(10^(-8));
	max_it = 50;
	tol = 1e-8;
	failure = 0;

}
