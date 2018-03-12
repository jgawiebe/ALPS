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
	frowvec row, col;
	vec val;
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
	//bool* failure;

	unsigned int n, x_size;
	const unsigned int max_iter = 500;

	double* parallel_sor(double *val, float *row, float *col, double *b, const unsigned int n, const unsigned int x_size, const unsigned int max_iter, double tol);
	extern double* serial_sor(double *val, float *row, float *col, double *b, const unsigned int n, const unsigned int x_size, const unsigned int max_iter);


int main() {
	//testing on matrix builder with new matlab done again, all good!
	//Creat new input vectors for test_SOR Then carry on

	init_variables();

	tie(row, col, val, b) = build_matrix(img2_dx, img2_dy, img_z, dxx, dxy, dyy, dxz, dyz, e_data, e_smooth, u, v, gam);
	cout << "Out of build_matrix" << endl;

	n = (int)val.n_elem;
	x_size = (int)b.n_elem;

	cout << n << endl;
	//cout << row << endl;
	//tie(A, duv) = build_sparse(img2_dx, img2_dy, img_z, dxx, dxy, dyy, dxz, dyz, e_data, e_smooth, u, v, gam);
	//cout << "Out of build_sparse" << endl;
	//cin.get();

	//double* x = serial_sor(val.memptr(), row.memptr(), col.memptr(), b.memptr(), n, x_size, max_iter);
	double* x = parallel_sor(val.memptr(), row.memptr(), col.memptr(), b.memptr(), n, x_size, max_iter, tol);
	cout << "Out of cuda - hit enter to save duv" << endl;
	//cout << "Fail: " << *failure << endl;
	
	vec duv(x, x_size); //may want to try the other parameters
	cin.get();
	duv.save("duv_cuda-c.txt", raw_ascii);


	//tie(A, duv) = build_sparse(img2_dx, img2_dy, img_z, dxx, dxy, dyy, dxz, dyz, e_data, e_smooth  , u, v, gam);
	//cout<<"Out of build_sparse"<<endl;
	////cin.get();

	//tie(duv, failure) = successive_overrelaxation(A, duv, omega, max_it, tol);
 //   cout<<"Out of successive_overrelaxation in test main() - hit enter to save duv"<<endl;
	//cout << "Fail: "<<failure; 
	//cin.get();
	//duv.save("duv_serial-c.txt", raw_ascii);

	

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
	max_it = 500;
	tol = 1e-8;
	//*failure = true;

}
