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

int main() {

	mat du, dv;

	//**get images using opencv that is the ideal thinggg **
	mat img_vert_norm;
	umat img_vert_warp;

    /*Know size of horizontal matrix is 607x684 which is the correct resolution
	for the image.*/
	img_vert_warp.load("mats/HorizontalBar/myFile.txt");

	img_vert_warp.save("mats/HorizontalBar/img_vert_warpmyFile",raw_ascii);
	/*issue here, on inspection of the img_vert_warp matrix it was found that
	 * it was not the 607x684 matrix it should have been - just like it is in
	 * matlab. Instead it was a 2965x1 column vector.
	 *
	 * IF you simply .load the actual bmp file then what is produced is a
	 * 155,702x1 column vector. This is better but the 607x684 is still far
	 * off as that would be a 415,188 row column vector.
	 *
	 * Note the myFile.txt method works at least for dimensions and creating a
	 * matrix structure, img_ver_warp has dimensions 607x684 which is correct.
	 * This myFile.txt was created in matlab using the command:
	 * M: c0 = rgb2gray(imread('Horizontal1.bmp'));
	 * M: dlmwrite('myFile.txt',c0) */

	//get size of img_vert_warp data structure and print to screen
	uword height = img_vert_warp.n_rows;
	uword width = img_vert_warp.n_cols;
	cout << "height of img_vert_warp: " << height << "  ";
	cout << "width img_vert_warp: " << width << endl;



	return 0;
}


