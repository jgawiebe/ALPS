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

	//**get images using opencv**
	mat img_vert_norm;
	mat img_vert_warp;
    //Know size of horizontal matrix is 607x684

	img_vert_norm.load("mats/HorizontalBar/hori_norm.mat");
//	image2.load("mats/main/img2-m.txt");
//
//	image1.save("mats/main/img1-c.txt", csv_ascii);

	//get size of image
//	uword height = image1.n_rows;
//	uword width = image1.n_cols;

//	cout << "height: " << height << endl;
//	cout << "width: " << width << endl;

//
//	img1.save("out/img1.txt", arma_ascii);
//	img2.save("out/img2.txt", arma_ascii);

	return 0;
}

//void simple_readwrite(mat input) {
//	input.load("mats/img1.mat", raw_ascii);
//
//	string buf;
//	uword i = 0;
//
//	ifstream in_mat("in/main_img1.txt");
//	ofstream out_mat("test_output.txt");
//
//	while (getline(in_mat, buf)) {
//		if (i < input.n_cols) {
//			out_mat << buf;
//		} else {
//			out_mat << endl << buf;
//		}
//	}
//	out_mat.close();
//}
