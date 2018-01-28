/*
 * help.cpp
 *
 *  Created on: Jan 18, 2018
 *      Author: s27508
 */

#include <iostream>
#include <fstream>

#include <armadillo>

#include "gaussian_smooth.hpp"

using namespace std;
using namespace arma;

int main() {
//	double alpha = 30.0, gamma = 80.0, omega = 1.8; //check omega value
	int num_levels = 40, outer_iter = 3, inner_iter = 500;

	double scale_factor = pow(0.95, num_levels);

	mat du, dv;

	//**get images using opencv**
	mat image1, image2;
//	image1.load("mats/main/img1-m.txt");
//	image2.load("mats/main/img2-m.txt");
//
//	image1.save("mats/main/img1-c.txt", csv_ascii);

	//get size of image
	uword height = image1.n_rows;
	uword width = image1.n_cols;

	cout << "height: " << height << endl;
	cout << "width: " << width << endl;

	//perform gaussian scaling on images
	mat img1 = g_smooth(image1, scale_factor);
//	mat img2 = g_smooth(image2, scale_factor);
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
