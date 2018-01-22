#include <iostream>
#include <fstream>

#include <armadillo>

#include "gaussian_smooth.hpp"
#include "gradient.hpp"
#include "successive_overrelaxation.hpp"

using namespace std;
using namespace arma;

void gauss_test();
void derivative_test();
void sor_test();

int main() {

	gauss_test();

	//derivative_test();

	//sor_test();

	return 0;
}

void gauss_test() {

	int num_levels = 40;
	double scale_factor = pow(0.95, num_levels);

	mat image1, image2;

	image1.load("mats/main/img1-m.txt");
	image2.load("mats/main/img2-m.txt");

	mat img1 = g_smooth(image1, scale_factor);
//	mat img2 = g_smooth(image2, scale_factor);

	img1.save("mats/g_smooth/img1_smooth-c.txt", raw_ascii);
//	img2.save("mats/g_smooth/img2_smooth-c.txt", raw_ascii);

}

void derivative_test() {
	mat img1_dx, img1_dy;
	mat img2_dx, img2_dy;
	mat img1, img2;

	img1.load("mats/g_smooth/img1_smooth-m.txt");
	img2.load("mats/g_smooth/img2_smooth-m.txt");

	tie(img1_dx, img1_dy) = gradient(img1_dx, img1_dy, img1);

	img1_dx.save("mats/derivatives/img1_dx-c.txt", raw_ascii);
	img1_dy.save("mats/derivatives/img1_dy-c.txt", raw_ascii);
}

void sor_test() {
	double omega = 1.8, tolerance = 1e-8;
	int inner_iter = 500;
	uword fail_flag = 0;

	mat img_z;
	img_z.load("mats/derivatives/img_z-m.txt");

	vec b;
	b.load("mats/sor/b-m.txt");

	uword ht = img_z.n_rows;
	uword wt = img_z.n_cols;
	uword side_length = (ht * wt * 2);

	sp_mat A(side_length, side_length);

	//A.load("mats/sor/A-m.txt");


	mat du, dv;
	vec duv(ht * wt * 2, fill::zeros);


	duv = successive_overrelaxation(&fail_flag, A, duv, b, omega, inner_iter,
			tolerance);

	duv.save("mats/sor/duv-c.txt", raw_ascii);
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
