#include <iostream>
#include <fstream>

#include <armadillo>

#include "gaussian_smooth.hpp"
#include "gradient.hpp"
#include "successive_overrelaxation.hpp"
#include "energy_calc.hpp"
#include "red-black-SOR-v2.hpp"
using namespace std;
using namespace arma;

void gauss_test();
void gradient_test();
//void sor_test();
void rb_sor_test();
void psi_test();
void energy_test();
void resize_test();

int main() {

	//gauss_test();
	//derivative_test();
	//sor_test();
	//psi_test();
	//energy_test();
	resize_test();

	//rb_sor_test();

	return 0;
}

void gauss_test() {

	int num_levels = 40;
	double scale = pow(0.95, num_levels);

	mat image1, image2;

	image1.load("mats/g_smooth/gs_input-m.txt");
	//image2.load("mats/main/Horizontal1.txt");

	mat img1 = g_smooth(image1, scale);
//	mat img2 = g_smooth(image2, scale_factor);

	img1.save("mats/g_smooth/img_smooth-c.txt", raw_ascii);
//	img2.save("mats/g_smooth/img2_smooth-c.txt", raw_ascii);

}

void gradient_test() {
	mat img1_dx, img1_dy;
	mat img2_dx, img2_dy;
	mat img1, img2;

	img1.load("mats/g_smooth/img1_smooth-m.txt");
	img2.load("mats/g_smooth/img2_smooth-m.txt");

	tie(img1_dx, img1_dy) = gradient(img1_dx, img1_dy, img1);

	img1_dx.save("mats/derivatives/img1_dx-c.txt", raw_ascii);
	img1_dy.save("mats/derivatives/img1_dy-c.txt", raw_ascii);
}

void rb_sor_test() {
	double omega = 1.8, tolerance = 1e-8;
	int inner_iter = 500;
	bool fail_flag = true;

	mat img_z;
	img_z.load("mats/derivatives/img_z-m.txt");


	vec b;
	b.load("mats/sor/b-m.txt");


	uword ht = img_z.n_rows;
	uword wt = img_z.n_cols;
	uword side_length = (ht * wt * 2);

	mat A(side_length, side_length, fill::zeros);

	//A.load("mats/sor/du-m.txt");

	mat du, dv;
	vec duv(ht * wt * 2, fill::zeros);

	tie(duv, fail_flag) = redblack_sor(A, duv, b, omega,
			inner_iter, tolerance);

	cout << "status: " << fail_flag;

	duv.save("mats/sor/duv-c.txt", raw_ascii);
}

void psi_test() {
	mat x;
	x.load("mats/energy/psi_input.txt");

	x = psi_function(x);

	x.save("mats/energy/e_data-c.txt", raw_ascii);
}

void energy_test() {
	mat u, v, e_smooth;
	//vec du, dv; //zero for test 1
	u.load("mats/energy/u-m.txt");
	v.load("mats/energy/v-m.txt");

	e_smooth = generate_esmooth(u, v);

	e_smooth.save("mats/energy/e_smooth-c.txt", raw_ascii);
}

void resize_test() {
	mat u, v, img1;

	img1.load("mats/main/img1-m.txt");
	u.load("mats/main/u-in-m.txt");
	v.load("mats/main/v-in-m.txt");

	u.resize(img1.n_rows, img1.n_cols);
	v.resize(size(u));

	u.save("mats/main/u-c.txt", raw_ascii);
	v.save("mats/main/v-c.txt", raw_ascii);
}

//void sor_test() {
//	double omega = 1.8, tolerance = 1e-8;
//	int inner_iter = 500;
//	uword fail_flag = 0;
//
//	mat img_z;
//	img_z.load("mats/derivatives/img_z-m.txt");
//
//	vec b;
//	b.load("mats/sor/b-m.txt");
//
//	uword ht = img_z.n_rows;
//	uword wt = img_z.n_cols;
//	uword side_length = (ht * wt * 2);
//
//	mat A(side_length, side_length);
//
//	A.load("mats/sor/A-m.txt");
//
//	mat du, dv;
//	vec duv(ht * wt * 2, fill::zeros);
//
//	duv = successive_overrelaxation(A, fail_flag, duv, b, omega, inner_iter,
//			tolerance);
//
//	duv.save("mats/sor/duv-c.txt", raw_ascii);
//}
