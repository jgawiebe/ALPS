#include <iostream>
#include <fstream>

//#include <armadillo>
#include <opencv2/opencv.hpp>
#include <cublas.h>
#include <cusparse_v2.h>

//#include "gaussian_smooth.hpp"
//#include "gradient.hpp"
//#include "successive_overrelaxation.hpp"
//#include "energy_calc.hpp"
//#include "optical_flow.hpp"
//#include "red-black-SOR.hpp"
//#include "bilinear_resize.hpp"
#include "red_black_wrapper.hpp"

using namespace std;
//using namespace arma;
//using namespace cv;

//void gauss_test();
//void gradient_test();
void sor_test();
//void rb_sor_test();
//void psi_test();
//void energy_test();
//void resize_test();

int main() {

	//gauss_test();
	//derivative_test();
	sor_test();
	//psi_test();
	//energy_test();
	//resize_test();

	//rb_sor_test();


	//red_black_sor();
	cin.get();

	return 0;
}

//void gauss_test() {
//
//	int num_levels = 40;
//	double scale = pow(0.95, num_levels);
//
//	mat image1, image2;
//
//	image1.load("mats/g_smooth/gs_input-m.txt");
//	//image2.load("mats/main/Horizontal1.txt");
//
//	mat img1 = g_smooth(image1, scale);
////	mat img2 = g_smooth(image2, scale_factor);
//
//	img1.save("mats/g_smooth/img_smooth-c.txt", raw_ascii);
////	img2.save("mats/g_smooth/img2_smooth-c.txt", raw_ascii);
//
//}
//
//void gradient_test() {
//	mat img1_dx, img1_dy;
//	mat img2_dx, img2_dy;
//	mat img1, img2;
//
//	img1.load("mats/g_smooth/img1_smooth-m.txt");
//	img2.load("mats/g_smooth/img2_smooth-m.txt");
//
//	/*tie(img1_dx, img1_dy) = gradient(img1_dx, img1_dy, img1);*/
//
//	img1_dx.save("mats/derivatives/img1_dx-c.txt", raw_ascii);
//	img1_dy.save("mats/derivatives/img1_dy-c.txt", raw_ascii);
//}
//
//void rb_sor_test() {
//	double omega = 1.8, tolerance = 1e-8;
//	int inner_iter = 500;
//	bool fail_flag = true;
//
//	mat img_z;
//	img_z.load("mats/derivatives/img_z-m.txt");
//
//
//	vec b;
//	b.load("mats/sor/b-m.txt");
//
//
//	uword ht = img_z.n_rows;
//	uword wt = img_z.n_cols;
//	uword side_length = (ht * wt * 2);
//
//	mat A(side_length, side_length, fill::zeros);
//
//	A.load("mats/sor/du-m.txt");
//
//	mat du, dv;
//	vec duv(ht * wt * 2, fill::zeros);
//
//	/*tie(duv, fail_flag) = redblack_sor(A, duv, b, omega,
//			inner_iter, tolerance);*/
//
//	cout << "status: " << fail_flag;
//
//	duv.save("mats/sor/duv-c.txt", raw_ascii);
//}
//
//void psi_test() {
//	mat x;
//	x.load("mats/energy/psi_input.txt");
//
//	x = psi_function(x);
//
//	x.save("mats/energy/e_data-c.txt", raw_ascii);
//}
//
//void energy_test() {
//	mat u, v, e_smooth;
//	//vec du, dv; //zero for test 1
//	u.load("mats/energy/u-m.txt");
//	v.load("mats/energy/v-m.txt");
//
//	e_smooth = generate_esmooth(u, v);
//
//	e_smooth.save("mats/energy/e_smooth-c.txt", raw_ascii);
//}
//
//void resize_test() {
//	mat u, img_size;
//
//	u.load("mats/resize/u_in.txt");
//	img_size.load("mats/resize/img_size.txt");
//
//	u = bilinear_resize(u, img_size);
//
//	/*cv::Mat cv_u = to_cvmat(u);
//	cv::Mat cv_img = to_cvmat(img_size);
//
//	resize(cv_u, cv_u, cv_img.size());
//
//	arma::mat u_rx = to_arma(cv_u);*/
//
//
//
//	//cv::Mat cv_u(cv_u.rows, cv_u.cols, CV_64FC1, u.memptr());
//	//cv::Mat cv_img(cv_img.rows, cv_img.cols, CV_64FC1, img1.memptr());
//
//	//cv::Mat_<double>{int(u.n_cols), int(u.n_rows), const_cast<double*>(u.memptr())};
//
//
//	
//
//	//cv::Mat cv_img = to_cvmat(img_size);
//
//	//int cols = img1.n_cols;
//	//int rows = img1.n_rows;
//
//	//cv::FileStorage storage("test.yml", cv::FileStorage::WRITE);
//	//storage << "img" << img;
//	//storage.release();
//
//	//resize(cv_img, cv_u, cv_u.size());
//
//	//arma::mat img_rx(reinterpret_cast<double*>(cv_img.data), cv_img.cols, cv_img.rows);
//	//arma::mat u_rx(reinterpret_cast<double*>(cv_u.data), cv_u.cols, cv_u.rows);
//
//	
//
//	/*cout << "input col >" << u.n_cols << " row >" << u.n_rows << endl;
//	cout << "img col >" << img_size.n_cols << " row >" << img_size.n_rows << endl;
//	cout << "output col >" << u_rx.n_cols << " row >" << u_rx.n_rows << endl;
//*/
//	u.save("mats/resize/u_c.txt", raw_ascii);
//	//v.save("mats/main/v-c.txt", raw_ascii);
//}