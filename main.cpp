/*
 energy_calc.hpp
 Jacob Wiebe & James Dolman
 Rev1: Nov 2017
 */

#include <iostream>
#include <armadillo>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "bilinear_resize.hpp"
#include "optical_flow.hpp"
#include "gaussian_smooth.hpp"

using namespace std;
using namespace arma;
using namespace cv;

//compile with: g++ <funtion>.cpp -o test -O2 -larmadillo; ./test

//M: optic_flow_brox
int main() {
	//set constants
	const double num_levels = 10.0, alpha = 10.0, gamma = 80.0, omega = 1.8;
	const int outer_iter = 1, inner_iter = 50;
	
	//set intial scale factor
	double scale_factor = pow(0.95, num_levels);
	
	//initialize input matrices
	mat image1 = to_arma(imread("cars1.png", CV_LOAD_IMAGE_GRAYSCALE));
	mat image2 = to_arma(imread("cars2.png", CV_LOAD_IMAGE_GRAYSCALE));

	//FOR TESTING MEMORY ERROR
	image1 = image1.submat(0, 0, 60, 60);
	image2 = image2.submat(0, 0, 60, 60);



	if (image1.n_elem == 0 || image2.n_elem == 0) {
		cout << "error: images could not be opened" << endl;
		system("pause"); //wait for any key press
		return -1;
	}
	else {
		cout << "Images loaded successfully" << endl;

		//DON'T USE THIS WHEN BENCHMARKING
		//String window_name1 = "Initial Image 1";
		//String window_name2 = "Initial Image 2";

		//namedWindow(window_name1);
		//namedWindow(window_name2);

		//cv::imshow(window_name1, to_cvmat(image1));
		//cv::imshow(window_name2, to_cvmat(image2));
	}

	//perform guassian scaling on images
	mat img1 = g_smooth(image1, scale_factor);
	mat img2 = g_smooth(image2, scale_factor);

	img1 = bilinear_resize(img1, img1*scale_factor);
	img2 = bilinear_resize(img2, img2*scale_factor);

	//define u and v matrices
	mat u(img1.n_rows, img1.n_cols, fill::zeros);
	mat v(img1.n_rows, img1.n_cols, fill::zeros);

	//define du and dv matrices
	mat du, dv;

	//levels of the guassian pyramid
	for (int i = 0; i < num_levels; i++) {
		
		//resolution increases with each pyramid level
		scale_factor = pow(0.95, i);
		
		//derivatives are stored as optical_flow globals
		compute_derivatives(img1, img2);
		
		//perform optical_flow to get du and dv
		tie(u, v, du, dv) = optical_flow(alpha, gamma, omega, u, v, outer_iter,
				inner_iter);

		//add incremental change in x and y domain
		u = u + du;
		v = v + dv;
		
		//scale original images to next level of the pyramid
 		img1 = g_smooth(image1, scale_factor);
		img2 = g_smooth(image2, scale_factor);

		//resize to current resolution
		u = bilinear_resize(u, img1);
		v = bilinear_resize(v, img1);


		//compare with u-m, v-m
		u.save("mats/fin/u-c.txt", raw_ascii);
		v.save("mats/fin/v-c.txt", raw_ascii);

		img1.save("mats/main/img1-output.txt", raw_ascii);
		img2.save("mats/main/img2-output.txt", raw_ascii);

		cout << "--------" << i + 1 << "/" << num_levels << " levels of guassian pyramid complete--------\n\n";
	}
	cout << "ALGORTIHM COMPLETE: SUCCESS\nPress any key to exit" << endl;
	cin.get();
	exit(EXIT_SUCCESS);
}

