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

int main(void) {
	//set constants
	const double num_levels = 15.0, alpha = 30.0, gamma = 80.0, omega = 1.8, sor_tolerance = 1e-8;
	const int outer_iter = 3, inner_iter = 500;
	const bool parallel_mode = true;
	
	//set intial scale factor
	double scale_factor = pow(0.95, num_levels);

	mat image1, image2;

	cout << "ALPS implementation of Brox Algorithm - Jacobi version" << endl;
	cin.get();

	cv::Mat c1 = imread("ch0.png");
	cv::cvtColor(c1, c1, cv::COLOR_BGR2GRAY);
	cv::Mat c2 = imread("ch1.png");
	cv::cvtColor(c2, c2, cv::COLOR_BGR2GRAY);
	
	image1 = to_arma(c1);
	image2 = to_arma(c2);
	
	//PRINT IMAGE 1 AND 2 AND VERIFY LOAD IN CORRECT 
	//LOAD IN PASS, GREY SCALE VALUES THE SAME  

	if (image1.n_elem == 0 || image2.n_elem == 0) {
		cout << "error: images could not be opened" << endl;
		system("pause"); //wait for any key press
		return -1;
	}else {
		cout << "Images loaded successfully" << endl;
	}

	//perform guassian scaling on images
	mat img1(image1.n_rows, image1.n_cols);
	mat img2(image2.n_rows, image2.n_cols);

	img1 = g_smooth(image1, scale_factor);
	img2 = g_smooth(image2, scale_factor);
	//EXACT MATCH OF MATLAB

	img1 = bilinear_resize(img1, scale_factor);
	img2 = bilinear_resize(img2, scale_factor);

	//IMAGES NEEDED TO BE ROTATED 90 DEGREES, 
	//NOW I HAVE GOOD DATA FOR FIRST OPTICAL FLOW ATTACK
	//define u and v matrices
	mat u(img1.n_rows, img1.n_cols, fill::zeros);
	mat v(img1.n_rows, img1.n_cols, fill::zeros);

	//define du and dv matrices
	mat du, dv;

	for (int i = num_levels; 0 < i; i--) {

		//resolution increases with each pyramid level
		scale_factor = pow(0.95, i-1);
		
		//derivatives are stored as optical_flow globals
		compute_derivatives(img1, img2);
		//COMPUTE DERIVATIVES PASSES THE CHAIN
		//perform optical_flow to get du and dv
		
		tie(u, v, du, dv) = optical_flow(alpha, gamma, omega, sor_tolerance, u, v, outer_iter,
			inner_iter, parallel_mode);
		//OPTICAL FLOW PASS FOR A SINGLE LOOP THROUGH, ERROR IS DUE TO 
		//THE FACT THAT WE DON'T ERROR CHECK RIGHT NOW SO WE DO 500 
		//ITERATIONS NO MATTER WHAT
		
		//add incremental change in x and y domain
		u = u + du;
		v = v + dv;

		//smooth and scale original images to next level of the pyramid
		img1 = g_smooth(image1, scale_factor);
		img2 = g_smooth(image2, scale_factor);
		img2 = bilinear_resize(img2, scale_factor);
		img1 = bilinear_resize(img1, scale_factor);
		
		//resize to current resolution
		u = bilinear_resize(u, img1);
		v = bilinear_resize(v, img1);
		//print statement for program execution progress
		cout << "--------" << i - 1 << " of " << num_levels << " levels of guassian pyramid remaining--------\n";
		
	}
	
	cout << "ALGORTIHM COMPLETE: SUCCESS\nPress any key to exit" << endl;

	//not part of tesing below
	u.save("u-jacobi.dat", raw_ascii);
	v.save("v-jacobi.dat", raw_ascii);

	cin.get();
	cout << "Hit Enter to Exit";
	exit(EXIT_SUCCESS); 
	return 0;
}

