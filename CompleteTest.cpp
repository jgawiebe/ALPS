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

void cli() {
	cout << "Welcome to PSEAG: Would you like to specify parameters(Y)? Otherwise press enter.";
}

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
		//image1.save("FULLTEST-img1.txt", raw_ascii);
		//image2.save("FULLTEST-img2.txt", raw_ascii);
		//cin.get(); //manual break
	}

	//perform guassian scaling on images
	mat img1(image1.n_rows, image1.n_cols);
	mat img2(image2.n_rows, image2.n_cols);

	img1 = g_smooth(image1, scale_factor);
	img2 = g_smooth(image2, scale_factor);
	//EXACT MATCH OF MATLAB
	//img1.save("FULLTEST-1smooth.txt", raw_ascii);
	//img2.save("FULLTEST-2smooth.txt", raw_ascii);


	img1 = bilinear_resize(img1, scale_factor);///////////////////////////////////////////
	img2 = bilinear_resize(img2, scale_factor);///////////////////////////////////////////

	//IMAGES NEEDED TO BE ROTATED 90 DEGREES, 
	//This was done. However the scaling is not perfecto.
	//Done with Bilin also, not good either
	//img1.save("FULLTEST-1resizeAREA.txt", raw_ascii);
	//img2.save("FULLTEST-2resizeAREA.txt", raw_ascii);
	///////////////////////////////////INJECT OF MATLAB DATA TO SYNCH THE EXECUTION DATA
	//img1.load("I1.txt");
	//img2.load("I2.txt");
	//img1.save("FULLTEST-I1inject.txt", raw_ascii);
	//img2.save("FULLTEST-I2inject.txt", raw_ascii);
	//NOW I HAVE GOOD DATA FOR FIRST OPTICAL FLOW ATTACK
	//cin.get();
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
		//OPTICAL FLOW PASS FOR A SINGLE LOOP THROUGH, ERROR IS DUDE TO 
		//THE FACT THAT WE DON'T ERROR CHECK RIGHT NOW SO WE DO 500 
		//ITERATIONS NO MATTER WHAT
		//du.save("FULLTEST-du_after_flow.txt", raw_ascii);
		//dv.save("FULLTEST-dv_after_flow.txt", raw_ascii);

		
		
		//add incremental change in x and y domain
		u = u + du;
		v = v + dv;

		//smooth and scale original images to next level of the pyramid
		img1 = g_smooth(image1, scale_factor);
		img2 = g_smooth(image2, scale_factor);

		
		//cout << "height:" << img1.n_rows << endl;
		//cout << "width:" << img2.n_cols << endl;
		//cout << scale_factor << endl;
		//cout << img1.n_rows*scale_factor << endl;
		//cout << img2.n_cols*scale_factor << endl;

		//cout << "addr1" << img2.memptr() << endl;
		img2 = bilinear_resize(img2, scale_factor);
		//cout << "addr2" << img2.memptr() << endl;
		//cout << scale_factor;
		
		img1 = bilinear_resize(img1, scale_factor);
		
		//resize to current resolution
		u = bilinear_resize(u, img1);
		v = bilinear_resize(v, img1);

		//cout << "u mat:" << endl;
		//cout << u << endl;

		//cout << "img1:" << endl;
		//cout << img1 << endl;
		//cin.get();

		

		//cout << img1.n_rows << 'x' << img1.n_cols << "  img dim" << endl;
		//cout << u.n_rows << 'x' << u.n_cols << "  u dim" << endl;
		//cout << i << endl;
		//cout << scale_factor << endl;

		//cin.get();


		//compare with u-m, v-m
		//cout << u << "\r";
		u.save("intermediate-u.txt", raw_ascii);
		v.save("intermediate-v.txt", raw_ascii);

		//img1.save("FULLTEST-I1", raw_ascii);
		//img2.save("FULLTEST-I2", raw_ascii);

		cout << "--------" << i - 1 << " of " << num_levels << " levels of guassian pyramid remaining--------\n";
		//cout << "Hit enter to continue\n\n" << endl;
		//cin.get();
		//Display image
		/*namedWindow("Image 1");
		imshow("Image 1", to_cvmat(img1));*/
	}
	cout << "ALGORTIHM COMPLETE: SUCCESS\nPress any key to exit" << endl;


	//not part of tesing below
	u.save("u-jacobi.dat", raw_ascii);
	v.save("v-jacobi.dat", raw_ascii);


	cin.get();
	exit(EXIT_SUCCESS); 
	return 0;
}

