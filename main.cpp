#include <iostream>
#include <armadillo>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "nvblas.h"

#include "bilinear_resize.hpp"
#include "optical_flow.hpp"
#include "gaussian_smooth.hpp"

using namespace std;
using namespace arma;
using namespace cv;

void cli() {
	cout << "Welcome to PSEAG: Would you like to specify parameters(Y)? Otherwise press enter.";
	
}

int main() {
	//set constants
	const double num_levels = 40.0, alpha = 10.0, gamma = 80.0, omega = 1.8, sor_tolerance = 1e-8;
	const int outer_iter = 3, inner_iter = 500;
	
	//set intial scale factor
	double scale_factor = pow(0.95, num_levels);

	//initialize input matrices
	mat image1 = to_arma(imread("cars1.png", CV_LOAD_IMAGE_GRAYSCALE));
	//mat image2 = to_arma(imread("cars2.png", CV_LOAD_IMAGE_GRAYSCALE));

	mat image1, image2;
		
	image1.load("cars1.txt");
	image2.load("cars2.txt");

	//FOR TESTING MEMORY ERROR
	image1 = image1.submat(0, 0, 9, 9);
	image2 = image2.submat(0, 0, 9, 9);

	//cout << image1 << endl;



	if (image1.n_elem == 0 || image2.n_elem == 0) {
		cout << "error: images could not be opened" << endl;
		system("pause"); //wait for any key press
		return -1;
	}
	else {
		cout << "Images loaded successfully" << endl;
		image1.save("img1-in.txt", raw_ascii);
		image2.save("img2-in.txt", raw_ascii);
		//cin.get(); //manual break
	}
	

	//perform guassian scaling on images
	mat img1 = g_smooth(image1, scale_factor);
	mat img2 = g_smooth(image2, scale_factor);

	//int scale_height = (int)img1.n_cols*scale_factor;

	img1 = bilinear_resize(img1, (int) (img1.n_cols*scale_factor + 1), (int) (img1.n_rows*scale_factor + 1));
	img2 = bilinear_resize(img2, (int) (img2.n_cols*scale_factor + 1), (int) (img2.n_rows*scale_factor + 1));


	img1.save("1resize.txt", raw_ascii);
	img1.save("2resize.txt", raw_ascii);

	//define u and v matrices
	mat u(img1.n_rows, img1.n_cols, fill::zeros);
	mat v(img1.n_rows, img1.n_cols, fill::zeros);

	//define du and dv matrices
	mat du, dv;

	//levels of the guassian pyramid
	for (int i = num_levels; 0 < i; i--) {

		//resolution increases with each pyramid level
		scale_factor = pow(0.95, i);
		
		//img1.save("img1-gaus.txt", raw_ascii);
		//img2.save("img2-gaus.txt", raw_ascii);
		//cin.get(); //manual break

		//derivatives are stored as optical_flow globals
		compute_derivatives(img1, img2);
		
		//perform optical_flow to get du and dv
		tie(u, v, du, dv) = optical_flow(alpha, gamma, omega, sor_tolerance, u, v, outer_iter,
				inner_iter);

		//add incremental change in x and y domain
		u = u + du;
		v = v + dv;
		
		//smooth and scale original images to next level of the pyramid
 		img1 = g_smooth(image1, scale_factor);
		img2 = g_smooth(image2, scale_factor);

		img1 = bilinear_resize(img1, (int)(img1.n_cols*scale_factor + 1), (int)(img1.n_rows*scale_factor + 1));
		img2 = bilinear_resize(img2, (int)(img2.n_cols*scale_factor + 1), (int)(img2.n_rows*scale_factor + 1));

		//resize to current resolution
		u = bilinear_resize(u, img1.n_cols, img1.n_rows);
		v = bilinear_resize(v, img1.n_cols, img1.n_rows);


		//compare with u-m, v-m
		/*u.save("mats/fin/u-c.txt", raw_ascii);
		v.save("mats/fin/v-c.txt", raw_ascii);

		img1.save("mats/fin/img1-output.txt", raw_ascii);
		img2.save("mats/fin/img2-output.txt", raw_ascii);*/

		cout << "--------" << i - 1 << " of " << num_levels << " levels of guassian pyramid remaining--------\n";
		cout << "Hit enter to continue\n\n" << endl;

		//Display image
		/*namedWindow("Image 1");
		imshow("Image 1", to_cvmat(img1));*/
	}
	cout << "ALGORTIHM COMPLETE: SUCCESS\nPress any key to exit" << endl;

	
	u.save("u-c.txt", raw_ascii);
	v.save("v-c.txt", raw_ascii);
	//cin.get();
	exit(EXIT_SUCCESS);
}

