/*
gradient.hpp
Jacob Wiebe & James Dolman
Rev1: Nov 2017
*/

#include <iostream>
#include <fstream>

#include <armadillo>

//#include <cv>

using namespace std;
using namespace arma;
//using namespace cv;
//g++ test.cpp -o test -O2 -larmadillo `pkg-config --cflags --libs opencv`; ./test

//compile with: g++ test.cpp -o test -O2 -larmadillo; ./test
int main(){

	ofstream out_mat;
	out_mat.open ("test_output.txt");
	out_mat << "Writing this to a file.\n";
	out_mat.close();


  //TESTING ARMA COMMANDS
//  mat A(54, 36, fill::randu);
//
//  SizeMat my_size = size(A);
//  cout << "Size is: " << my_size << endl;
//
//  int height = A.n_rows;
//	int width = A.n_cols;
//
//  cout << "Size is: " << height << 'x' << width << endl;


//TESTING MATRIX_BUILDER
  // mat u(50, 50, fill::randu);
  // mat e_smooth(8, 16, fill::randu);



 // cv::Mat image;
 // image = imread("Car.bmp", CV_LOAD_IMAGE_UNCHANGED);


  // mat B, C;
  // B.ones(50);
  // C.ones(50);

  //gradient(A, B, C);


  return 0;
}
