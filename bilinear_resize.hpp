//This file performs bilinear resizing on the selected matrix given a scale 
//factor or a destination size. It also contains the functions for conversions
//between the Armadillo and OpenCV matrix types.

//NOTE THAT THE BILINEAR RESIZING METHOD PROVIDED BY OPENCV IS NOT EQUIVALENT
//TO THAT IN MATLAB. THIS CODE DOES NOT PRODUCE THE CORRECT RESULT.
#pragma once

#include <iostream>

#include <armadillo>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace std;
using namespace arma;

//Conversion functions
cv::Mat_<double> to_cvmat(const arma::Mat<double> &src);
arma::mat to_arma(const cv::Mat_<double> &src);

//Resizing functions
mat bilinear_resize(mat data, int scalefactor);
mat bilinear_resize(mat data, mat frame);

//This function converts data from Armadillo matrix type to OpenCV matrix type.
cv::Mat_<double> to_cvmat(const arma::Mat<double> &src) {

	int ht = src.n_rows;
	int wt = src.n_cols;

	cv::Mat dst = cv::Mat_<double>{ int(src.n_rows), int(src.n_cols), const_cast<double*>(src.memptr()) };
	cv::Mat dst_t(dst.rows, dst.cols, CV_64F);

	int i = 0, j = 0;
	//This loop rotates the image 90 degrees
	int n = 0;
	for (i; i < dst.cols; i++) {
		for (j; j < dst.rows; j++) {
			dst_t.at<double>(j, i) = dst.at<double>(n);
			n++;
		}
		j = 0;
	}
	return dst_t;
}

//This function converts data from OpenCV matrix type to Armadillo matrix type.
arma::mat to_arma(const cv::Mat_<double> &src) {
	arma::mat dst(reinterpret_cast<double*>(src.data), src.cols, src.rows);
	mat dstF(dst.n_cols, dst.n_rows);
	uword i = 0, j = 0;
	//This loop rotates the image 90 degrees
	for (i; i < dst.n_rows; i++) {
		for (j; j < dst.n_cols; j++) {
			dstF(j, i) = dst(i, j);
		}
		j = 0;
	}
	return dstF;
}

//This function uses a scale factor to resize a matrix using the OpenCV resize
//function. It must convert to an OpenCV matrix type to perform this. The 
//CV_INTER_AREA parameter was found to be as close to the MATLAB bilinear
//resize as possible.
mat bilinear_resize(mat data, double scalefactor) {

	cv::Mat src = to_cvmat(data);
	cv::Mat dst;

	cv::resize(src, dst, cv::Size(), scalefactor, scalefactor, CV_INTER_AREA);
	src.release();

	return to_arma(dst);
}

//Overloaded resize function uses a destination "frame" size instead of a scale
//factor.
mat bilinear_resize(mat data, mat frame) {

	cv::Mat temp = to_cvmat(frame);

	int ht = temp.rows;
	int wt = temp.cols;

	cv::Mat src = to_cvmat(data);
	cv::Mat dst;

	cv::resize(src, dst, temp.size(), 0, 0, CV_INTER_AREA);
	src.release();

	return to_arma(dst);
}

