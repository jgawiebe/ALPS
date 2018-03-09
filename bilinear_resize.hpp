#pragma once

#include <iostream>

#include <armadillo>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace arma;
//using namespace cv;

cv::Mat_<double> to_cvmat(const arma::Mat<double> &src);
arma::mat to_arma(const cv::Mat_<double> &src);
mat bilinear_resize(mat data, int width, int height);


cv::Mat_<double> to_cvmat(const arma::Mat<double> &src) {
	return cv::Mat_<double>{int(src.n_cols), int(src.n_rows), const_cast<double*>(src.memptr())};
}

arma::mat to_arma(const cv::Mat_<double> &src) {
	arma::mat dst(reinterpret_cast<double*>(src.data), src.cols, src.rows);
	return dst;
}

mat bilinear_resize(mat data, int width, int height) {

	//convert from arma mat to cv mat
	cv::Mat src = to_cvmat(data);
	cv::Mat temp(height, width, CV_32F); //this needs to change

	//bilinear resize using openCV
	cv::resize(src, src, temp.size());

	//convert back into arma mat
	//mat dst = to_arma(src);

	return to_arma(src);
}

