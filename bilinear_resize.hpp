#pragma once

#include <iostream>

#include <armadillo>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace std;
using namespace arma;


cv::Mat_<double> to_cvmat(const arma::Mat<double> &src);
arma::mat to_arma(const cv::Mat_<double> &src);
mat bilinear_resize(mat data,int scalefactor);
mat bilinear_resize(mat data, mat frame);


cv::Mat_<double> to_cvmat(const arma::Mat<double> &src) {

	int ht = src.n_rows;
	int wt = src.n_cols;
	/*
	arma::mat src2(src.n_cols, src.n_rows );

	cout << "src2 in the to_cvmat routine pre rotation" << endl << src2 << endl;

	uword i = 0, j = 0;
	//This loop rotates the image 90 degrees
	for (i; i<src.n_rows; i++) {
		for (j; j < src.n_cols; j++) {
			src2.at(j, i) = src.at(i, j);
		}
		j = 0;
	}
	cout << "src2 in the to_cvmat routine post rotation" << endl << src2 << endl;*/
	
	

	cv::Mat dst = cv::Mat_<double>{int(src.n_rows), int(src.n_cols), const_cast<double*>(src.memptr())};
	cv::Mat dst_t(dst.rows, dst.cols, CV_64F);
	//cout << "dst pre " << endl << dst << endl;
	//int inc = max(dst.cols, dst.rows);

	

	int i = 0, j = 0;
	//This loop rotates the image 90 degrees
	int n = 0;
	for (i; i< dst.cols; i++) {
		for (j; j < dst.rows; j++) {
			dst_t.at<double>(j, i) = dst.at<double>(n);
			n++;
		}
		j = 0;
	}
	//cout << "dst in the to_cvmat routine post rotation" << endl << dst_t << endl;

	return dst_t;
}

arma::mat to_arma(const cv::Mat_<double> &src) {
	arma::mat dst(reinterpret_cast<double*>(src.data), src.cols, src.rows);
	mat dstF(dst.n_cols, dst.n_rows); 
	uword i = 0, j = 0;
	//This loop rotates the image 90 degrees
	for (i ; i<dst.n_rows ; i++) {
		for (j; j < dst.n_cols; j++) {
				dstF(j, i) = dst(i, j);
		}
		j = 0;
	}
	return dstF;
}

mat bilinear_resize(mat data, double scalefactor) {

	cv::Mat src = to_cvmat(data);
	cv::Mat dst;

	cv::resize(src, dst, cv::Size(), scalefactor, scalefactor, CV_INTER_AREA);

	//cout << "done resizing" << endl;
	src.release();

	return to_arma(dst);

	
	//convert from arma mat to cv mat
	//double ht = data.n_rows;
	//double wt = data.n_cols;

	//cout << "height:" << ht << endl;
	//cout << "width:" << wt << endl;
	//cout << scalefactor << endl;
	//cout << ht*scalefactor << endl;
	//cout << wt*scalefactor << endl;
	//cin.get();

	//for (int i = 0; i < 10; i++) {
	//	cout << data.at(i, i) << endl;
	//}

	
	//for (int i = 0; i < 10; i++) {
	//	cout << src.at<double>(i, i) << endl;
	//}
	
	
	//cv::Mat frame(ht*scalefactor, wt*scalefactor, CV_64F);
	//cv::Mat dst(ht*scalefactor, wt*scalefactor, CV_64F);
	

	//uchar* pt = cv::Mat::ptr(&src);

	//cv::resize(src, dst, frame.size(), 0, 0, CV_INTER_LINEAR);
	
	//cv::resize(src, src, Size(), scalefactor, scalefactor, INTER_LINEAR);
	//for (int i = 0; i < 10; i++) {
	//	cout << dst.at<double>(i, i) << endl;
	//}
	
	//frame.release();

	
	//cin.get();
	
}

mat bilinear_resize(mat data, mat frame) {

	/*cout << "frame is " << frame.n_rows << "x" << frame.n_cols << endl;
	cout << frame << endl;*/

	cv::Mat temp = to_cvmat(frame);

	int ht = temp.rows;
	int wt = temp.cols;

	/*cout << "temp is " << ht << "x" << wt << endl;
	cout << temp << endl;*/

	cv::Mat src = to_cvmat(data);
	cv::Mat dst;

	cv::resize(src, dst, temp.size(), 0, 0, CV_INTER_AREA);

	/*cout << "done resizing 2" << endl;
	src.release();*/

	return to_arma(dst);
	
	/*uword ht = frame.n_rows;
	uword wt = frame.n_cols;*/

	//cout << "my mat" << endl;
	//cout << data.at(3) << endl; //<< data.submat(0, 0, ht-1, wt-1) << endl;

	//convert from arma mat to cv mat
	//cv::Mat src = to_cvmat(data);         //This is the problem
	//
	//cv::Mat temp = to_cvmat(frame);
	//
	//cv::resize(src, src, temp.size(), 0, 0, CV_INTER_AREA);

	//return to_arma(src);
}

