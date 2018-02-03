/*
 * help.cpp
 *
 *  Created on: Jan 30, 2018
 *      Author: s27280 OCdt James Dolman
 */

#include <iostream>
#include <fstream>
#include <armadillo>
#include "gaussian_smooth.hpp"

using namespace std;
using namespace arma;
void init_variables();

	//variables for energy_calc
    mat img;
    mat gaussian;
    double scale;

int main() {

	init_variables();
	cout<<"IMG width: "<<img.n_cols<<endl;
	cout<<"IMG hight: "<<img.n_rows<<endl;
	//pass. Fix through addition of submat on line 85.
  	gaussian = g_smooth(img,scale);
	cout<<"Out of in test main()"<<endl;

    gaussian.save("mats/test_gaussian_smooth/Outputs/gaussian-c", raw_ascii);
	cout<<"Gaussian width: "<<gaussian.n_cols<<endl;
	cout<<"Gaussian hight: "<<gaussian.n_rows<<endl;

	return 0;
}



void init_variables(){
	//variables loaded in for energy_calc
	img.load("mats/test_gaussian_smooth/Inputs/imgv2grey.txt");
	scale = 0.1285;



}

