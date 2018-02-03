/*
 * help.cpp
 *
 *  Created on: Jan 30, 2018
 *      Author: s27280 OCdt James Dolman
 */

#include <iostream>
#include <fstream>
#include <armadillo>
#include "energy_calc.hpp"

using namespace std;
using namespace arma;
void init_variables();

	//variables for energy_calc
    mat e_smooth;

	mat u;
	mat v;
	//mat x;

int main() {

	init_variables();
    //  tested and works. psi_function
    //  x = psi_function(x);

	//tested and works. generate_esmooth
  	e_smooth = generate_esmooth(u,v);
	cout<<"Out of in test main()"<<endl;

    //e_smooth.save("mats/test_energy_calc/Outputs/e_smoothv2-c", raw_ascii);
	cout<<"E_smooth width: "<<e_smooth.n_cols<<endl;
	cout<<"E_smooth hight: "<<e_smooth.n_rows<<endl;

	return 0;
}



void init_variables(){
	//variables loaded in for energy_calc
	//x.load("mats/test_energy_calc/Inputs/x_psi.txt");
	u.load("mats/test_energy_calc/Inputs/u.txt");
	v.load("mats/test_energy_calc/Inputs/v.txt");



}

