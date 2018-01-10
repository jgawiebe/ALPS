/*
 energy_calc.hpp
 Jacob Wiebe & James Dolman
 Rev1: Nov 2017
 */

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

//compile with: g++ <funtion>.cpp -o test -O2 -larmadillo; ./test

//M: optic_flow_brox
int main() {

	//non-matrices: alpha, dt, gamma, ht, i, num_levels, wt
	double alpha = 30.0, gamma = 80.0;
	int num_levels = 40, outer_iter = 3, inner_iter = 500;

	//get size of image
	//int ht = 0;
	
	

	//define u and v

	//check this loop
	for (int i = 0; i < num_levels; i++) {

		//perform calc derivatives from optic_flow.hpp
		//perform optic flow to get du and dv

		u = u + du;
		v = v + dv;
	}

	return 1;
}
