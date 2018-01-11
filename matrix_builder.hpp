/*
energy_calc.hpp
Jacob Wiebe & James Dolman
Rev1: Nov 2017
*/

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

// mat u ;
// mat e_smooth;
// int height = 0, width = 0;

//M: constructMatrix
void build_matrix (mat& A, vec& b, mat img2_dx, mat img2_dy, mat img_z,
		mat dxx, mat dxy, mat dyy, mat img_dxz, mat e_data,
		mat e_smooth, mat u, mat v, double gamma){

	int height = u.n_rows;
	int width = u.n_cols;

	int e_height = e_smooth.n_rows;
	int e_width = e_smooth.n_cols;

  //top and bottom row to zero
  e_smooth.head_rows(1) = zeros<rowvec>(e_width);
  e_smooth.tail_rows(1) = zeros<rowvec>(e_width);

  //left and right col to zero
  e_smooth.head_cols(1) = zeros<vec>(e_height);
  e_smooth.tail_cols(1) = zeros<vec>(e_height);

  cout << "e_smooth after zero frame: " << endl << e_smooth << endl;

	// M: tmp = repmat( 1 : 2 * ht * wt, 6, 1 ) ;
	// M: ros = tmp(:);
	int temp4repmat = 2*height*width*6;
	int tempPop = 1;
	vec rows = (temp4repmat); //column vector of size temp4repmat
	for (int i = 0; i<temp4repmat ; i++){
		rows(i) = tempPop;
		if(i%6 == 0){
			tempPop++;
		}
	}

	//M:cols = rows

  return;
}
