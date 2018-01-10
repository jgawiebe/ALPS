/*
energy_calc.hpp
Jacob Wiebe & James Dolman
Rev1: Nov 2017
*/

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

mat u ;
mat E_smooth;
int height = 0;
int width = 0;

//M: constructMatrix
void build_matrix (mat& A, vec& b, img2_dx,mat img2_dy,mat img_z,
		mat dxx,mat dxy,mat dyy,mat img_dxz,mat e_data,
		mat E_smooth,mat u,mat v,double gamma){
	//Question here, what if matrix u is not a square?
	//how do we use size() to get height and width then.

	height = size(u);
	width = size(u);
	//or do height = u.n_rows; //Question: does that work?

	//way to test what size u produces
	//mat u (5,9);
	//cout << "size of u: " << size(u) << endl;
	//if this size doesn't work, try doing this and test again
	//height = u.n_rows;
	//width = u.n_cols;

	//next frame off the matrix E_smooth with zeroes (not sure if correct)
	//please test
	int E_smoothheight = E_smooth.n_rows;
	int E_smoothwidth = E_smooth.n_cols;
	//top and bottom row to zero
	E_smooth.row(0).zeros;
	E_smooth.row(E_smoothheight-1).zeros;
	//left and right cols to zero
	E_smooth.col(0).zeros;
	E_smooth.col(E_smoothwidth-1).zeros;
	//see if this frames it off with zeros
	//main concern is if the .zeros() function can be used this way.
	//cout << "E_smooth after zero frame: " << E_smooth << endl;



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
