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


	vec rows(temp4repmat); //column vector of size temp4repmat

	for (int i = 0; i<temp4repmat ; i++){
		rows(i) = tempPop;
		if(i%6 == 0){
			tempPop++;
		}
	}
   /////////////////Where I started 11 Jan 2018////////////////////
	//M:cols = rows
	vec cols = rows;

	//M:vals = zeros( size( rows ) )
	mat vals(size(rows), fill::zeros);

	//MatLab is 1 indexed and C++ is 0 indexed. So thats why i in the
	//loop is 1 less than in the Matlab comment.
	//M:cols(1:6:end) = rows(1:6:end) - 2 * ht ;	% x-1
	for (int i = 0; i<temp4repmat ; i=i+6){
			cols(i) = rows(i) - (2*height);
	}
	//M:cols(2:6:end) = rows(2:6:end) - 2 ;			% y-1
	for (int i = 1; i<temp4repmat ; i=i+6){
			cols(i) = rows(i) - 2;
	}

	//M:cols(9:12:end) = rows(9:12:end) - 1 ;		% v
	for (int i = 8; i<temp4repmat ; i=i+12){
			cols(i) = rows(i) - 1;
	}

	//M:cols(4:12:end) = rows(4:12:end) + 1 ;		% u
	for (int i = 3; i<temp4repmat ; i=i+12){
			cols(i) = rows(i) + 1;
	}

	//M:cols(5:6:end) = rows(5:6:end) + 2 ;			% y+1
	for (int i = 4; i<temp4repmat ; i=i+6){
			cols(i) = rows(i) + 2;
	}

	//M:cols(6:6:end) = rows(6:6:end) + 2 * ht ;	% x+1
	for (int i = 5; i<temp4repmat ; i=i+6){
			cols(i) = rows(i) + (2*height);
	}

	//M:E_sum = aE_smooth( 1 : 2 : 2 * ht, 2 : 2 : end ) + aE_smooth( 3 : 2 : end, 2 : 2 : end ) +
	//aE_smooth( 2 : 2 : end, 1 : 2 : 2 * wt ) + aE_smooth( 2 : 2 : end, 3 : 2 : end ) ;
	mat e_sum(size(e_smooth));

	//wow this is tough... im done for today


  return;
}

