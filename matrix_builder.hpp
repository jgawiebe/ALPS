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
		mat dxx, mat dxy, mat dyy, mat dxz, mat dyz, mat e_data,
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

	for (uword i = 0; i<temp4repmat ; i++){
		if(i%6 == 0){
			tempPop++;
		}
		rows(i) = tempPop;
	}
   /////////////////Where I started 11 Jan 2018////////////////////
	//M:cols = rows
	vec cols = rows;

	//M:vals = zeros( size( rows ) )
	vec vals(size(rows), fill::zeros);

	//MatLab is 1 indexed and C++ is 0 indexed. So thats why i in the
	//loop is 1 less than in the Matlab comment.
	//M:cols(1:6:end) = rows(1:6:end) - 2 * ht ;	% x-1
	for (uword i = 0; i<temp4repmat ; i=i+6){
			cols(i) = rows(i) - (2*height);
	}
	//M:cols(2:6:end) = rows(2:6:end) - 2 ;			% y-1
	for (uword i = 1; i<temp4repmat ; i=i+6){
			cols(i) = rows(i) - 2;
	}

	//M:cols(9:12:end) = rows(9:12:end) - 1 ;		% v
	for (uword i = 8; i<temp4repmat ; i=i+12){
			cols(i) = rows(i) - 1;
	}

	//M:cols(4:12:end) = rows(4:12:end) + 1 ;		% u
	for (uword i = 3; i<temp4repmat ; i=i+12){
			cols(i) = rows(i) + 1;
	}

	//M:cols(5:6:end) = rows(5:6:end) + 2 ;			% y+1
	for (uword i = 4; i<temp4repmat ; i=i+6){
			cols(i) = rows(i) + 2;
	}

	//M:cols(6:6:end) = rows(6:6:end) + 2 * ht ;	% x+1
	for (uword i = 5; i<temp4repmat ; i=i+6){
			cols(i) = rows(i) + (2*height);
	}

        //start for 14 Jan 17

	//M:E_sum = (1) aE_smooth( 1 : 2 : 2 * ht, 2 : 2 : end ) + (2)aE_smooth( 3 : 2 : end, 2 : 2 : end ) +...
	// (3)aE_smooth( 2 : 2 : end, 1 : 2 : 2 * wt ) + (4)aE_smooth( 2 : 2 : end, 3 : 2 : end ) ;
	//loops built around aE_smooth( 3 : 2 : end, 2 : 2 : end ) so the loops indexing matches that
	//the initial values are 1 less at initialization due to 0 indexing vs 1 indexing

	//through workspace break point analysis it was determined that height and width (components of matrix u)
	//are never larger than half of the dimensions of e smooth. Therefore, we can use the dimensions of
	//smooth for e_sum.

	//THIS WILL NEED MORE WORK AS E_SUM NEEDS TO BE INDEXED DIFFERENTLY THAN E_SMOOTH
	//I think the ik jk stuff fixes our problem...
	mat e_sum(size(e_smooth), fill::zeros);
	uword ik = 0;
	uword jk = 0;
	//NOTE MAY HAVE TO REVISE THE IF STATEMENTS TO MAKE THEM HAVE 3 IF STATEMENTS
		//ONE FOR HEIGHT, ONE FOR WIDTH THEN ONE FOR NEITHER
	for (uword i = 2; i<e_height ; i=i+2){
		for (uword j = 1; j<e_width ; j=j+2){
				if ((i < 2*height) && (j <2*width)){
					e_sum(ik,jk) += (e_smooth(i,j) + e_smooth(i-1,j+1) + e_smooth(i-2,j) + e_smooth(i-1,j-1)) ; //(2)+(4)+(1)+(3)
				} else{
					e_sum(ik,jk) += (e_smooth(i,j) + e_smooth(i-1,j+1)) ; //(2)+(4)
				}
				jk++;
			}//for loop for the columns
		ik++;
	}//for loop for the rows

	//M:uapp = E_Data .* ( Ikx .^ 2 + gamma * ( Ixx .^ 2 + Ixy .^ 2 ) ) + E_sum ;
	mat uapp = e_data % (square(img2_dx) + gamma * (square(dxx)+square(dxy))) + e_sum;

	//M:vapp = E_Data .* ( Iky .^ 2 + gamma * ( Iyy .^ 2 + Ixy .^ 2 ) ) + E_sum ;
	mat vapp = e_data % (square(img2_dy) + gamma * (square(dyy)+square(dxy))) + e_sum;

	//M:uvapp = E_Data .* ( Ikx .* Iky + gamma * ( Ixx .* Ixy + Iyy .* Ixy ) ) ;
	mat uvapp = e_data % ((img2_dx % img2_dy) + gamma * ((dxx % dxy) + (dyy % dxy)));

	//M:vuapp = E_Data .* ( Ikx .* Iky + gamma * ( Ixx .* Ixy + Iyy .* Ixy ) ) ;
	// vuapp declaration is the same as the uvapp, therefore vuapp is a duplicate
	// and will not be used in our C++ implementation of the code.

	//end for 14 Jan 17

	//start 15 Jan 17
	//M:vals( 3 : 12 : end ) = uapp(:) ;
	//M:vals( 10 : 12 : end ) = vapp(:) ;
	//M:vals( 4 : 12 : end ) = uvapp(:) ;
	//M:vals( 9 : 12 : end ) = vuapp(:) ;

	//4 temp matrixes used in initializing vals,
	mat tmp1(size(e_smooth), fill::zeros);
	mat tmp2 = tmp1;
	mat tmp3 = tmp1;
	mat tmp4 = tmp1;
    ik = 0;
	jk = 0;
	//NOTE MAY HAVE TO REVISE THE IF STATEMENTS TO MAKE THEM HAVE 3 IF STATEMENTS
	//ONE FOR HEIGHT, ONE FOR WIDTH THEN ONE FOR NEITHER
	
	for (uword i = 2; i<e_height ; i=i+2){
			for (uword j = 2; j<e_width ; j=j+2){
					if ((i < 2*height) && (j <2*width)){
						tmp1(ik,jk) = e_smooth(i-1,j-2); //M:tmp = aE_smooth( 2 : 2 : end, 1 : 2 : 2 * wt ) ;
						tmp2(ik,jk) = e_smooth(i-1,j); //M:tmp = aE_smooth( 2 : 2 : end, 3 : 2 : end ) ;
					} else {
						tmp2(ik,jk) = e_smooth(i-1,j);
					}
					jk++;
				}//for loop for the columns
			ik++;
		}//for loop for the rows

	vectorise(uapp); //turns matrix into column vector
	vectorise(vapp);
	vectorise(uvapp);
	ik = 0;
	
	for (uword i = 9; i<temp4repmat ; i=i+12){
		vals(i-7) = uapp(ik);
		vals(i)   = vapp(ik);
		vals(i-6) = uvapp(ik);
		vals(i-1) = uvapp(ik);
		ik++;
	}

	//READ ME: there is three notes that need to be done. they are in CAPS
	//both have to do with indexing the loops to make sure we are getting the right stuff

	//NOTE THIS WILL NEED TO BE APPLIED TO ALL LOOPS WHEN YOU WORK ON THEM
	/*The indices of elements are specified via the uword type, which is a
	 * typedef for an unsigned integer type. When using loops to access
	 * elements, it's best to use uword instead of int. For example:
	 * for(uword i=0; i<X.n_elem; ++i) { X(i) = ... }*/

	//STOPPING HERE 15 Jan, LOOK AT THE NOTES AND BUILD THE LAST TWO TEMP MATRIXES
	
	//Just did a quick touch up on the uword loops. Then will also go through the notes on
	//the if statements when dealing with width and height. 
	//Finally then calc temp 3 and temp 4. Then finish off creating the vals in the
	//about for loop. 
  return;
}


