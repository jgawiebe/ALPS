/*
energy_calc.hpp
Jacob Wiebe & James Dolman
Rev1: Nov 2017
Rev2: Feb 2018 TESTING PASSED
*/

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;
//e_smooth.save("mats/test_matrix_builder/Outputs/e_smooth-c", raw_ascii);
// mat u ;
// mat e_smooth;
// int height = 0, width = 0;

//M: constructMatrix
tuple<sp_mat,vec> build_matrix (mat img2_dx, mat img2_dy, mat img_z,
		mat dxx, mat dxy, mat dyy, mat dxz, mat dyz, mat e_data,
		mat e_smooth, mat u, mat v, double gamma){
	cout<<"In build_matrix"<<endl;
	sp_mat A;
	vec b;
	uword height = u.n_rows;
	uword width = u.n_cols;

	uword e_height = e_smooth.n_rows;
	uword e_width = e_smooth.n_cols;
//-----------------------------------------------------------------------------------------------
  //top and bottom row to zero
  e_smooth.head_rows(1) = zeros<rowvec>(e_width);
  e_smooth.tail_rows(1) = zeros<rowvec>(e_width);

  //left and right col to zero
  e_smooth.head_cols(1) = zeros<vec>(e_height);
  e_smooth.tail_cols(1) = zeros<vec>(e_height);
//-------------------------------------------------------------------------------------------------


	// M: tmp = repmat( 1 : 2 * ht * wt, 6, 1 ) ;
	// M: ros = tmp(:);
	uword temp4repmat = 2*height*width*6;
	int tempPop = 1;


	urowvec rows(temp4repmat); //column vector of size temp4repmat

	for (uword i = 0; i<temp4repmat ; i++){
		if(i%6 == 0 && i != 0){
			tempPop++;
		}
		rows(i) = tempPop;
	}
   /////////////////Where I started 11 Jan 2018////////////////////

	//M:cols = rows
	urowvec cols = rows;

	//M:vals = zeros( size( rows ) )
	vec vals(temp4repmat, fill::zeros);

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

	//THIS WILL NEED MORE WORK AS E_SUM NEEDS TO BE INDEXED DIFFERENTLY THAN E_SMOOTH>DONE

	mat e_sum(height , width , fill::zeros);

	//NOTE MAY HAVE TO REVISE THE IF STATEMENTS TO MAKE THEM HAVE 3 IF STATEMENTS> DONE
	uword ik = 0;
	uword jk = 0;
	for (uword i = 2; i<e_height ; i=i+2){
		for (uword j = 2; j<e_width ; j=j+2){
				if (i < (2*height)){//ht good
					if(j <(2*width)){//ht good //width good

						e_sum(ik,jk) += e_smooth(i,j-1); //(2)
						e_sum(ik,jk) += e_smooth(i-1,j); //(4)
						e_sum(ik,jk) += e_smooth(i-2,j-1); //(1)
						e_sum(ik,jk) += e_smooth(i-1,j-2); //(3)
					}else{ //ht good
						e_sum(ik,jk) += e_smooth(i,j-1);
						e_sum(ik,jk) += e_smooth(i-1,j);
						e_sum(ik,jk) += e_smooth(i-2,j-1); //(2)+(4)+(1), just height good
					}
				} else if(j <(2*width)){//width good
					e_sum(ik,jk) += e_smooth(i,j-1);
					e_sum(ik,jk) += e_smooth(i-1,j);
					e_sum(ik,jk) += e_smooth(i-1,j-2); //(2)+(4)+(3), just width good
				}else{
					e_sum(ik,jk) += e_smooth(i,j-1);
					e_sum(ik,jk) += e_smooth(i-1,j); //(2)+(4), neither height or width good
				}
				jk++;
			}//for loop for the columns
		jk = 0;
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

	//4 temp matrixes used in initializing vals,
	mat tmp1(height , width , fill::zeros); //NOTE NO IT WONT IT WILL BE HALF THE SIZE AS YOU TAKE EVERY 2ND
	mat tmp2 = tmp1;
	mat tmp3 = tmp1;
	mat tmp4 = tmp1;
	ik = 0;
    jk = 0;
	//NOTE MAY HAVE TO REVISE THE IF STATEMENTS TO MAKE THEM HAVE 3 IF STATEMENTS>Done
	//INDEXING CORRECT HERE AS TMP JUST TAKES EVERY SECOND ONE OF E_SMOOTH
	for (uword i = 2; i<e_height ; i=i+2){
			for (uword j = 2; j<e_width ; j=j+2){
				 tmp1(ik,jk) = e_smooth(i-1,j-2); //M:tmp = aE_smooth( 2 : 2 : end, 1 : 2 : 2 * wt ) ;
				 tmp2(ik,jk) = e_smooth(i-1,j);   //M:tmp = aE_smooth( 2 : 2 : end, 3 : 2 : end ) ;
				 tmp3(ik,jk) = e_smooth(i-2,j-1);//M:tmp = aE_smooth( 1 : 2 : 2 * ht, 2 : 2 : end ) ;
				 tmp4(ik,jk) = e_smooth(i,j-1);       //M:tmp = aE_smooth( 3 : 2 : end, 2 : 2 : end ) ;

					jk++;
				}//for loop for the columns
			jk = 0;
			ik++;
		}//for loop for the rows


	//M:vals( 3 : 12 : end ) = uapp(:) ;
	//M:vals( 10 : 12 : end ) = vapp(:) ;
	//M:vals( 4 : 12 : end ) = uvapp(:) ;
	//M:vals( 9 : 12 : end ) = vuapp(:) ;

	//M:vals( 1 : 12 : end ) = -tmp(:) ;
	//M:vals( 7 : 12 : end ) = -tmp(:) ;
	//M:vals( 6 : 12 : end ) = -tmp(:) ;
	//M:vals( 12 : 12 : end ) = -tmp(:) ;

	//M:vals( 2 : 12 : end ) = -tmp(:) ;
	//M:vals( 8 : 12 : end ) = -tmp(:) ;
	//M:vals( 5 : 12 : end ) = -tmp(:) ;
	//M:vals( 11 : 12 : end ) = -tmp(:) ;

	vectorise(uapp); //turns matrix into column vector
	vectorise(vapp);
	vectorise(uvapp);
	ik = 0;

	for (uword i = 11; i<temp4repmat ; i=i+12){
		vals(i-9) = uapp(ik);
		vals(i-2) = vapp(ik);
		vals(i-8) = uvapp(ik);
		vals(i-3) = uvapp(ik);

		vals(i-11)= -(tmp1(ik));
		vals(i-5) = -(tmp1(ik));
		vals(i-6) = -(tmp2(ik));
		vals(i)   = -(tmp2(ik));

		vals(i-10)= -(tmp3(ik));
		vals(i-4) = -(tmp3(ik));
		vals(i-7) = -(tmp4(ik));
		vals(i-1) = -(tmp4(ik));

		ik++;
	}
	//vals.save("mats/test_matrix_builder/Outputs/valsm-c", raw_ascii);
	//STOPPING HERE 15 Jan, LOOK AT THE NOTES AND BUILD THE LAST TWO TEMP MATRIXES
	//Jan 17 Start: Fixed all notes I had
	//M:upad = padarray( u, [1 1] ) ;
	//M:vpad = padarray( v, [1 1] ) ;
	mat upad(height + 2, width + 2, fill::zeros);
	mat vpad(v.n_rows + 2, v.n_cols + 2, fill::zeros);
//--------------------------------------------------------------------------------------
	upad.submat(1,1,height,width)=u;
	vpad.submat(1,1,v.n_rows,v.n_cols)=v;
//---------------------------------------------------------------------------------------
	/*% Computing the constant terms for the first of the Euler Lagrange equations
	pdfaltsumu = aE_smooth(2:2:end, 1:2:2*wt) .* ( upad(2:ht+1, 1:wt) -  upad(2:ht+1, 2:wt+1) ) + (1)
				aE_smooth( 2:2:end, 3:2:end) .*  ( upad(2:ht+1, 3:end) - upad(2:ht+1, 2:wt+1) ) + (2)
				aE_smooth( 1:2:2*ht, 2:2:end) .* ( upad(1:ht, 2:wt+1) -  upad(2:ht+1, 2:wt+1) ) + (3)
				aE_smooth( 3:2:end, 2:2:end) .*  ( upad(3:end, 2:wt+1) - upad(2:ht+1, 2:wt+1) ) ; (4) */

	mat pdfaltsumu( height , width , fill::zeros);
	mat pdfaltsumv  = pdfaltsumu;

	uvec rowstk1(e_height/2,fill::zeros) ;
	uvec colstk1(e_width/2,fill::zeros) ;
	uvec rowstk2(e_height/2,fill::zeros) ;
	uvec colstk2(e_width/2,fill::zeros) ;
	uvec rowstk3(e_height/2,fill::zeros) ;
	uvec colstk3(e_width/2,fill::zeros) ;
	uvec rowstk4(e_height/2,fill::zeros) ;
	uvec colstk4(e_width/2,fill::zeros) ;

	tmp1.resize(upad.n_rows,upad.n_cols);//tmp size upad
	tmp2 =tmp1;
	tmp3 =tmp1;
	tmp4 =tmp1;
mat tmp5 =tmp1;

	tmp1 = upad.submat(1,1,height,width) ; //(0)
	tmp2 = upad.submat(1,0,height,width-1);// (1)
	tmp3 = upad.submat(1,2,height,(upad.n_cols)-1);// (2)
	tmp4 = upad.submat(0,1,height-1,width);// (3)
	tmp5 = upad.submat(2,1,(upad.n_rows)-1,width);//(4)

	//----------------------------CURRENTLY WORKING IN THIS LOOP, TRYING TO BUILD THE INDEX VECTORS FOR THE TWO LOOPS BELOW

	ik = 0;
	for (uword i = 2; i<e_height ; i=i+2){
		rowstk1(ik) = i-1;
		rowstk2(ik) = i-1;
		rowstk3(ik) = i-2;
		rowstk4(ik) = i;
		ik++;
	}//for loop for the rows

	jk = 0;
	for (uword j = 2; j<e_width ; j=j+2){
		colstk1(jk) = j-2;
		colstk2(jk) = j;
		colstk3(jk) = j-1;
		colstk4(jk) = j-1;
		jk++;
	}//for loop for the columns

	pdfaltsumu += e_smooth.submat(rowstk1,colstk1)%(tmp2 - tmp1);
	pdfaltsumu += e_smooth.submat(rowstk2,colstk2)%(tmp3 - tmp1);
	pdfaltsumu += e_smooth.submat(rowstk3,colstk3)%(tmp4 - tmp1);
	pdfaltsumu += e_smooth.submat(rowstk4,colstk4)%(tmp5 - tmp1);

/*	% Computing the constant terms for the second of the Euler Lagrange equations
	pdfaltsumv = aE_smooth(2:2:end, 1:2:2*wt) .* ( vpad(2:ht+1, 1:wt) - vpad(2:ht+1, 2:wt+1) ) + ...
				aE_smooth( 2:2:end, 3:2:end) .* ( vpad(2:ht+1, 3:end) - vpad(2:ht+1, 2:wt+1) ) + ...
				aE_smooth( 1:2:2*ht, 2:2:end) .* ( vpad(1:ht, 2:wt+1) - vpad(2:ht+1, 2:wt+1) ) + ...
				aE_smooth( 3:2:end, 2:2:end) .* ( vpad(3:end, 2:wt+1) - vpad(2:ht+1, 2:wt+1) ) ;
*/
	tmp1 = vpad.submat(1,1,height,width) ; //(0)
	tmp2 = vpad.submat(1,0,height,width-1);// (1)
	tmp3 = vpad.submat(1,2,height,(upad.n_cols)-1);// (2)
	tmp4 = vpad.submat(0,1,height-1,width);// (3)
	tmp5 = vpad.submat(2,1,(upad.n_rows)-1,width);//(4)
    //use same row and col vectors as they are the same as pdfaltsumu
	pdfaltsumv += e_smooth.submat(rowstk1,colstk1)%(tmp2 - tmp1);
	pdfaltsumv += e_smooth.submat(rowstk2,colstk2)%(tmp3 - tmp1);
	pdfaltsumv += e_smooth.submat(rowstk3,colstk3)%(tmp4 - tmp1);
	pdfaltsumv += e_smooth.submat(rowstk4,colstk4)%(tmp5 - tmp1);

	//M:constu = E_Data .* ( Ikx .* Ikz + gamma * ( Ixx .* Ixz + Ixy .* Iyz ) ) - pdfaltsumu ;
	//M:constv = E_Data .* ( Iky .* Ikz + gamma * ( Ixy .* Ixz + Iyy .* Iyz ) ) - pdfaltsumv ;
	mat constu = e_data % ((img2_dx %  img_z) + gamma * ((dxx % dxz) + (dxy % dyz))) - pdfaltsumu;
	mat constv = e_data % ((img2_dy %  img_z) + gamma * ((dxy % dxz) + (dyy % dyz))) - pdfaltsumv;

	//M:b = zeros( 2 * ht * wt, 1 ) ;
	//M:b(1:2:end) = -constu(:) ;
	//M:b(2:2:end) = -constv(:) ;

	b.zeros((2*height*width));
	vectorise(constu);
	vectorise(constv);

	ik = 0;

	for (uword i = 1; i<b.n_rows ; i=i+2){
		b(i) = -(constv(ik));
		b(i-1) = -(constu(ik));
		ik++;
	}

	//M:ind = find(cols > 0) ;
	//M:rows = rows( ind ) ;
	//M:cols = cols( ind ) ;
	//M:vals = vals( ind ) ;


	for (uword i = 0; i<cols.n_cols ; i=i+1){
			if(cols(i) > 0){

			}else{
				cols.shed_col(i);
				rows.shed_col(i);
				vals.shed_row(i);
			}
	}
	for (uword i = 0; i<cols.n_cols ; i=i+1){
			if(cols(i) > 0){

			}else{
				cols.shed_col(i);
				rows.shed_col(i);
				vals.shed_row(i);
			}
	}

	//M:ind = find(cols < ( 2 * ht * wt + 1 ) ) ;
	//M:rows = rows( ind ) ;
	//M:cols = cols( ind ) ;
	//M:vals = vals( ind ) ;
	for (uword i = 0; i<cols.n_cols ; i=i+1){
			if(cols(i) < (2*height*width+1)){

			}else{
				cols.shed_col(i);
				rows.shed_col(i);
				vals.shed_row(i);
			}
	}
	for (uword i = 0; i<cols.n_cols ; i=i+1){
			if(cols(i) < (2*height*width+1)){

			}else{
				cols.shed_col(i);
				rows.shed_col(i);
				vals.shed_row(i);
			}
	}

	for (uword i = 0; i<cols.n_cols ; i=i+1){

		cols(i) = cols(i) -1;
		rows(i) = rows(i) -1;


	}

	//cols.save("mats/test_matrix_builder/Outputs/cols-c", raw_ascii);
	//M:A = sparse (rows,cols,vals) ;

	//---------------Create Square Matrix-------------------------------------------
	cout<<" SIZE OF COLS"<<cols.n_cols<<endl;
	cout<<" MAX OF COL S"<<max(cols)<<endl;
	cout<<" SIZE OF B   "<<b.n_rows<<endl;
	//do that here
    umat locations = join_cols(rows,cols);

    sp_mat C1(locations,vals);

    A = C1;
    cout<<" SIZE OF A  "<<A.n_cols<<endl;
	//end Jan 17 - just need to do pdfaltsumv's > Done morning of 19th.
	return make_tuple(A,b);
}
