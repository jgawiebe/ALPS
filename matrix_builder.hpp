//matrix_builder.hpp holds two functions, build_sparse() and 
//build_matrix(). These functions are heavily mathematical and 
//produce the same result. This result is a matrix A, called 
//the coefficient matrix, and a b vector that are returned for
//the linear solver to solve the system of linear equations
//A(x) = b where x is the solution. 

//The difference between the two functions is how they package 
//the data they output. build()_sparse stores the outputted 
//coeffient matrix in an sp_mat (sparse matrix) format and 
//build_matrix() stores the coeffient matrix as two uvec
//(unsigned vectors) and a mat. 


#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

//M: constructMatrix
tuple<sp_mat, vec> build_sparse(mat img2_dx, mat img2_dy, mat img_z,
	mat dxx, mat dxy, mat dyy, mat dxz, mat dyz, mat e_data,
	mat e_smooth, mat u, mat v, double gamma) {
	
	sp_mat A;
	vec b;
	uword height = u.n_rows;
	uword width = u.n_cols;
	uword e_height = e_smooth.n_rows;
	uword e_width = e_smooth.n_cols;

	e_smooth.head_rows(1) = zeros<rowvec>(e_width);
	e_smooth.tail_rows(1) = zeros<rowvec>(e_width);

	e_smooth.head_cols(1) = zeros<vec>(e_height);
	e_smooth.tail_cols(1) = zeros<vec>(e_height);

	uword temp4repmat = 2 * height*width * 6;
	uword tempPop = 1;
	urowvec rows(temp4repmat); //column vector of size temp4repmat

	for (uword i = 0; i<temp4repmat; i++) {
		if (i % 6 == 0 && i != 0) {
			tempPop++;
		}
		rows(i) = tempPop;
	}

	urowvec cols = rows;
	vec vals(temp4repmat, fill::zeros);

	for (uword i = 0; i<temp4repmat; i = i + 6) {
		cols(i) = rows(i) - (2 * height);
	}
	
	for (uword i = 1; i<temp4repmat; i = i + 6) {
		cols(i) = rows(i) - 2;
	}

	for (uword i = 8; i<temp4repmat; i = i + 12) {
		cols(i) = rows(i) - 1;
	}

	for (uword i = 3; i<temp4repmat; i = i + 12) {
		cols(i) = rows(i) + 1;
	}

	for (uword i = 4; i<temp4repmat; i = i + 6) {
		cols(i) = rows(i) + 2;
	}

	for (uword i = 5; i<temp4repmat; i = i + 6) {
		cols(i) = rows(i) + (2 * height);
	}

	mat e_sum(height, width, fill::zeros);

	uword ik = 0;
	uword jk = 0;

	for (uword i = 2; i<e_height; i = i + 2) {
		for (uword j = 2; j<e_width; j = j + 2) {
			e_sum(ik, jk) += e_smooth(i, j - 1); //(2)
			e_sum(ik, jk) += e_smooth(i - 1, j); //(4)
			e_sum(ik, jk) += e_smooth(i - 2, j - 1); //(1)
			e_sum(ik, jk) += e_smooth(i - 1, j - 2); //(3) 
											
			jk++;
		}//for loop for the columns
		jk = 0;
		ik++;
	}//for loop for the rows

	mat uapp = e_data % (square(img2_dx) + gamma * (square(dxx) + square(dxy))) + e_sum;
	mat vapp = e_data % (square(img2_dy) + gamma * (square(dyy) + square(dxy))) + e_sum;
	mat uvapp = e_data % ((img2_dx % img2_dy) + gamma * ((dxx % dxy) + (dyy % dxy)));
	
	//4 temp matrixes used in initializing vals
	mat tmp1(height, width, fill::zeros); 
	mat tmp2 = tmp1;
	mat tmp3 = tmp1;
	mat tmp4 = tmp1;
	ik = 0;
	jk = 0;

	for (uword i = 2; i<e_height; i = i + 2) {
		for (uword j = 2; j<e_width; j = j + 2) {
			tmp1(ik, jk) = e_smooth(i - 1, j - 2); 
			tmp2(ik, jk) = e_smooth(i - 1, j);   
			tmp3(ik, jk) = e_smooth(i - 2, j - 1);
			tmp4(ik, jk) = e_smooth(i, j - 1);    

			jk++;
		}//for loop for the columns
		jk = 0;
		ik++;
	}//for loop for the rows

	vec uappv = vectorise(uapp); //turns matrix into column vector
	vec vappv = vectorise(vapp);
	vec uvappv = vectorise(uvapp);

	ik = 0;

	for (uword i = 11; i<temp4repmat; i = i + 12) {
		vals(i - 9) = uappv(ik);
		vals(i - 2) = vappv(ik);
		vals(i - 8) = uvappv(ik);
		vals(i - 3) = uvappv(ik);

		vals(i - 11) = -(tmp1(ik));
		vals(i - 5)  = -(tmp1(ik));
		vals(i - 6)  = -(tmp2(ik));
		vals(i)      = -(tmp2(ik));

		vals(i - 10) = -(tmp3(ik));
		vals(i - 4)  = -(tmp3(ik));
		vals(i - 7)  = -(tmp4(ik));
		vals(i - 1)  = -(tmp4(ik));

		ik++;
	}
	
	mat upad(height + 2, width + 2, fill::zeros);
	mat vpad(v.n_rows + 2, v.n_cols + 2, fill::zeros);

	upad.submat(1, 1, height, width) = u;
	vpad.submat(1, 1, v.n_rows, v.n_cols) = v;

	mat pdfaltsumu(height, width, fill::zeros);
	mat pdfaltsumv = pdfaltsumu;

	uvec rowstk1(e_height / 2, fill::zeros);
	uvec colstk1(e_width / 2, fill::zeros);
	uvec rowstk2(e_height / 2, fill::zeros);
	uvec colstk2(e_width / 2, fill::zeros);
	uvec rowstk3(e_height / 2, fill::zeros);
	uvec colstk3(e_width / 2, fill::zeros);
	uvec rowstk4(e_height / 2, fill::zeros);
	uvec colstk4(e_width / 2, fill::zeros);

	tmp1.resize(upad.n_rows, upad.n_cols);//tmp size upad
	tmp2 = tmp1;
	tmp3 = tmp1;
	tmp4 = tmp1;
	mat tmp5 = tmp1;

	tmp1 = upad.submat(1, 1, height, width); //(0)
	tmp2 = upad.submat(1, 0, height, width - 1);// (1)
	tmp3 = upad.submat(1, 2, height, (upad.n_cols) - 1);// (2)
	tmp4 = upad.submat(0, 1, height - 1, width);// (3)
	tmp5 = upad.submat(2, 1, (upad.n_rows) - 1, width);//(4)

	ik = 0;
	for (uword i = 2; i<e_height; i = i + 2) {
		rowstk1(ik) = i - 1;
		rowstk2(ik) = i - 1;
		rowstk3(ik) = i - 2;
		rowstk4(ik) = i;
		ik++;
	}//for loop for the rows

	jk = 0;
	for (uword j = 2; j<e_width; j = j + 2) {
		colstk1(jk) = j - 2;
		colstk2(jk) = j;
		colstk3(jk) = j - 1;
		colstk4(jk) = j - 1;
		jk++;
	}//for loop for the columns

	pdfaltsumu += e_smooth.submat(rowstk1, colstk1) % (tmp2 - tmp1);
	pdfaltsumu += e_smooth.submat(rowstk2, colstk2) % (tmp3 - tmp1);
	pdfaltsumu += e_smooth.submat(rowstk3, colstk3) % (tmp4 - tmp1);
	pdfaltsumu += e_smooth.submat(rowstk4, colstk4) % (tmp5 - tmp1);
	
	tmp1 = vpad.submat(1, 1, height, width); //(0)
	tmp2 = vpad.submat(1, 0, height, width - 1);// (1)
	tmp3 = vpad.submat(1, 2, height, (upad.n_cols) - 1);// (2)
	tmp4 = vpad.submat(0, 1, height - 1, width);// (3)
	tmp5 = vpad.submat(2, 1, (upad.n_rows) - 1, width);//(4)
								 
	//use same row and col vectors as they are the same as pdfaltsumu
	pdfaltsumv += e_smooth.submat(rowstk1, colstk1) % (tmp2 - tmp1);
	pdfaltsumv += e_smooth.submat(rowstk2, colstk2) % (tmp3 - tmp1);
	pdfaltsumv += e_smooth.submat(rowstk3, colstk3) % (tmp4 - tmp1);
	pdfaltsumv += e_smooth.submat(rowstk4, colstk4) % (tmp5 - tmp1);

	mat constu = e_data % ((img2_dx %  img_z) + gamma * ((dxx % dxz) + (dxy % dyz))) - pdfaltsumu;
	mat constv = e_data % ((img2_dy %  img_z) + gamma * ((dxy % dxz) + (dyy % dyz))) - pdfaltsumv;

	b.zeros((2 * height*width));
	vec constuv = vectorise(constu);
	vec constvv = vectorise(constv);

	ik = 0;

	for (uword i = 1; i<b.n_rows; i = i + 2) {
		b(i) = -(constvv(ik));
		b(i - 1) = -(constuv(ik));
		ik++;
	}

	ik = 0;
	uvec templ(cols.n_cols);
	for (uword i = 0; i<cols.n_cols; i = i + 1) {
		if ((cols(i) > 0) && (cols(i) < (2 * height*width + 1))) {
			templ(ik) = i;
			ik++;
		}
	}
	
	templ.resize(ik);
	for (uword i = 0; i<templ.n_rows; i = i + 1) {
		ik = templ(i);
		cols(i) = cols(ik);
		rows(i) = rows(ik);
	}

	cols.resize(templ.n_rows);
	rows.resize(templ.n_rows);
	vals = vals(templ);

	for (uword i = 0; i<cols.n_cols; i = i + 1) {
		cols(i) = cols(i) - 1;
		rows(i) = rows(i) - 1;
	}

	umat locations = join_cols(rows, cols);
	sp_mat C1(locations, vals);
	A = C1;
	
	return make_tuple(A, b);
}





tuple<urowvec, urowvec, vec, vec> build_matrix(mat img2_dx, mat img2_dy, mat img_z,
	mat dxx, mat dxy, mat dyy, mat dxz, mat dyz, mat e_data,
	mat e_smooth, mat u, mat v, double gamma) {
	sp_mat A;
	vec b;
	int height = u.n_rows;
	int width = u.n_cols;
	int e_height = e_smooth.n_rows;
	int e_width = e_smooth.n_cols;

	e_smooth.head_rows(1) = zeros<rowvec>(e_width);
	e_smooth.tail_rows(1) = zeros<rowvec>(e_width);

	e_smooth.head_cols(1) = zeros<vec>(e_height);
	e_smooth.tail_cols(1) = zeros<vec>(e_height);

	int temp4repmat = 2 * height*width * 6;
	int tempPop = 1;
	urowvec rows(temp4repmat); //column vector of size temp4repmat

	for (int i = 0; i<temp4repmat; i++) {
		if (i % 6 == 0 && i != 0) {
			tempPop++;
		}
		rows(i) = tempPop;
	}

	urowvec cols = rows;
	vec vals(temp4repmat, fill::zeros);

	for (int i = 0; i<temp4repmat; i = i + 6) {
		cols(i) = rows(i) - (2 * height);
	}

	for (int i = 1; i<temp4repmat; i = i + 6) {
		cols(i) = rows(i) - 2;
	}

	for (int i = 8; i<temp4repmat; i = i + 12) {
		cols(i) = rows(i) - 1;
	}

	for (int i = 3; i<temp4repmat; i = i + 12) {
		cols(i) = rows(i) + 1;
	}

	for (int i = 4; i<temp4repmat; i = i + 6) {
		cols(i) = rows(i) + 2;
	}

	for (int i = 5; i<temp4repmat; i = i + 6) {
		cols(i) = rows(i) + (2 * height);
	}

	mat e_sum(height, width, fill::zeros);

	int ik = 0;
	int jk = 0;

	for (int i = 2; i<e_height; i = i + 2) {
		for (int j = 2; j<e_width; j = j + 2) {
			e_sum(ik, jk) += e_smooth(i, j - 1); 
			e_sum(ik, jk) += e_smooth(i - 1, j); 
			e_sum(ik, jk) += e_smooth(i - 2, j - 1); 
			e_sum(ik, jk) += e_smooth(i - 1, j - 2); 

			jk++;
		}//for loop for the columns
		jk = 0;
		ik++;
	}//for loop for the rows

	mat uapp = e_data % (square(img2_dx) + gamma * (square(dxx) + square(dxy))) + e_sum;
	mat vapp = e_data % (square(img2_dy) + gamma * (square(dyy) + square(dxy))) + e_sum;
	mat uvapp = e_data % ((img2_dx % img2_dy) + gamma * ((dxx % dxy) + (dyy % dxy)));

	//4 temp matrixes used in initializing vals
	mat tmp1(height, width, fill::zeros); 
	mat tmp2 = tmp1;
	mat tmp3 = tmp1;
	mat tmp4 = tmp1;
	ik = 0;
	jk = 0;

	for (int i = 2; i<e_height; i = i + 2) {
		for (int j = 2; j<e_width; j = j + 2) {
			tmp1(ik, jk) = e_smooth(i - 1, j - 2); 
			tmp2(ik, jk) = e_smooth(i - 1, j);   
			tmp3(ik, jk) = e_smooth(i - 2, j - 1);
			tmp4(ik, jk) = e_smooth(i, j - 1);      

			jk++;
		}//for loop for the columns
		jk = 0;
		ik++;
	}//for loop for the rows

	vec uappv = vectorise(uapp); //turns matrix into column vector
	vec vappv = vectorise(vapp);
	vec uvappv = vectorise(uvapp);

	ik = 0;

	for (int i = 11; i<temp4repmat; i = i + 12) {
		vals(i - 9) = uappv(ik);
		vals(i - 2) = vappv(ik);
		vals(i - 8) = uvappv(ik);
		vals(i - 3) = uvappv(ik);

		vals(i - 11) = -(tmp1(ik));
		vals(i - 5) = -(tmp1(ik));
		vals(i - 6) = -(tmp2(ik));
		vals(i) = -(tmp2(ik));

		vals(i - 10) = -(tmp3(ik));
		vals(i - 4) = -(tmp3(ik));
		vals(i - 7) = -(tmp4(ik));
		vals(i - 1) = -(tmp4(ik));

		ik++;
	}

	mat upad(height + 2, width + 2, fill::zeros);
	mat vpad(v.n_rows + 2, v.n_cols + 2, fill::zeros);

	upad.submat(1, 1, height, width) = u;
	vpad.submat(1, 1, v.n_rows, v.n_cols) = v;

	mat pdfaltsumu(height, width, fill::zeros);
	mat pdfaltsumv = pdfaltsumu;

	uvec rowstk1(e_height / 2, fill::zeros);
	uvec colstk1(e_width / 2, fill::zeros);
	uvec rowstk2(e_height / 2, fill::zeros);
	uvec colstk2(e_width / 2, fill::zeros);
	uvec rowstk3(e_height / 2, fill::zeros);
	uvec colstk3(e_width / 2, fill::zeros);
	uvec rowstk4(e_height / 2, fill::zeros);
	uvec colstk4(e_width / 2, fill::zeros);

	tmp1.resize(upad.n_rows, upad.n_cols);
	tmp2 = tmp1;
	tmp3 = tmp1;
	tmp4 = tmp1;
	mat tmp5 = tmp1;

	tmp1 = upad.submat(1, 1, height, width); //(0)
	tmp2 = upad.submat(1, 0, height, width - 1);// (1)
	tmp3 = upad.submat(1, 2, height, (upad.n_cols) - 1);// (2)
	tmp4 = upad.submat(0, 1, height - 1, width);// (3)
	tmp5 = upad.submat(2, 1, (upad.n_rows) - 1, width);//(4)

	ik = 0;
	for (int i = 2; i<e_height; i = i + 2) {
		rowstk1(ik) = i - 1;
		rowstk2(ik) = i - 1;
		rowstk3(ik) = i - 2;
		rowstk4(ik) = i;
		ik++;
	}//for loop for the rows

	jk = 0;
	for (int j = 2; j<e_width; j = j + 2) {
		colstk1(jk) = j - 2;
		colstk2(jk) = j;
		colstk3(jk) = j - 1;
		colstk4(jk) = j - 1;
		jk++;
	}//for loop for the columns

	pdfaltsumu += e_smooth.submat(rowstk1, colstk1) % (tmp2 - tmp1);
	pdfaltsumu += e_smooth.submat(rowstk2, colstk2) % (tmp3 - tmp1);
	pdfaltsumu += e_smooth.submat(rowstk3, colstk3) % (tmp4 - tmp1);
	pdfaltsumu += e_smooth.submat(rowstk4, colstk4) % (tmp5 - tmp1);

	tmp1 = vpad.submat(1, 1, height, width); 
	tmp2 = vpad.submat(1, 0, height, width - 1);
	tmp3 = vpad.submat(1, 2, height, (upad.n_cols) - 1);
	tmp4 = vpad.submat(0, 1, height - 1, width);
	tmp5 = vpad.submat(2, 1, (upad.n_rows) - 1, width);
													  
	pdfaltsumv += e_smooth.submat(rowstk1, colstk1) % (tmp2 - tmp1);
	pdfaltsumv += e_smooth.submat(rowstk2, colstk2) % (tmp3 - tmp1);
	pdfaltsumv += e_smooth.submat(rowstk3, colstk3) % (tmp4 - tmp1);
	pdfaltsumv += e_smooth.submat(rowstk4, colstk4) % (tmp5 - tmp1);

	mat constu = e_data % ((img2_dx %  img_z) + gamma * ((dxx % dxz) + (dxy % dyz))) - pdfaltsumu;
	mat constv = e_data % ((img2_dy %  img_z) + gamma * ((dxy % dxz) + (dyy % dyz))) - pdfaltsumv;

	b.zeros((2 * height*width));
	vec constuv = vectorise(constu);
	vec constvv = vectorise(constv);

	ik = 0;

	for (int i = 1; i<b.n_rows; i = i + 2) {
		b(i) = -(constvv(ik));
		b(i - 1) = -(constuv(ik));
		ik++;
	}

	ik = 0;
	uvec templ(cols.n_cols);
	for (int i = 0; i<cols.n_cols; i = i + 1) {
		if ((cols(i) > 0) && (cols(i) < (2 * height*width + 1))) {
			templ(ik) = i;
			ik++;
		}
	}
	templ.resize(ik);
	for (int i = 0; i<templ.n_rows; i = i + 1) {
		ik = templ(i);
		cols(i) = cols(ik);
		rows(i) = rows(ik);
	}

	cols.resize(templ.n_rows);
	rows.resize(templ.n_rows);
	vals = vals(templ);

	for (int i = 0; i<cols.n_cols; i = i + 1) {
		cols(i) = cols(i) - 1;
		rows(i) = rows(i) - 1;
	}

	return make_tuple(rows, cols, vals, b);
}
