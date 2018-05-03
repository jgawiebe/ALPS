
#include <iostream>
#include <fstream>

#include <armadillo>
#include "matrix_builder.hpp"
#include "successive_overrelaxation.hpp"

using namespace std;
using namespace arma;
void init_variables();

	//variables for matrix_builder
    sp_mat A;//inspect matlab code, see if these are edited
	vec b;
	urowvec lrow, lcol;
	vec val;
	mat img2_dx;
	mat img2_dy;
	mat img_z;
	mat dxx;
	mat dxy;
	mat dyy;
	mat dxz;
	mat dyz;
    mat e_data;
	mat e_smooth;
	mat u;
	mat v;
	double gam;

	//variables for successive_overrelaxation (less A and b)
	double omega;
	double max_it;
	double tol;
	vec duv;
	int failure = 1;
	int* fail = &failure;

	unsigned int n, x_size;
	const unsigned int max_iter = 500;

	extern double* parallel_solve_v2(double *val, long long *lrow, long long *lcol, double *b, const unsigned int n, const unsigned int x_size, const unsigned int max_iter, double tol);
	extern double* parallel_solve_v1(double *val, long long *lrow, long long *lcol, double *b, int *fail, const unsigned int n, const unsigned int x_size, const unsigned int max_iter, double tol);
	double* jacobi_serial(double *val, long long *lrow, long long *lcol, double *b, const unsigned int n, const unsigned int x_size, const unsigned int max_iter);
	extern double* jacobi_serial_v2(double *val, long long *lrow, long long *lcol, double *b, const unsigned int n, const unsigned int x_size, const unsigned int max_iter);


int main() {

	init_variables();

	tie(lrow, lcol, val, b) = build_matrix(img2_dx, img2_dy, img_z, dxx, dxy, dyy, dxz, dyz, e_data, e_smooth, u, v, gam);
	//cout << "Out of build_matrix" << endl;

	n = (int)val.n_elem;
	x_size = (int)b.n_elem;

	//int *row = (int*)calloc(n, sizeof(int));
	//int *col = (int*)calloc(n, sizeof(int));

	////cast all elements of vector into int array
	////put this in a second kernel eventually
	//for (int i = 0; i < 10; i++) {
	//	//row[i] = (int)lrow[i];
	//	//col[i] = (int)lcol[i];
	//cout << lrow[i] << endl;
	//}

	//MUST SELECT ONE OF THESE FOR TESTING
	//double* x = jacobi_serial_v2(val.memptr(), (long long*)lrow.memptr(), (long long*)lcol.memptr(), b.memptr(), n, x_size, max_iter);
	//double* x = parallel_solve_v1(val.memptr(), (long long*)lrow.memptr(), (long long*)lcol.memptr(), b.memptr(), fail, n, x_size, max_iter, tol);
	//double* x = parallel_solve_v2(val.memptr(), (long long*)lrow.memptr(), (long long*)lcol.memptr(), b.memptr(), n, x_size, max_iter, tol);
	
	//vec duv(x, x_size); //may want to try the other parameters

	tie(A, duv) = build_sparse(img2_dx, img2_dy, img_z, dxx, dxy, dyy, dxz, dyz, e_data, e_smooth  , u, v, gam);

	printf("Solving using SOR method\nCurrent matrix contains %d coefficients, solving for %d values\n1 thread assigned\n", n, x_size);
	tie(duv, failure) = successive_overrelaxation(A, duv, omega, max_it, tol);
	cout<<"Out of successive_overrelaxation - hit enter to save duv"<<endl;
	cin.get();
	duv.save("duv_cuda-c.txt", raw_ascii);
	return 0;
}

void init_variables(){
	//variables loaded in for matrix_builder
	img2_dx.load("mats/SmallHori/IkxHoriSmall.txt");
	img2_dy.load("mats/SmallHori/IkyHoriSmall.txt");
	img_z.load("mats/SmallHori/IkzHoriSmall.txt");
	dxx.load("mats/SmallHori/IxxHoriSmall.txt");
	dxy.load("mats/SmallHori/IxyHoriSmall.txt");
	dyy.load("mats/SmallHori/IyyHoriSmall.txt");
	dxz.load("mats/SmallHori/IxzHoriSmall.txt");
    dyz.load("mats/SmallHori/IyzHoriSmall.txt");
	e_data.load("mats/SmallHori/E_DataHoriSmall.txt");
	e_smooth.load("mats/SmallHori/aE_smoothHoriSmall.txt");
	u.load("mats/SmallHori/uHoriSmall.txt");
	v.load("mats/SmallHori/vHoriSmall.txt");
    gam = 80; //value stored in gamma.txt

    //variables loaded in for successive_overrelaxation.
    omega = 1.8;
	max_it = 500;
	tol = 1e-8;
	//*failure = true;

}
