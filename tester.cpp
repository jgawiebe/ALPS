#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

vec val = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
vec row_ix = { 0, 0, 1, 1, 2, 6, 3, 4, 5 };
vec col_ix = { 2, 4, 3, 4, 5, 0, 0, 1, 1};
vec b = { 10, 9, 8, 6, 7, 5, 4, 3, 2 };

const unsigned int n = val.n_elem, x_size = b.n_elem;
const unsigned int max_iter = 500;

extern double* serial_sor(double *val, double *row, double *col, double *b, const unsigned int n, const unsigned int x_size, const unsigned int max_iter);

int main() {

	double* x = serial_sor(val.memptr(), row_ix.memptr(), col_ix.memptr(), b.memptr(), n, x_size, max_iter);

	vec duv(x, x_size); //may want to try the other parameters

	cout << duv << endl;
	cin.get();

}
