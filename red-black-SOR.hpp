/*
energy_calc.hpp
Jacob Wiebe & James Dolman
Rev1: Nov 2017
Rev2: Jan 2018
*/

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

//A: sor
//the solution 'x' is the vector 'duv'
//failure flag can also be used outside this function
tuple<vec, bool> redblack_sor(mat A, vec x, vec b, double omega,
	int inner_iter, double tolerance) {

	cout << "performing Red / Black SOR >" << endl;

	double error = 0.0;
	bool failure = true;

	mat U(A.n_rows, A.n_cols, fill::zeros);
	mat L(A.n_rows, A.n_cols, fill::zeros);
	//mat A(A.n_rows + 2, A.n_cols + 2, fill::zeros);
	//mat A(500, 500, fill::zeros);
	double approx = 0.01;
	//A = (A * x) + b;

	//from matlab
	double norml = norm(b);
	if (norml == 0) {
		norml = 1;
	}

	//RED LOOP

	for (int iter = 0; iter <= inner_iter; iter++) {
		cout << "Complete iteration " << iter << endl;

		//RED odd cells
		for (uword i = 2; i < A.n_rows; i += 2) { //start at first odd row, step by 2
			for (uword j = 2; j < A.n_cols; j += 2) { //start at first odd col, step by 2

				approx = (A.at(i - 1, j) + A.at(i, j) + A.at(i - 2, j - 1)
					+ A.at(i - 1, j - 2)) / 4;

				error = max(abs(omega * (approx - A.at(i - 1, j - 1))), error);

				//trailing loop
				U.at(i - 1, j - 1) = ((1 - omega) * A.at(i - 1, j - 1))
					+ (omega * approx); //Sr
			}
		}

		//RED even cells
		for (uword i = 2; i < A.n_rows; i += 2) { //start at first even row, step by 2
			for (uword j = 2; j < A.n_cols; j += 2) { //start at first even col, step by 2

				approx = (A.at(i - 1, j) + A.at(i, j) + A.at(i - 2, j - 1)
					+ A.at(i - 1, j - 2)) / 4;

				error = max(abs(omega * (approx - A.at(i - 1, j - 1))), error);

				//trailing loop
				U.at(i - 1, j - 1) = ((1 - omega) * A.at(i - 1, j - 1))
					+ (omega * approx);
			}
		}

		//BLACK LOOP

		//BLACK odd cells
		for (uword i = 2; i < A.n_rows; i += 2) { //start at first odd row, step by 2
			for (uword j = 3; j < A.n_cols; j += 2) { //start at second even col, step by 2

				approx = (A.at(i - 1, j) + A.at(i, j - 1) + A.at(i - 2, j - 1)
					+ A.at(i - 1, j - 2)) / 4;

				error = max(abs(omega * (approx - A.at(i - 1, j - 1))), error);

				//trailing loop
				L.at(i - 1, j - 1) = ((1 - omega) * A.at(i - 1, j - 1))
					+ (omega * approx);
			}
		}

		//BLACK even cells
		for (uword i = 3; i < A.n_rows; i += 2) { //start at first even row, step by 2
			for (uword j = 2; j < A.n_cols; j += 2) { //start at first even col, step by 2
													  //to avoid out of bounds erros, both loops are start ahead by 1

				approx = (A.at(i - 1, j) + A.at(i, j - 1) + A.at(i - 2, j - 1)
					+ A.at(i - 1, j - 2)) / 4;

				error = max(abs(omega * (approx - A.at(i - 1, j - 1))), error);

				//trailing loop
				L.at(i - 1, j - 1) = ((1 - omega) * A.at(i - 1, j - 1))
					+ (omega * approx);
				//printf("L is %.01f\n", L.at(i - 1, j - 1));
			}
		}

		cout << "Error is:" <<error<< endl;
		//if (error < tolerance) {
			cout << "Error is under tolerance level"<< endl;
			mat diagA = diagmat(A);


			//black cell matrix
			L = (diagA - omega * L);
			cout << "BLACK: " << L << endl;

			//red cell matrix
			U = ((omega * U) + ((1 - omega) * diagA)) * x;
			cout << "    RED: " << U;

			x = inv(L) * (U + 1);
			cout << "THIS: "<< x.at(0) << endl;
			//x = L * (U + 1);

			//failure = false;
			//break; //convergence reached
		//}


	}

	return make_tuple(x, failure);
}