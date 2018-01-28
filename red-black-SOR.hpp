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

//M: sor
//the solution 'x' is the vector 'duv'
//failure flag can also be used outside this function
tuple<vec, bool> redblack_sor(bool failure, mat A, vec x, vec b, double omega,
		int inner_iter, double tolerance) {

	double error;
	failure = true;

	mat M(A);
	//mat M(A.n_rows + 2, A.n_cols + 2, fill::zeros);
	//mat M(500, 500, fill::zeros);
	double approx = 0.0;
	//M = (A * x) + b;

	for (uword i = 0; i < M.n_rows; i++) {
		for (uword j = 0; j < M.n_cols; j++) {
			if (i == 0) {
				M.at(i, j) = 1.0;
			} else {
				M.at(i, j) = approx;
			}
		}
	}

//RED LOOP
	double sum = 0.0;

	for (int iter = 0; iter <= inner_iter; iter++) {
		error = 0.0;

		//RED odd cells
		for (uword i = 2; i < M.n_rows; i += 2) { //start at first odd row, step by 2
			for (uword j = 2; j < M.n_cols; j += 2) { //start at first odd col, step by 2

				sum = (M.at(i - 1, j) + M.at(i, j) + M.at(i - 2, j - 1)
						+ M.at(i - 1, j - 2)) * 0.25;

				error = max(abs(omega * (sum - M.at(i - 1, j - 1))), error);

				//trailing loop
				M.at(i - 1, j - 1) = ((1 - omega) * M.at(i - 1, j - 1))
						+ (omega * sum);
			}
		}

		//RED even cells
		for (uword i = 2; i < M.n_rows; i += 2) { //start at first even row, step by 2
			for (uword j = 2; j < M.n_cols; j += 2) { //start at first even col, step by 2

				sum = (M.at(i - 1, j) + M.at(i, j) + M.at(i - 2, j - 1)
						+ M.at(i - 1, j - 2)) * 0.25;

				error = max(abs(omega * (sum - M.at(i - 1, j - 1))), error);

				//trailing loop
				M.at(i - 1, j - 1) = ((1 - omega) * M.at(i - 1, j - 1))
						+ (omega * sum);
			}
		}

		//BLACK LOOP

		//BLACK odd cells
		for (uword i = 2; i < M.n_rows; i += 2) { //start at first odd row, step by 2
			for (uword j = 3; j < M.n_cols; j += 2) { //start at second even col, step by 2

				sum = (M.at(i - 1, j) + M.at(i, j - 1) + M.at(i - 2, j - 1)
						+ M.at(i - 1, j - 2)) * 0.25;

				error = max(abs(omega * (sum - M.at(i - 1, j - 1))), error);

				//trailing loop
				M.at(i - 1, j - 1) = ((1 - omega) * M.at(i - 1, j - 1))
						+ (omega * sum);
			}
		}

		//BLACK even cells
		for (uword i = 3; i < M.n_rows; i += 2) { //start at first even row, step by 2
			for (uword j = 2; j < M.n_cols; j += 2) { //start at first even col, step by 2
				//to avoid out of bounds erros, both loops are start ahead by 1

				sum = (M.at(i - 1, j) + M.at(i, j - 1) + M.at(i - 2, j - 1)
						+ M.at(i - 1, j - 2)) * 0.25;

				error = max(abs(omega * (sum - M.at(i - 1, j - 1))), error);

				//trailing loop
				M.at(i - 1, j - 1) = ((1 - omega) * M.at(i - 1, j - 1))
						+ (omega * sum);
			}
		}

		if (error < tolerance) {
			failure = false;
			break; //convergence reached
		}

	}

	return make_tuple(x, failure);
}

