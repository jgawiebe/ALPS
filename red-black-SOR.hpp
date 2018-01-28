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
tuple<vec, bool> redblack_sor(mat A, vec x, vec b, double omega,
		int inner_iter, double tolerance) {

	double error;
	double approx = 0, sum = 0;
	bool failure = true;

	mat M(A.n_rows + 2, A.n_cols + 2, fill::zeros);
	mat Mp(M); //M prime
	M.at(0, 0) = 1.0;

	//M = (A * x) + b;

	for (int iter = 1; iter <= inner_iter; iter++) {
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

