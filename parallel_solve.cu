#include <stdlib.h>  
#include <stdio.h>
#include <algorithm> 
#include <cstdlib>  
#include <cuda_runtime.h>

//used for norm function
//#include <cublas.h>
#include <cublas_v2.h>

//used for conversion from coo to csr sparse matrix formats
#include <cusparse_v2.h>

#define BLOCKS 512
#define THREADS 1024
#define THREAD_INC 32


__global__ void kernel_solve_v0(double *val, int *row, int *col, double* x0, double* sum, double* aii, unsigned int ix, unsigned int n) {
	*sum = 0.0;
	
	for (int ir = 0; ir < n; ir++) {
		
		if (row[ir] == ix) {
			if (col[ir] == ix) {
				*aii = val[ir];
				//printf("ix is: %d   col: %d  \n", ix, col[ir]);
			}
			else {
				*sum += (val[ir] * x0[col[ir]]); //sum non-zero values in row vec
				//printf("sum is: %f  ", *sum);
			}
		}
	}
}

__global__ void kernel_solve_v1(double *val, int *row, int *col, double* x0, double* sum, double* aii, unsigned int ix, unsigned int n) {
	extern __shared__ double thread_sum[]; //access dynamic shared mem


	//initialize thread_sum;
	for (int ir = threadIdx.x; ir < n; ir = ir + blockDim.x) { thread_sum[ir] = 0; }

	__syncthreads();


	for (int ir = threadIdx.x; ir < n; ir = ir + blockDim.x) {
		if (row[ir] == ix) {

			if (col[ir] == ix) {
				*aii = val[ir];
			}
			else
			{
				thread_sum[ir] = ((val[ir]) * x0[col[ir]]); //sum non-zero values in row vec
				//printf("this is row %0.f at thread %d: val at index is %f\n", row[ir], ir, val[ir]);
			}
		}
	}

	__syncthreads();

	//sum all threads and save them to the sum pointer
	if (threadIdx.x == 0) {
		double total = 0.0;
		for (int i = 0; i < n; i++) {
			total += thread_sum[i];
		}
		*sum = total; //return result
	}
}

__global__ void kernel_solve_v2(double *val, int *row, int *col, double* b, double* x0, double* x, double* xdiff, unsigned int n, unsigned int m, cublasHandle_t blas_handle, bool *end_flag) {
	double sum = 0.0;
	double aii = 0.0;
	int indx = 0, next = 0;

	//for (int ix = threadIdx.x; ix < m; ix = ix + blockDim.x) {
	for (int ix = 0; ix < m; ix++) {

			indx = row[ix];
			next = row[ix + 1];
			//printf("%d:%d      ", ix, indx);
			//printf("at x %d is %d\n", ix, col[indx]);
			sum = 0.0;
			while (indx < next) {
				if (col[indx] == ix) { //condition is never met
					aii = val[indx];
					//printf("here sum shuld be 0 %d\n", sum);
				}
				else
				{
					sum += ((val[indx]) * x0[col[indx]]); //sum non-zero values in row vec
														  //printf("AT x %d, sum is %d\n", ix, sum);
				}
				indx++;
			}

			x[ix] = (b[ix] - sum) / aii;
			xdiff[ix] = x[ix] - x0[ix];
			printf("Jacobi iteration %d - error is: %f\r", ix, xdiff[ix]); //ADD THIS BACK IN
			//if (ix < 10) {
			//	printf("row: %d, x is: %f, sum is: %d\n", ix, x[ix], sum);
			//}
		} //end for (vals of x)

	//if (threadIdx.x == 0) {
	//	
	//	//compute the norm of vec x
		//cublasDnrm2(blas_handle, m, x, m, x_norm);
	//	//compute the norm of vec x_diff
	//	cublasDnrm2(blas_handle, m, xdiff, m, x_diffnorm);

	//	if (*x_diffnorm / *x_norm < 1e-8) {
	//		*end_flag = true;
	//		printf("CONVERGENCE");
	//		return;
	//	}
	//}
	
	
}

__global__ void init_coo(long long *lrow, long long *lcol, int* row, int* col, const unsigned int n) {
	for (int i = threadIdx.x; i < n; i = i + blockDim.x) {

		row[i] = (int)lrow[i];
		col[i] = (int)lcol[i];
		//if (i < 20) {
		//	printf("%d, ", row[i]);
		//}
	}
}


int* convert_csr(cusparseHandle_t sparse_handle, const int *cooRowInd, int *csrRowPtr, int n, int m) {
	
	int *d_rowPtr;

	//cudaMalloc((void **)&d_rowPtr, sizeof(int) * m+1);
	//cudaMemcpy(d_rowPtr, csrRowPtr, sizeof(int) * m+1, cudaMemcpyHostToDevice);


	//for (int i = 0; i < 10; i++) {
	//	printf("%d\n", cooRowInd[i]);
	//	//csrRowPtr[i] = 1;
	//}
	cusparseStatus_t s;
	s = cusparseXcoo2csr(sparse_handle, cooRowInd, n, m, csrRowPtr, CUSPARSE_INDEX_BASE_ZERO);
	cudaMemcpy(csrRowPtr, d_rowPtr, sizeof(double) * m+1, cudaMemcpyDeviceToHost);
	//printf("%d\n", *csrRowPtr);
	//if (s == CUSPARSE_STATUS_SUCCESS) {
	//	for (int i = 0; i < 20; i++) {
	//		printf("%d, ", csrRowPtr[i]);
	//	}
	//}
	return csrRowPtr;
}

int determine_threads(int n, int m) {
	int threads = n / THREAD_INC;
	if ((threads*THREAD_INC) + THREAD_INC <= 1024) {
		threads = (threads*THREAD_INC) + THREAD_INC;
	}
	else {
		threads = 1024;
	}

	if (n * sizeof(double) > 48 * 1024) {
		printf("IMAGE TOO LARGE > EXITING");
		exit(EXIT_FAILURE);
	}

	return threads;
}

double compute_error(double* x, double* x_diff, double* d_x, double* d_xdiff, int x_size, cublasHandle_t blas_handle) {
	double* x_norm = (double*)calloc(1, sizeof(double));
	double* x_diffnorm = (double*)calloc(1, sizeof(double));

	cudaMemcpy(d_x, x, (sizeof(double) * x_size), cudaMemcpyHostToDevice);
	cudaMemcpy(d_xdiff, x_diff, (sizeof(double) * x_size), cudaMemcpyHostToDevice);

	//compute the norm of vec x
	cublasDnrm2(blas_handle, x_size, d_x, x_size, x_norm);
	//compute the norm of vec x_diff
	cublasDnrm2(blas_handle, x_size, d_xdiff, x_size, x_diffnorm);

	return *x_diffnorm / *x_norm;
}

double* parallel_solve_v2(double *val, long long *lrow, long long *lcol, double *b, const unsigned int n, const unsigned int x_size, const unsigned int max_iter, double tol) {
	//for norm function
	cublasHandle_t blas_handle;
	cusparseHandle_t sparse_handle;
	cusparseMatDescr_t descr = 0;

	cusparseCreate(&sparse_handle);
	cublasCreate(&blas_handle);
	cusparseCreateMatDescr(&descr);

	cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);

	//allocate memory of local variables
	bool *end_flag = (bool*)calloc(1, sizeof(bool));
	int *row = (int*)calloc(n, sizeof(int));
	int *col = (int*)calloc(n, sizeof(int));
	//double* sum = (double*)calloc(1, sizeof(double));
	//double* aii = (double*)calloc(1, sizeof(double));
	double* x = (double*)calloc(x_size, sizeof(double));
	double* x0 = (double*)calloc(x_size, sizeof(double));
	double* x_diff = (double*)calloc(x_size, sizeof(double));
	double* x_norm = (double*)calloc(1, sizeof(double));
	double* x_diffnorm = (double*)calloc(1, sizeof(double));

	//declare pointers to device memory
	bool *d_flag;
	long long *d_lrow, *d_lcol;
	double *d_val, *d_x0, *d_b, *d_x, *d_xdiff;
	int *d_row, *d_col;

	//allocate device memory
	cudaMalloc((void **)&d_flag, sizeof(double));
	cudaMalloc((void **)&d_val, sizeof(double) * n);
	cudaMalloc((void **)&d_row, sizeof(int) * n);
	cudaMalloc((void **)&d_col, sizeof(int) * n);
	cudaMalloc((void **)&d_lrow, sizeof(long long) * n);
	cudaMalloc((void **)&d_lcol, sizeof(long long) * n);
	cudaMalloc((void **)&d_x0, sizeof(double) * x_size);
	//cudaMalloc((void **)&d_sum, sizeof(double));
	//cudaMalloc((void **)&d_aii, sizeof(double));
	cudaMalloc((void **)&d_b, sizeof(double) * x_size);
	cudaMalloc((void **)&d_x, sizeof(double) * x_size);
	cudaMalloc((void **)&d_xdiff, sizeof(double) * x_size);

	//copy sparse matrix onto device (INPUTS)
	cudaMemcpy(d_val, val, (sizeof(double) * n), cudaMemcpyHostToDevice);
	cudaMemcpy(d_lcol, lcol, (sizeof(long long) * n), cudaMemcpyHostToDevice);
	cudaMemcpy(d_lrow, lrow, (sizeof(long long) * n), cudaMemcpyHostToDevice);
	cudaMemcpy(d_b, b, (sizeof(double) * x_size), cudaMemcpyHostToDevice);
	//cudaMemcpy(d_x, x, (sizeof(double) * n), cudaMemcpyHostToDevice);
	//cudaMemcpy(d_xdiff, x_diff, (sizeof(double) * n), cudaMemcpyHostToDevice);

	//determine optimal number of threads while preventing overflow
	//printf("Size of X: %d\n", x_size);
	int threads = determine_threads(x_size, x_size);
	//int threads = x_size; //what do you think Jacob?
	printf("Solving using parallel Jacobi method\nCurrent matrix contains %d coefficients, solving for %d values\n%d threads assigned\n", n, x_size, threads);
	//printf("threads: %d\n", threads);

	//convert input row & col vectors into integer vectors
	init_coo << <1, n >> > (d_lrow, d_lcol, d_row, d_col, n);


	//this is just for testing
	//cudaMemcpy(row, d_row, (sizeof(int) * n), cudaMemcpyDeviceToHost);
	//cudaMemcpy(col, d_col, (sizeof(int) * n), cudaMemcpyDeviceToHost);

	cudaDeviceSynchronize();


	//int *csrRowPtr = (int*)calloc(x_size + 1, sizeof(int));
	int *d_rowPtr = 0;

	//for (int i = 0; i < (x_size + 1); i++) {
	//	printf(":%d %d\n", row[i],i);
	//}

	cudaMalloc((void**)&d_rowPtr, (x_size + 1) * sizeof(int));
	cusparseXcoo2csr(sparse_handle, d_row, n, x_size, d_rowPtr, CUSPARSE_INDEX_BASE_ZERO);
	cudaMemcpy(row, d_rowPtr, (sizeof(int) * (x_size + 1)), cudaMemcpyDeviceToHost);

	//for (int i = 0; i < (x_size + 1); i++) {
	//	printf(":%d\n", row[i]);
	//}

	//row = convert_csr(sparse_handle, row, n, x_size);
	cudaFree(d_lrow);
	cudaFree(d_lcol);
	cudaFree(d_rowPtr);

	for (int i = 0; i < max_iter; i++) {
		*end_flag = false;
		//printf("x at iter %d is %f : %f\n", i, *x, *x0);

		std::swap(x0, x);
		cudaMemcpy(d_x0, x0, (sizeof(double) * x_size), cudaMemcpyHostToDevice);
		
		kernel_solve_v2 << <1, threads>> > (d_val, d_row, d_col, d_b, d_x0, d_x, d_xdiff, n, x_size, blas_handle, end_flag);

		cudaMemcpy(x, d_x, sizeof(double) * x_size, cudaMemcpyDeviceToHost);
		cudaMemcpy(end_flag, d_flag, sizeof(bool), cudaMemcpyDeviceToHost);

		if (*end_flag) { break; };
	} //end for (iterations)
	printf("\n");

	  //wait for all processes to complete
	cudaDeviceSynchronize();

	cublasDestroy(blas_handle);

	//free device memory
	cudaFree(d_val);
	cudaFree(d_row);
	cudaFree(d_col);
	cudaFree(d_x0);
	cudaFree(d_x);
	cudaFree(d_xdiff);
	//cudaFree(d_sum);
	//cudaFree(d_aii);

	return x;
}

double* parallel_solve_v1(double *val, long long *lrow, long long *lcol, double *b, int *fail, const unsigned int n, const unsigned int x_size, const unsigned int max_iter, double tol) {
	//for norm function
	cublasHandle_t blas_handle;
	cublasCreate(&blas_handle);


	double* x_norm = (double*)calloc(1, sizeof(double));
	double* x_diffnorm = (double*)calloc(1, sizeof(double));

	//allocate memory of local variables
	int *row = (int*)calloc(n, sizeof(int));
	int *col = (int*)calloc(n, sizeof(int));
	double* sum = (double*)calloc(1, sizeof(double));
	double* aii = (double*)calloc(1, sizeof(double));
	double* x = (double*)calloc(x_size, sizeof(double));
	double* x0 = (double*)calloc(x_size, sizeof(double));
	double* x_diff = (double*)calloc(x_size, sizeof(double));

	//declare pointers to device memory
	double *d_val, *d_x0, *d_sum, *d_aii, *d_x, *d_xdiff;
	int *d_row, *d_col;

	//allocate device memory
	cudaMalloc((void **)&d_val, sizeof(double) * n);
	cudaMalloc((void **)&d_row, sizeof(int) * n);
	cudaMalloc((void **)&d_col, sizeof(int) * n);
	cudaMalloc((void **)&d_x0, sizeof(double) * x_size);
	cudaMalloc((void **)&d_sum, sizeof(double));
	cudaMalloc((void **)&d_aii, sizeof(double));
	cudaMalloc((void **)&d_x, sizeof(double) * x_size);
	cudaMalloc((void **)&d_xdiff, sizeof(double) * x_size);


	for (int i = 0; i < n; i++) {
		row[i] = (int)lrow[i];
		col[i] = (int)lcol[i];
		//printf("col: %d\n", row[i]);
	}

	//copy sparse matrix onto device (INPUTS)
	cudaMemcpy(d_val, val, (sizeof(double) * n), cudaMemcpyHostToDevice);
	cudaMemcpy(d_col, col, (sizeof(int) * n), cudaMemcpyHostToDevice);
	cudaMemcpy(d_row, row, (sizeof(int) * n), cudaMemcpyHostToDevice);

	//determine optimal number of threads while preventing overflow
	int threads = determine_threads(n, x_size);
	printf("Solving using parallel Jacobi method\nCurrent matrix contains %d coefficients, solving for %d values\n%d threads assigned\n", n, x_size, threads);

	//convert input row & col vectors into integer vectors
	//init_coo << <1, threads >> > ( d_lrow, d_lcol, d_row, d_col, n);
	

	//cudaMemcpy(row, d_row, (sizeof(int) * n), cudaMemcpyDeviceToHost);
	//cudaMemcpy(col, d_col, (sizeof(int) * n), cudaMemcpyDeviceToHost);

	////cudaDeviceSynchronize();

	//for (int ib = 0; ib < 20; ib++) {
	//	printf("%d, ", row[ib]);
	//}

	//int *csrRowPtr = 0;

	//csrRowPtr = (int*)calloc(x_size + 1, sizeof(int));

	//
	//cusparseXcoo2csr(sparse_handle, d_row, n, x_size, d_rowPtr, CUSPARSE_INDEX_BASE_ZERO);
	//cudaMemcpy(csrRowPtr, d_rowPtr, (sizeof(int) * (x_size + 1)), cudaMemcpyDeviceToHost);


	
	//cusparseXcoo2csr(sparse_handle, d_row, n, x_size, csrRowPtr, CUSPARSE_INDEX_BASE_ZERO);
	//csrRowPtr = convert_csr(sparse_handle, row, csrRowPtr, n, x_size);

	//for (int ib = 0; ib < 20; ib++) {
	//	printf("(%d, ", col[ib]);
	//	printf("%d)", row[ib]);
	//}
	//cudaFree(d_lrow);
	//cudaFree(d_lcol);


	/* Here is full program that is set for parallel running. It always produces sum = 0. There may be an issue with copying 
	mem from device in a loop. We may need to re-copy back to device for each iteration. currently testing with full-fat but
	I need a testing strategy to break it up and check the individual components.*/

	/*At the 3rd iteration, (x=2) the graphics crashes and the outputs start being wrong (I think the actual kernel stops
	executing and the loop continues printing and updating x[ix]*/

	/*first 2 iterations of x are correct. From there x goes off some threshold. This is a memory problem. I need some way to
	free up x0 memory when I'm done each iteration because its eating up all the space.*/

	/*Results vary from totally correct to compeletely random with every run. I have no idea what is causing this.*/

	/*Problem was the order that the memcpy happened*/
	
	for (int i = 0; i < max_iter; i++) {
		//printf("x at iter %d is %f : %f\n", i, *x, *x0);

		//x0 = x;
		std::swap(x0, x);
		cudaMemcpy(d_x0, x0, (sizeof(double) * x_size), cudaMemcpyHostToDevice);

		for (int ix = 0; ix < x_size; ix++) {
			*sum = 0.0;

				//execute solve in parallel and allocate n doubles of shared memory
				kernel_solve_v1 << <1, threads, n * sizeof(double) >> > (d_val, d_row, d_col, d_x0, d_sum, d_aii, ix, n);

				//copy shared memory to host (OUTPUTS)
				cudaMemcpy(sum, d_sum, sizeof(double), cudaMemcpyDeviceToHost);
				cudaMemcpy(aii, d_aii, sizeof(double), cudaMemcpyDeviceToHost);
				
			x[ix] = (b[ix] - *sum) / *aii;
			x_diff[ix] = x[ix] - x0[ix];
			if (ix < 10) {
				//printf("row: %d, x is: %f:%f, sum is: %f\n", ix, x[ix], x0[ix], *sum);
			}
		} //end for (vals of x)

		//calculate error between x and x_diff
		double error = compute_error(x, x_diff, d_x, d_xdiff, x_size, blas_handle);
		printf("Jacobi iteration %d - error is: %f\r", i, error); //ADD THIS BACK IN
		////printf("iteration: %d error: %f\n", i, error);
		if (error <= tol) {
			*fail = 0;
			printf("Convergence reached at iter %d, val %f\n", i, *x);
			break; //convergence reached
		}
	} //end for (iterations)

	//wait for all processes to complete
	cudaDeviceSynchronize();

	cublasDestroy(blas_handle);

	//free device memory
	cudaFree(d_val);
	cudaFree(d_row);
	cudaFree(d_col);
	cudaFree(d_x0);
	cudaFree(d_x);
	cudaFree(d_xdiff);
	cudaFree(d_sum);
	cudaFree(d_aii);

	return x;
}

double* jacobi_serial(double *val, long long *lrow, long long *lcol, double *b, const unsigned int n, const unsigned int x_size, const unsigned int max_iter) {
	//allocate output vectors x and x0
	double *x = new double[x_size];
	double *x0 = new double[x_size];
	double sum = 0.0, aii = 0.0;

	int *row = (int*)calloc(n, sizeof(int));
	int *col = (int*)calloc(n, sizeof(int));

	for (int i = 0; i < n; i++) {
		row[i] = (int)lrow[i];
		col[i] = (int)lcol[i];
		//printf("%d\n", row[i]);
	}

	for (int i = 0; i < max_iter; i++) {
		x0 = x;
		//printf("x at iter %d is %f\n", i, *x);

		for (int ix = 0; ix < x_size; ix++) {
			sum = 0; //reset sum
			
			for (int ir = 0; ir < n; ir++) {
				if (row[ir] == ix) {
					if (col[ir] == ix) {
						aii = val[ir];
					}
					else {
						sum += (val[ir] * x0[col[ir]]); //sum non-zero values in row vec
					}
				}
			}
			x[ix] = (b[ix] - sum) / aii; //update position of x
		}
		printf("Jacobi iteration %d - x[0] is: %1.9f\r", i, x[0]); //ADD THIS BACK IN
	}
	
	return x;
}

double* jacobi_serial_v2(double *val, long long *lrow, long long *lcol, double *b, const unsigned int n, const unsigned int x_size, const unsigned int max_iter) {
	double sum = 0.0;
	double aii = 0.0;
	int indx = 0, next = 0;

	//cusparseHandle_t sparse_handle;

	//double *x = new double[x_size];
	//double *x0 = new double[x_size];
	//double *xdiff = new double[x_size];
	//double sum = 0.0;
	//double aii = 0.0;
	//int indx = 0, next = 0;

	//for norm function
	cublasHandle_t blas_handle;
	cusparseHandle_t sparse_handle;
	cusparseMatDescr_t descr = 0;

	cusparseCreate(&sparse_handle);
	cublasCreate(&blas_handle);
	cusparseCreateMatDescr(&descr);

	cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);

	//allocate memory of local variables
	bool *end_flag = (bool*)calloc(1, sizeof(bool));
	int *row = (int*)calloc(n, sizeof(int));
	int *col = (int*)calloc(n, sizeof(int));
	//double* sum = (double*)calloc(1, sizeof(double));
	//double* aii = (double*)calloc(1, sizeof(double));
	double* x = (double*)calloc(x_size, sizeof(double));
	double* x0 = (double*)calloc(x_size, sizeof(double));
	double* x_diff = (double*)calloc(x_size, sizeof(double));
	double* x_norm = (double*)calloc(1, sizeof(double));
	double* x_diffnorm = (double*)calloc(1, sizeof(double));

	//declare pointers to device memory
	bool *d_flag;
	long long *d_lrow, *d_lcol;
	double *d_val, *d_x0, *d_b, *d_x, *d_xdiff;
	int *d_row, *d_col;

	//allocate device memory
	cudaMalloc((void **)&d_flag, sizeof(double));
	cudaMalloc((void **)&d_val, sizeof(double) * n);
	cudaMalloc((void **)&d_row, sizeof(int) * n);
	cudaMalloc((void **)&d_col, sizeof(int) * n);
	cudaMalloc((void **)&d_lrow, sizeof(long long) * n);
	cudaMalloc((void **)&d_lcol, sizeof(long long) * n);
	cudaMalloc((void **)&d_x0, sizeof(double) * x_size);
	//cudaMalloc((void **)&d_sum, sizeof(double));
	//cudaMalloc((void **)&d_aii, sizeof(double));
	cudaMalloc((void **)&d_b, sizeof(double) * x_size);
	cudaMalloc((void **)&d_x, sizeof(double) * x_size);
	cudaMalloc((void **)&d_xdiff, sizeof(double) * x_size);

	//copy sparse matrix onto device (INPUTS)
	cudaMemcpy(d_val, val, (sizeof(double) * n), cudaMemcpyHostToDevice);
	cudaMemcpy(d_lcol, lcol, (sizeof(long long) * n), cudaMemcpyHostToDevice);
	cudaMemcpy(d_lrow, lrow, (sizeof(long long) * n), cudaMemcpyHostToDevice);
	cudaMemcpy(d_b, b, (sizeof(double) * x_size), cudaMemcpyHostToDevice);
	//cudaMemcpy(d_x, x, (sizeof(double) * n), cudaMemcpyHostToDevice);
	//cudaMemcpy(d_xdiff, x_diff, (sizeof(double) * n), cudaMemcpyHostToDevice);

	//determine optimal number of threads while preventing overflow
	printf("Size of X: %d\n", x_size);
	int threads = determine_threads(x_size, x_size);
	//int threads = x_size; //what do you think Jacob?
	printf("threads: %d\n", threads);

	//convert input row & col vectors into integer vectors
	init_coo << <1, n >> > (d_lrow, d_lcol, d_row, d_col, n);


	//this is just for testing
	//cudaMemcpy(row, d_row, (sizeof(int) * n), cudaMemcpyDeviceToHost);
	cudaMemcpy(col, d_col, (sizeof(int) * n), cudaMemcpyDeviceToHost); //serial only!

	cudaDeviceSynchronize();


	//int *csrRowPtr = (int*)calloc(x_size + 1, sizeof(int));
	int *d_rowPtr = 0;

	//for (int i = 0; i < (x_size + 1); i++) {
	//	printf(":%d %d\n", row[i],i);
	//}

	cudaMalloc((void**)&d_rowPtr, (x_size + 1) * sizeof(int));
	cusparseXcoo2csr(sparse_handle, d_row, n, x_size, d_rowPtr, CUSPARSE_INDEX_BASE_ZERO);
	cudaMemcpy(row, d_rowPtr, (sizeof(int) * (x_size + 1)), cudaMemcpyDeviceToHost);

	//row = convert_csr(sparse_handle, row, n, x_size);
	cudaFree(d_lrow);
	cudaFree(d_lcol);
	cudaFree(d_rowPtr);

	//for (int i = 0; i < (x_size + 1); i++) {
	//	printf(":%d col:%d at %d\n", row[i], col[i], i);
	//}
	for (int i = 0; i < 2; i++) {
		*end_flag = false;
		//printf("x at iter %d is %f : %f\n", i, *x, *x0);

		std::swap(x0, x);
		

		//THIS IS KERNEL CALL CODE///////////////////////////////////////////////////
		for (int ix = 0; ix < x_size; ix++) {
			
			indx = row[ix];
			next = row[ix + 1];
			
			//printf("at x %d is %d\n", ix, col[indx]);
			sum = 0.0;
			while (indx < next) {
				printf("%d:%d      ", ix, indx);
				if (col[indx] == ix) { //condition is never met
					aii = val[indx];
					//printf("here sum shuld be 0 %d\n", sum);
				}
				else
				{
					sum += ((val[indx]) * x0[col[indx]]); //sum non-zero values in row vec
					//printf("AT x %d, sum is %d\n", ix, sum);
				}
				indx++;
			}

			x[ix] = (b[ix] - sum) / aii;
			x_diff[ix] = x[ix] - x0[ix];
			//if (ix < 10) {
			//	printf("row: %d, x is: %f, sum is: %d\n", ix, x[ix], sum);
			//}
		} //end for (vals of x)
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////
	  //wait for all processes to complete
	cudaDeviceSynchronize();

	cublasDestroy(blas_handle);

	//free device memory
	cudaFree(d_val);
	cudaFree(d_row);
	cudaFree(d_col);
	cudaFree(d_x0);
	cudaFree(d_x);
	cudaFree(d_xdiff);

	return x;
}