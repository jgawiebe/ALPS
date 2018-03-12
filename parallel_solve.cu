#include <stdlib.h>  
#include <stdio.h>
#include <cstdlib>  
#include <cuda_runtime.h>

//used for norm function
#include <cublas.h>
#include <cublas_v2.h>

#define BLOCKS 512
#define THREADS 1024


__global__ void kernel_solve1(double *val, float *row, float *col, double* x0, double* sum, double* aii, unsigned int ix, unsigned int n) {
	*sum = 0;
	for (int ir = 0; ir < n; ir++) {
		if ((int)row[ir] == ix) {
			if ((int)col[ir] == ix) {
				*aii = val[ir];
			}
			else {
				*sum += (val[ir] * x0[(int)col[ir]]); //sum non-zero values in row vec
			}
		}
	}
}


__global__ void kernel_solve(double *val, float *row, float *col, double* x0, double* sum, double* aii, unsigned int ix, unsigned int n) {
	extern __shared__ double thread_sum[]; //access dynamic shared mem
	//thread_sum = new double[n]; //dynamically allocated array of sums

	//index for each element in the row vector
	//unsigned int row_idx = threadIdx.x;


	for (int row_idx = threadIdx.x; row_idx < n; row_idx = row_idx + blockDim.x)
	{
		thread_sum[row_idx] = 0;
	}
	__syncthreads();


	for (int row_idx = threadIdx.x; row_idx < n; row_idx = row_idx + blockDim.x) {
		if ((int)row[row_idx] == ix) {

			if ((int)col[row_idx] == ix) {
				*aii = val[row_idx];
			}
			else
			{
				thread_sum[row_idx] = ((val[row_idx]) * x0[(int)col[row_idx]]); //sum non-zero values in row vec
				//printf("this is row %0.f at thread %d: val at index is %f\n", row[row_idx], row_idx, val[row_idx]);
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

double* parallel_sor(double *val, float *row, float *col, double *b, const unsigned int n, const unsigned int x_size, const unsigned int max_iter, double tol) {
	//for norm function
	cublasHandle_t handle;
	cublasCreate(&handle);

	//declare memory size for vectors
	const unsigned int val_mem = sizeof(double) * n;
	const unsigned int vec_mem = sizeof(float) * n;
	const unsigned int x_mem = sizeof(double) * x_size;
	const unsigned int single = sizeof(double);

	//declare and allocate solver outputs
	double* sum = (double*)malloc(single);
	double* aii = (double*)malloc(single);
	double* x = (double*)malloc(x_mem);
	double* x0 = (double*)calloc(x_size, sizeof(double));
	

	//move these to where they're needed??
	double* x_norm = (double*)malloc(single);
	double* x_diff = (double*)malloc(x_mem);

	//declare pointers to device memory
	double *d_val, *d_x0, *d_sum, *d_aii;
	float *d_row, *d_col;

	//allocate device memory for each vector
	cudaMalloc((void **)&d_val, sizeof(double) * n);
	cudaMalloc((void **)&d_row, vec_mem);
	cudaMalloc((void **)&d_col, vec_mem);
	cudaMalloc((void **)&d_x0, x_mem);
	cudaMalloc((void **)&d_sum, single);
	cudaMalloc((void **)&d_aii, single);

	//copy memory and link pointers (INPUTS)
	cudaMemcpy(d_val, val, val_mem, cudaMemcpyHostToDevice);
	cudaMemcpy(d_row, row, vec_mem, cudaMemcpyHostToDevice);
	cudaMemcpy(d_col, col, vec_mem, cudaMemcpyHostToDevice);
	cudaMemcpy(d_x0, x0, x_mem, cudaMemcpyHostToDevice);

	//are these needed?
	cudaMemcpy(d_aii, aii, single, cudaMemcpyHostToDevice);
	cudaMemcpy(d_sum, aii, single, cudaMemcpyHostToDevice);


	/* Here is full program that is set for parallel running. It always produces sum = 0. There may be an issue with copying 
	mem from device in a loop. We may need to re-copy back to device for each iteration. currently testing with full-fat but
	I need a testing strategy to break it up and check the individual components.*/

	/*At the 3rd iteration, (x=2) the graphics crashes and the outputs start being wrong (I think the actual kernel stops
	executing and the loop continues printing and updating x[ix]*/

	/*first 2 iterations of x are correct. From there x goes off some threshold. This is a memory problem. I need some way to
	free up x0 memory when I'm done each iteration because its eating up all the space.*/

	/*Results vary from totally correct to compeletely random with every run. I have no idea what is causing this.*/
	for (int i = 0; i < 200; i++) {
		x0 = x;
		printf("x at iter %d is %f\n", i, *x);

		if (i > 0) {
			//freeing mem doesn't seem to make a difference
			//cudaFree(d_x0);
			//cudaMalloc((void **)&d_x0, x_mem);

			cudaMemcpy(d_x0, x0, x_mem, cudaMemcpyHostToDevice);
		}
		
		for (int ix = 0; ix < x_size; ix++) {
			//*sum = 0.0;

			//execute solve in parallel and allocate n doubles of shared memory
			//check shared mem size n*sizeof(double) <= 48*1024
			//kernel_solve << <1, 128, n*sizeof(double) >> > (d_val, d_row, d_col, d_x0, d_sum, d_aii, ix, n);
			int option = 0;
			if (option == 0) {
				kernel_solve1 << <1, 1 >> > (d_val, d_row, d_col, d_x0, d_sum, d_aii, ix, n);
				//copy shared memory to host (OUTPUTS)
				cudaMemcpy(sum, d_sum, single, cudaMemcpyDeviceToHost);
				cudaMemcpy(aii, d_aii, single, cudaMemcpyDeviceToHost);
			}
			else
			{
				*sum = 0; //reset sum

				for (int ir = 0; ir < n; ir++) {
					if ((int)row[ir] == ix) {
						if ((int)col[ir] == ix) {
							*aii = val[ir];
						}
						else {
							*sum += (val[ir] * x0[(int)col[ir]]); //sum non-zero values in row vec
						}
					}
				}
			}

			x[ix] = (b[ix] - *sum) / *aii;
			//x_diff[ix] = x[ix] - x0[ix];
			//printf("x at %d is %f. sum is %0.1f, aii is %0.1f\n", ix, x[ix], *sum, *aii);

		} //end for (vals of x)

		//compute the norm of vec x
		cublasDnrm2(handle, x_size, x, 0, x_norm);
		//compute the norm of x_diff
		cublasDnrm2(handle, x_size, x_diff, 0, x_diff); //put norm in first location of x_diff

		//check if error is within tolerance
		double error = *x_diff / *x_norm;
		if (error <= tol) {
			printf("convergence...");
			break; //convergence reached
		} //end if

	} //end for (iterations)

	//wait for all processes to complete
	cudaDeviceSynchronize();

	cublasDestroy(handle);

	//free device memory
	cudaFree(d_val);
	cudaFree(d_row);
	cudaFree(d_col);
	cudaFree(d_x0);
	cudaFree(d_sum);
	cudaFree(d_aii);

	return x;
}

double* serial_sor(double *val, float *row, float *col, double *b, const unsigned int n, const unsigned int x_size, const unsigned int max_iter) {
	//allocate output vectors x and x0
	double *x = new double[x_size];
	double *x0 = new double[x_size];
	double sum = 0.0, aii = 0.0;


	for (int i = 0; i < 50; i++) {
		x0 = x;
		printf("x at iter %d is %f\n", i, *x);

		for (int ix = 0; ix < x_size; ix++) {
			sum = 0; //reset sum
			
			for (int ir = 0; ir < n; ir++) {
				if (row[ir] == ix) {
					if (col[ir] == ix) {
						aii = val[ir];
					}
					else {
						sum += (val[ir] * x0[(int)col[ir]]); //sum non-zero values in row vec
					}
				}
			}
			x[ix] = (b[ix] - sum) / aii; //update position of x
		}
	}
	return x;
}