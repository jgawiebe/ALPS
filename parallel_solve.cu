#include <stdlib.h>  
#include <stdio.h>
#include <cstdlib>  
#include <cuda_runtime.h>

#define BLOCKS 512
#define THREADS 1024

__global__ void kernel_solve(double *val, double *row_ix, double *col_ix, double* b, double* x, unsigned int n) {
	__shared__ int g_val, g_row, g_col, g_x;
	const unsigned int index = threadIdx.x + blockIdx.x * blockDim.x;
	if (index < n) {
		g_val = val[index];
		g_row = row_ix[index];
		g_col = col_ix[index];
		//result
		x[index] = g_row * g_col;

		printf("Thread %d - value: %d at (%d, %d)\n", index, g_val, g_row, g_col);
	}


}
double* parallel_sor(double *val_ix, double *row_ix, double *col_ix, double *b, const int n) {
	//allocate output vector x
	double *x = new double[n];

	//declare memory size for vectors
	const unsigned int memsize = sizeof(double) * n;

	//declare pointers to device memory
	double *d_val, *d_row, *d_col, *d_b, *d_x;

	//allocate device memory for each vector
	cudaMalloc((void **)&d_val, memsize);
	cudaMalloc((void **)&d_row, memsize);
	cudaMalloc((void **)&d_col, memsize);
	cudaMalloc((void **)&d_b, memsize);
	cudaMalloc((void **)&d_x, memsize);

	//copy memory and link pointers
	cudaMemcpy(d_val, val_ix, memsize, cudaMemcpyHostToDevice);
	cudaMemcpy(d_row, row_ix, memsize, cudaMemcpyHostToDevice);
	cudaMemcpy(d_col, col_ix, memsize, cudaMemcpyHostToDevice);
	cudaMemcpy(d_b, b, memsize, cudaMemcpyHostToDevice);
	cudaMemcpy(d_x, x, memsize, cudaMemcpyHostToDevice);

	dim3 threadsPerBlock(16, 16);
	dim3 numBlocks(n / threadsPerBlock.x, n / threadsPerBlock.y);
	
	//execute solve in parallel
	kernel_solve <<<16, 16>>> (d_val, d_row, d_col, d_b, d_x, n);

	//return result
	cudaMemcpy(x, d_x, memsize, cudaMemcpyDeviceToHost);

	//wait for all processes to complete
	cudaDeviceSynchronize();

	//free device memory
	cudaFree(d_val);
	cudaFree(d_row);
	cudaFree(d_col);
	cudaFree(d_b);
	cudaFree(d_x);

	return x;
}

double* serial_sor(double *val, double *row, double *col, double *b, const unsigned int n, const unsigned int x_size, const unsigned int max_iter) {
	//allocate output vectors x and x0
	double *x = new double[x_size];
	double *x0 = new double[x_size];
	double sum = 0.0, aii = 0.0;


	for (int i = 0; i < max_iter; i++) {
		x0 = x;

		for (int ix = 0; ix < x_size; ix++) {
			sum = 0; //reset sum
			
			for (int ir = 0; ir < n; ir++) {
				if (row[ir] == ix) {
					if (col[ir] == ix) {
						aii = val[ir];
					}
					else {
						sum += ( val[ir] * x0[(int) col[ir]] ); //sum non-zero values in row vec
					}
				}
			}
			x[ix] = (b[ix] - sum) / aii; //update position of x
		}
	}
	return x;
}