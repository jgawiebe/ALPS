#include <cuda_runtime.h>
#include <iostream>

//#include "red_black_wrapper.hpp"

__global__ void test() {
	printf("THIS IS A TEST");
	
}

int red_black_sor() {

	test<<<1, 1>>>();
	cudaDeviceSynchronize();
	return 0;
}