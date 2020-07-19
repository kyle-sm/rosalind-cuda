#include "rosalind.h"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
// Problem 6: Counting Point Mutations

__global__ void hammkernel(const char* first, const char* second,  int* hamming) {
	if (first[threadIdx.x] == second[threadIdx.x])
		atomicAdd(hamming, 1);
}

cudaError_t hammcuda(const char* first, const char* second, int size, int* hamming) {
	char* devfirst = 0;
	char* devsecond = 0;
	int* devhamming = 0;
	cudaError_t cudaStatus;

	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Couldn't find CUDA-compatible device.\n");
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&devfirst, size * sizeof(char));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Couldn't allocate memory for input on GPU.\n");
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&devsecond, size * sizeof(char));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Couldn't allocate memory for input on GPU.\n");
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&devhamming, sizeof(int));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Couldn't allocate memory for output on GPU.\n");
		goto Error;
	}

	cudaStatus = cudaMemcpy(devfirst, first, size * sizeof(char), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Couldn't copy input to GPU: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	cudaStatus = cudaMemcpy(devsecond, second, size * sizeof(char), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Couldn't copy input to GPU: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	hammkernel <<<1, size>>> (devfirst, devsecond, devhamming);

	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Couldn't execute kernel to calculate hamming distance: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Couldn't synchronize threads for hammkernel: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	cudaStatus = cudaMemcpy(hamming, devhamming, sizeof(int), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Couldn't copy output to host device.\n");
		goto Error;
	}

Error:
	cudaFree(devfirst);
	cudaFree(devsecond);
	cudaFree(devhamming);

	return cudaStatus;
}

