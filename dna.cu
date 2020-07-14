#include "rosalind.h"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
// Problem 1: Counting DNA Nucleotides

__global__ void ntcountkernel(const char* in, int* out)
{
	switch (in[threadIdx.x]) {
	case 'A':
		atomicAdd(&out[0], 1);
		break;
	case 'C':
		atomicAdd(&out[1], 1);
		break;
	case 'G':
		atomicAdd(&out[2], 1);
		break;
	case 'T':
		atomicAdd(&out[3], 1);
		break;
	default:
		break;
	}
}

cudaError_t ntcount(const char* in, int* out, unsigned int size)
{
	char* devin = 0;
	int* devout = 0;
	cudaError_t cudaStatus;

	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Couldn't find CUDA-compatible device.\n");
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&devin, size * sizeof(char));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Couldn't allocate memory for input on GPU.\n");
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&devout, 4 * sizeof(int));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Couldn't allocate memory for output on GPU.\n");
		goto Error;
	}

	cudaStatus = cudaMemcpy(devin, in, size * sizeof(char), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Couldn't copy input to GPU: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	ntcountkernel <<<1, size>>> (devin, devout);

	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Couldn't execute ntcountkernel: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Couldn't synchronize threads for ntcountkernel: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	cudaStatus = cudaMemcpy(out, devout, 4 * sizeof(int), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Couldn't copy output to host device.\n");
		goto Error;
	}


Error:
	cudaFree(devin);
	cudaFree(devout);

	return cudaStatus;
}