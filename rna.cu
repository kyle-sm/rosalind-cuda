#include "rosalind.h"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
// Problem 2: Transcribing DNA into RNA

__global__ void transcribekernel(const char* in, char* out)
{
	if (in[threadIdx.x] == 'T')
		out[threadIdx.x] = 'U';
	else
		out[threadIdx.x] = in[threadIdx.x];
}

cudaError_t rnatranscribe(const char* in, char* out, unsigned int size)
{
	char* devin = 0;
	char* devout = 0;
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

	cudaStatus = cudaMalloc((void**)&devout, size * sizeof(char));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Couldn't allocate memory for output on GPU.\n");
		goto Error;
	}

	cudaStatus = cudaMemcpy(devin, in, size * sizeof(char), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Couldn't copy input to GPU: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	transcribekernel <<<1, size>>> (devin, devout);

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

	cudaStatus = cudaMemcpy(out, devout, size * sizeof(char), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Couldn't copy output to host device.\n");
		goto Error;
	}

Error:
	cudaFree(devin);
	cudaFree(devout);

	return cudaStatus;
}