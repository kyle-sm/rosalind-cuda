#include "rosalind.h"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
// Problem 2: Transcribing DNA into RNA

__global__ void revckernel(const char* in, char* out)
{
	int oindex = blockDim.x - threadIdx.x;
	switch (in[threadIdx.x]) {
	case 'A':
		out[oindex] = 'T';
		break;
	case 'C':
		out[oindex] = 'G';
		break;
	case 'G':
		out[oindex] = 'C';
		break;
	case 'T':
		out[oindex] = 'A';
		break;
	default:
		break;
	}
}

cudaError_t dnacomplement(const char* in, char* out, unsigned int size)
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

	revckernel <<<1, size>>> (devin, devout);

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


	/* We do a bunch of terrible string hacks here because the null character keeps getting copied to the beginning of the string.
	I should come up with a more elegant solution, in case this comes up in future problems. */
	cudaStatus = cudaMemcpy(out, devout + 1, (size - 1) * sizeof(char), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Couldn't copy output to host device.\n");
		goto Error;
	}
	switch (in[0]) {
	case 'A':
		out[size - 1] = 'T';
		break;
	case 'C':
		out[size - 1] = 'G';
		break;
	case 'G':
		out[size - 1] = 'C';
		break;
	case 'T':
		out[size - 1] = 'A';
		break;
	default:
		break;
	}
	out[size] = '\0';

Error:
	cudaFree(devin);
	cudaFree(devout);

	return cudaStatus;
}