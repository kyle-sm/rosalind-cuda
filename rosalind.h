#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>
#include <fstream>
#include <vector>

class GeneSequence {
private:
	char* name;
	char* sequence;
public:
	GeneSequence(char* n, char* s);
	~GeneSequence(void);
	void setName(char* name);
	void setSequence(char* sequence);
	char* getName();
	char* getSequence();
};

std::vector<GeneSequence> parseFASTA(char* filename);

void dna(FILE *infile);
void rna(FILE *infile);
void revc(FILE *infile);
void fib(FILE* infile);
void gc(char* filename);
void hamm(char* filename);

cudaError_t ntcount(const char* in, int* out, unsigned int size);
cudaError_t rnatranscribe(const char* in, char* out, unsigned int size);
cudaError_t dnacomplement(const char* in, char* out, unsigned int size);
cudaError_t hammcuda(const char* first, const char* second, int size, int* hamming);
