#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>
#include <vector>

class GeneSequence {
private:
	char* name;
	char* sequence;
public:
	GeneSequence(char* n, char* s) :name(n), sequence(s) {}
	~GeneSequence(void);
	void setName(char* name);
	void setSequence(char* sequence);
	char* getName();
	char* getSequence();
};

void parseFASTA(FILE *infile, std::vector<GeneSequence>);

void dna(FILE *infile);
void rna(FILE *infile);
void revc(FILE *infile);
void fib(FILE* infile);

cudaError_t ntcount(const char* in, int* out, unsigned int size);
cudaError_t rnatranscribe(const char* in, char* out, unsigned int size);
cudaError_t dnacomplement(const char* in, char* out, unsigned int size);
