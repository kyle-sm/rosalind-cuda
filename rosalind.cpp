#include "rosalind.h"

int main(int argc, char** argv)
{
	if (argc < 3) {
		fprintf(stderr, "Proper usage: rosalind <problem number> <input file>\n");
	}
	
	int prob = atoi(argv[1]);
	FILE* infile;
	fopen_s(&infile, argv[2], "r");

	switch (prob) {
	case(1):
		dna(infile);
	case(2):
		rna(infile);
	case(3):
		revc(infile);
	case(4):
		fib(infile);
	case(5):
		gc(argv[2]);
	case(6):
		hamm(argv[2]);
	}


	return 0;
}

void dna(FILE *infile) {
	char* in = (char*)malloc(sizeof(char) * 1000);
	if (in == 0) {
		printf("Failed allocating memory for input.\n");
		exit(1);
	}
	fgets(in, 1000, infile);
	int length = strlen(in);
	int out1[4] = { 0 };
	cudaError_t status = ntcount(in, out1, length);
	free(in);
	if (status != cudaSuccess) {
		exit(1);
	}
	printf("%d %d %d %d\n", out1[0], out1[1], out1[2], out1[3]);
	status = cudaDeviceReset();
	if (status != cudaSuccess) {
		fprintf(stderr, "Couldn't reset device: %s", cudaGetErrorString(status));
		exit(1);
	}
}

void rna(FILE *infile) {
	char* in = (char*)malloc(sizeof(char) * 1000);
	if (in == 0) {
		printf("Failed allocating memory for input.\n");
		exit(1);
	}
	fgets(in, 1000, infile);
	int length = strlen(in);
	char* out = (char*)malloc(length);
	cudaError_t status = rnatranscribe(in, out, length);
	free(in);
	printf("%s", out);
	free(out);
	status = cudaDeviceReset();
	if (status != cudaSuccess) {
		fprintf(stderr, "Couldn't reset device: %s", cudaGetErrorString(status));
		exit(1);
	}
}

void revc(FILE *infile) {
	char* in = (char*)malloc(sizeof(char) * 1000);
	if (in == 0) {
		printf("Failed allocating memory for input.\n");
		exit(1);
	}
	fgets(in, 1000, infile);
	int length = strlen(in);
	char* out = (char*)malloc(length);
	cudaError_t status = dnacomplement(in, out, length);
	printf("%s", out);
	status = cudaDeviceReset();
	if (status != cudaSuccess) {
		fprintf(stderr, "Couldn't reset device: %s", cudaGetErrorString(status));
		exit(1);
	}
}

/* This one can't be parallelized, but it's here for completion's sake */
void fib(FILE *infile) {
	unsigned int k = 0;
	unsigned int n = 0;
	if (fscanf(infile, "%d %d", &n, &k) < 2) {
		fprintf(stderr, "Couldn't read input file.\n");
		exit(1);
	}

	if (n > 50 || k > 5) {
		fprintf(stderr, "n must be less than or equal to 50 and k must be less than or equal to 5.");
	}

	fprintf(stderr, "%d %d\n", n, k);
	/* We can cheat and avoid any calculations for values of n less than or equal to 3. I don't know if this is actually more efficient or
	more annoying. */
	if (n < 3) {
		printf("1\n");
		exit(0);
	}
	if (n == 3) {
		printf("%d\n", k + 1);
		exit(0);
	}

	/* We need only the most recent 2 values in the sequence. Since the first two values will always be one, we set them ahead of time.
	Since the third value is always k + 1, we set that too. */
	unsigned long long rabbits[2] = { 1, k + 1 };

	for (int i = 3; i < n; i++) {
		unsigned long long temp = (rabbits[0] * k) + rabbits[1];
		rabbits[0] = rabbits[1];
		rabbits[1] = temp;
	}

	printf("%llu\n", rabbits[1]);
}

void gc(char* filename) {
	int nucleotideCount[4];
	float gcmax = 0.0f;
	char name[100];
	std::vector<GeneSequence> genes = parseFASTA(filename);
	for (auto& gene : genes) {
		cudaError_t status = ntcount(gene.getSequence(), nucleotideCount, strlen(gene.getSequence()));
		/* Nucleotide counts are in the following indices: {A, G, C, T} */
		float gc = (float)(nucleotideCount[1] + nucleotideCount[2]) /
			(float)(nucleotideCount[0] + nucleotideCount[1] + nucleotideCount[2] + nucleotideCount[3]);
		if (gc > gcmax) {
			gcmax = gc;
			strcpy(name, gene.getName());
		}
	}
	printf("%s\n%.6f\n", name, gcmax * 100);
}

void hamm(char* filename) {
	std::ifstream infile;
	char first[1000];
	char second[1000];
	int hammingDistance;
	infile.open(filename);
	infile.getline(first, 1000);
	infile.getline(second, 1000);
	cudaError_t status = hammcuda(first, second, strlen(first), &hammingDistance);
	printf("%d\n", hammingDistance);
}

std::vector<GeneSequence> parseFASTA(char* filename) {
	std::ifstream infile;
	std::vector<GeneSequence> genes;
	char buf[100];
	char name[100];
	char seq[1000] = { '\0' };

	infile.open(filename);
	while (infile.getline(buf, 100)) {
		if (buf[0] == '>')
			strcpy(name, &buf[1]); /* chop off the leading '>' character */
		else
			strcat(seq, buf);
		/* If the next character is a >, it indicates the start of a new gene. If it's EOF, there are no more
		genes to process. Either way, we're done reading the current gene and can send it to the vector. */
		if (infile.peek() == '>' || infile.peek() == EOF) {
			genes.push_back(GeneSequence(name, seq));
			seq[0] = '\0';
		}
	}
	infile.close();

	return genes;
}