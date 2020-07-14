#include "rosalind.h"

GeneSequence::~GeneSequence(void) {
	free(name);
	free(sequence);
}

char* GeneSequence::getName() {
	return this->name;
}

char* GeneSequence::getSequence() {
	return this->sequence;
}

void GeneSequence::setName(char* name) {
	this->name = name;
}

void GeneSequence::setSequence(char* sequence) {
	this->sequence = sequence;
}