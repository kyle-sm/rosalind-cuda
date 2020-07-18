#include "rosalind.h"

GeneSequence::GeneSequence(char* n, char* s) {
	this->name = strdup(n);
	this->sequence = strdup(s);
}

GeneSequence::~GeneSequence(void) {
}

char* GeneSequence::getName() {
	return this->name;
}

char* GeneSequence::getSequence() {
	return this->sequence;
}

void GeneSequence::setName(char* name) {
	strcpy(this->name, name);
}

void GeneSequence::setSequence(char* sequence) {
	strcpy(this->sequence, sequence);
}