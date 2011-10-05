/*
 * genome.cpp
 *
 *  Created on: 5 Oct 2011
 *      Author: smith68
 */

#include "shortReadGenerator.h"
#include "chromosome.h"
#include "reads.h"

class Genome {
	int numChrom;
	vector<Chromosome> chromVec;
public:
	Genome();
	void addChromosome(string, string);
	void printChromosomes();
	Reads generateReads(int, int, int, int, int, int);
    double* chromosomeWeights();
};
