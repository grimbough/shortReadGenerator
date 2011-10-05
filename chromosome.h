/*
 * chromosome.h
 *
 *  Created on: 5 Oct 2011
 *      Author: smith68
 */

#include "shortReadGenerator.h"

class Chromosome {
    string name;
	string sequence;
public:
	Chromosome(string n, string s) { name = n; sequence = s; };
	void printName()  { cout << name << "\n"; };
	void printSequence() { cout << sequence << "\n"; };
	int length() { return(sequence.length()); };
	string returnSequence() { return(sequence); };
	string returnName() { return(name); };
};
