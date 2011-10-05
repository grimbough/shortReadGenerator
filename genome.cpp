/*
 * genome.cpp
 *
 *  Created on: 5 Oct 2011
 *      Author: smith68
 */

#include "genome.h"
#include "shortReadGenerator.h"

Genome::Genome() {
	numChrom = 0;
}

void Genome::addChromosome (string name, string sequence) {
	numChrom++;
	chromVec.push_back(Chromosome(name, sequence));
}

void Genome::printChromosomes() {
	for(int i = 0; i < numChrom; i++) {
		chromVec[i].printName();
	}
}

double* Genome::chromosomeWeights() {

        int i;
        vector<double> chromLengths;
        double *weights = (double *) calloc(sizeof(double), numChrom);

        /* get all the chromosome lengths */
        for(i = 0; i < numChrom; i++) {
            chromLengths.push_back( (double) chromVec[i].length());
        }

        double max = *max_element(chromLengths.begin(), chromLengths.end());

        for(i = 0; i < numChrom; i++) {
            chromLengths[i] /= max ;
        }

        double sum_of_elems = 0;
        //for(std::vector<double>::iterator j=chromLengths.begin(); j!=chromLengths.end(); ++j)
        //    sum_of_elems += *j;
        for(i = 0; i < numChrom; i++)
            sum_of_elems += chromLengths[i];

        for(i = 0; i < numChrom; i++) {
            weights[i] = chromLengths[i] / sum_of_elems;

        }

        return(weights);
}


Reads Genome::generateReads(int totalReads, int readLength1, int readLength2, int insertSizeMean, int insertSizeSD, int seed) {

	int i, chrom, pos;
	string end1, end2, start1, start2;
	size_t found1, found2;
	Reads reads = Reads(readLength1, readLength2, insertSizeMean, insertSizeSD);
        vector<int> chromLengths;

        std::cout << "Generating reads" << std::endl;

        double *weights = chromosomeWeights();
        //double weights[1] = {1.0};

        /* get all the chromosome lengths */
        for(i = 0; i < numChrom; i++) {
            chromLengths.push_back(chromVec[i].length());
        }


        for(chrom = 0; chrom < numChrom; chrom++) {

            int numReads = (int) floor( ( weights[chrom] * totalReads ) + 0.5);

            for(i = 0; i < numReads; i++) {

                int insertSize = SampleNormal((double) insertSizeMean, (double) insertSizeSD, seed);
                unsigned int internalSize = insertSize - readLength1 - readLength2;

                    /// choose whether the first read is from the +ve or -ve strand
                    int firstReadStrand = (int) SampleBernoulli(0.5, seed);

                    if(firstReadStrand) { // if 1 first read is on positive strand
                        pos = SampleUniform(0, chromVec[chrom].length() - (readLength1 + readLength2 + internalSize + 1), seed);
                        end1 = chromVec[chrom].returnSequence().substr(pos, readLength1);
                        end2 = chromVec[chrom].returnSequence().substr(pos + readLength1 + internalSize, readLength2);
                        end2 = reverseComplement(end2);
                    }
                    else { // if 0 first read -ve strand
                        pos = SampleUniform(readLength2 + internalSize, chromVec[chrom].length() - (readLength1 + 1), seed);
                        end1 = chromVec[chrom].returnSequence().substr(pos, readLength1);
                        end1 = reverseComplement(end1);
                        end2 = chromVec[chrom].returnSequence().substr(pos - internalSize - readLength2, readLength2);
                    }

                    found1 = end1.find("N");
                    found2 = end2.find("N");

                    while(found1 != string::npos || found2 != string::npos) {

                        if(firstReadStrand) { // if 1 first read is on positive strand
                            pos = SampleUniform(0, chromVec[chrom].length() - (readLength1 + readLength2 + internalSize + 1), seed);
                            end1 = chromVec[chrom].returnSequence().substr(pos, readLength1);
                            end2 = chromVec[chrom].returnSequence().substr(pos + readLength1 + internalSize, readLength2);
                            end2 = reverseComplement(end2);
                        }
                        else { // if 0 first read -ve strand
                            pos = SampleUniform(readLength2 + internalSize, chromVec[chrom].length() - (readLength1 + 1), seed);
                            end1 = chromVec[chrom].returnSequence().substr(pos, readLength1);
                            end1 = reverseComplement(end1);
                            end2 = chromVec[chrom].returnSequence().substr(pos - internalSize - readLength2, readLength2);
                        }

                        found1 = end1.find("N");
                        found2 = end2.find("N");
                    }

                    start1 = boost::lexical_cast< string >( pos + 1 );
                    if(firstReadStrand)
                        start2 = boost::lexical_cast< string >( pos + readLength1 + internalSize + 1 );
                    else
                        start2 = boost::lexical_cast< string >( pos - readLength2 - internalSize + 1 );

                    reads.addReads(end1, end2, chromVec[chrom].returnName(), chromVec[chrom].returnName(), start1, start2, firstReadStrand);

            }

        }

        free(weights);

	return(reads);
}

