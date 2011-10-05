/*
 * reads.cpp
 *
 *  Created on: 5 Oct 2011
 *      Author: smith68
 */

#include "reads.h"
#include "shortReadGenerator.h"

void Reads::addReads(string read1, string read2, string chr1, string chr2, string p1, string p2, int firstReadStrand) {

	end1.push_back(read1);
	end2.push_back(read2);
	chrom1.push_back(chr1);
	chrom2.push_back(chr2);
	pos1.push_back(p1);
	pos2.push_back(p2);
    strand1.push_back((bool) firstReadStrand);
}

void Reads::generateQualities(std::string tmpQuals1, std::string tmpQuals2, int seed) {

    int nlines, pos;
    unsigned i;
    int qualsLength1, qualsLength2;
    std::string buff;

    std::cout << "Generating Qualities" << std::endl;

    std::ifstream qualsf1(tmpQuals1.c_str());
    std::ifstream qualsf2(tmpQuals2.c_str());
    if (!qualsf1.is_open() | !qualsf2.is_open())
    {
        std::cerr << "error reading temporary qualities file" << std::endl;
        exit (EXIT_FAILURE); // error
    }
    /* find the number of lines */
    ///first read one line from each file and find their lengths
    getline(qualsf1, buff);
    qualsLength1 = buff.length();
    getline(qualsf2, buff);
    qualsLength2 = buff.length();

    /// now go to the end of the file
    qualsf1.seekg(0, std::ios::end);
    qualsf2.seekg(0, std::ios::end);

    if ( (int) (qualsf1.tellg() / (qualsLength1 + 1)) != (int) (qualsf2.tellg() / (qualsLength2 + 1)) ) {
        cerr << "Number of lines in qualities files do not match" << endl;
        exit(EXIT_FAILURE);
    }

    nlines = (int) (qualsf2.tellg() / (qualsLength2 + 1));

    for(i = 0; i < end1.size(); i++) {
        int line = SampleUniform(0, nlines - 1, seed);
        pos = (qualsLength1 + 1) * line;

        qualsf1.seekg(pos, std::ios::beg);
        getline(qualsf1, buff);
        if(buff.length() < readLength1) {
            buff = interpolateQualities(buff, readLength1);
        }
        quals1.push_back(buff);

        pos = (qualsLength2 + 1) * line;
        qualsf2.seekg(pos, std::ios::beg);
        getline(qualsf2, buff);
        if(buff.length() < readLength2) {
            buff = interpolateQualities(buff, readLength2);
        }
        quals2.push_back(buff);
    }
    //cout << "Made it" << endl;
    qualsf1.close();
    qualsf2.close();

}


void Reads::applyQualities(int seed, int qualityAdjust, double snpProb) {

    std::string bases = "ACGT";
    unsigned i, j;
    int q, baseChoice;
    double p;

    cout << "Applying qualities" << endl;

    for(i = 0; i < end1.size(); i++) {
        for(j = 0; j < end1[i].length(); j++) {

            // get the quality value
            q = (int) quals1[i][j] - qualityAdjust;
            // convert into a probability of a miscall
            p = pow(10, (double) -q / 10);
            if(p != 1)
                p += snpProb;

            if(p > 1) {
                cout << "End1: " << quals1[i] << endl;
                cout << i << "\t" << j << "\t" << p << "\t" << quals1[i][j] << endl;
                exit(EXIT_FAILURE);
            }

            if(SampleBernoulli(p, seed)) {
                baseChoice = SampleUniform(0, 3, seed);
                end1[i][j] = bases[baseChoice];
            }
        }
    }

    for(i = 0; i < end2.size(); i++) {
        for(j = 0; j < end2[i].length(); j++) {
            // repeat for the second read
            q = (int) quals2[i][j] - qualityAdjust;
            p = pow(10, (double) -q / 10);
            if(p != 1)
                p += snpProb;

            if(p > 1) {
                 cout << "End: " << quals2[i] << endl;
                 cout << i << "\t" << j << "\t" << p << "\t" << quals2[i][j] << endl;
                 exit(EXIT_FAILURE);
             }

            if(SampleBernoulli(p, seed)) {
                baseChoice = SampleUniform(0, 3, seed);
                end2[i][j] = bases[baseChoice];
            }

        }
    }

}

void Reads::writeFASTA() {
	unsigned i;
	for(i = 0; i < end1.size(); i++) {
		cout << chrom1[i] << ":" << pos1[i] << "/1\n" << end1[i] << "\n";
	}

	for(i = 0; i < end1.size(); i++) {
		cout << chrom2[i] << ":" << pos1[i] << "/2\n" << end2[i] << "\n";
	}
}

void Reads::writeFASTA(string fileName) {

        unsigned i;
        ofstream file1;
        ofstream file2;

        file1.open("test1.fa");
        for(i = 0; i < end1.size(); i++) {
                file1 << chrom1[i] << ":" << pos1[i] << "/1\n" << end1[i] << "\n";
        }
        file2.close();

        file2.open("test2.fa");
        for(i = 0; i < end1.size(); i++) {
                file2 << chrom2[i] << ":" << pos2[i] << "/2\n" << end2[i] << "\n";
        }
        file2.close();
}

void Reads::writeFASTQ(string fileName) {

        unsigned i;
        ofstream file1;
        ofstream file2;

        std::string name1 = fileName + "_1.fq";
        std::string name2 = fileName + "_2.fq";

        file1.open(name1.c_str());
        for(i = 0; i < end1.size(); i++) {
                file1 << "@" << chrom1[i].substr(1) << ":" << pos1[i] << ":" << pos2[i] << ":";

                if(strand1[i])
                    file1 << "+";
                else
                    file1 << "-";

                file1 << "/1" << std::endl;
                file1 << end1[i] << std::endl;
                file1 << "+" << std::endl;
                file1 << quals1[i] << std::endl;
        }
        file2.close();

        file2.open(name2.c_str());
        for(i = 0; i < end1.size(); i++) {
                file2 << "@" << chrom2[i].substr(1) << ":" << pos1[i] << ":" << pos2[i] << ":";

                if(strand1[i])
                    file2 << "+";
                else
                    file2 << "-";

                file2 << "/2" << std::endl;
                file2 << end2[i] << std::endl;
                file2 << "+" << std::endl;
                file2 << quals2[i] << std::endl;
        }
        file2.close();
}



