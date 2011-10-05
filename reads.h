/*
 * reads.cpp
 *
 *  Created on: 5 Oct 2011
 *      Author: smith68
 */

#include "shortReadGenerator.h"

class Reads {
	int readLength1;
    int readLength2;
	int insertSizeMean;
	int insertSizeSD;
	vector<string> end1;
	vector<string> end2;
    vector<string> quals1;
    vector<string> quals2;
	vector<string> chrom1;
	vector<string> chrom2;
	vector<string> pos1;
	vector<string> pos2;
    vector<bool> strand1;
public:
	Reads(int l1, int l2, int insM, int insSD) { readLength1 = l1; readLength2 = l2; insertSizeMean = insM; insertSizeSD = insSD;};
	void addReads(string, string, string, string, string, string, int);
    void generateQualities(std::string qualsFile1, std::string qualsFile2, int seed);
    void generateQualities(std::string qualsFile1, std::string qualsFile2, int seed, int readLength1, int readLength2);
    void applyQualities(int seed, int qualityAdjust, double snpProb);
    void applyQualities_MP(int seed, int qualityAdjust, double snpProb);
	void writeFASTA();
	void writeFASTA(string fileName);
    void writeFASTQ(string fileName);
};
