/*
 * shortReadGenerator.h
 *
 *  Created on: 5 Oct 2011
 *      Author: smith68
 */

#ifndef SHORTREADGENERATOR_H_
#define SHORTREADGENERATOR_H_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <ctime>
#include <sys/time.h>
#include <stdio.h>
#include <unistd.h>
#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
#include <boost/program_options.hpp>

using namespace std;
using namespace boost;
namespace po = boost::program_options;

/* sampling.cpp */
int SampleNormal (double, double, int);
int SampleUniform (int, int, int);
double SampleBernoulli (double, int);

/* Utility functions */
string reverseComplement(string) ;
string genRandomString(const int);
string trim(string&, const string&);
int processLine(string&);
string reverseComplement(string);
string interpolateQualities(string, int);
bool FileExists(string);
string createTmpQualitiesFile(string, string, int, string);
long wallTime(timeval, timeval);


#endif /* SHORTREADGENERATOR_H_ */
