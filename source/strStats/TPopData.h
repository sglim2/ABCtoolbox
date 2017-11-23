/*
 * TPopData.h
 *
 *  Created on: Feb 5, 2009
 *      Author: wegmannd
 */

#ifndef TPOPDATA_H_
#define TPOPDATA_H_

#include "my_cstring.h";
#include "TException.h"
#include <strstream>
#include <fstream>
#include <iostream.h>
#include <vector>
#include "math.h"

class TLocus{
private:
	int* alleles;
	int readAlleles;
	int size;
	double mean;
	double var;
	int numNotMissing;
	bool meanCalculated, varCalculated;

public:
	TLocus(int gotSize);
	void add(int allele);
	double getVariance();
	double getMean();
};

class TPopData{
private:
	my_string name;
		int sampleSize;
		bool isGenotypic;
		bool nextPhase;
		int readLines;
		vector<TLocus*> loci;
		int numLoci;

public:
	TPopData(){};
	TPopData(my_string gotName, int gotSize, bool gotIsGenotypic);
	void readArpfileLine(my_string line);
	double getMean();
	double getVariance();
	int getNumLoci();
	double getMeanOneLocus(int locus);
	double getVarianceOneLocus(int locus);

};

class TPopDataVector{
private:
	my_string missingDataString;
	bool isGenotypic;
	vector<TPopData*> myPopData;

public:
	TPopDataVector(){};
	TPopDataVector(my_string filename);
	void writeStats(ofstream* output);
	void writeStatsForEachLocusIndividually(ofstream* output);
	void writeHeader(ofstream* output);
};




#endif /* TPOPDATA_H_ */
