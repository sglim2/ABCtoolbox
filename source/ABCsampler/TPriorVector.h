//---------------------------------------------------------------------------

#ifndef TPriorVectorH
#define TPriorVectorH
#include "global.h"
#include "TRule.h"
#include "TException.h"
#include "my_cstring.h"
#include <strstream>
#include <vector>
#include <map>
#include "newmat.h"
#include "newmatap.h"
#include "newmatio.h"
//---------------------------------------------------------------------------

class TPriorVector{
   friend istream& operator >> (istream& is, TPriorVector& S);
	friend ostream& operator << (ostream& os, const TPriorVector& S);
	public:
		map<my_string,TSimplePrior*> mapSimplePrior;
		map<my_string,TCombinedPrior*> mapCombinedPrior;
		map<my_string,TSimplePrior*>::iterator curSimplePrior, endSimplePrior;
		map<my_string,TCombinedPrior*>::iterator curCombinedPrior, endCombinedPrior;
		TRuleVector rules;
		int numPrior;
		int numSimplePrior;
		int numCombinedPrior;
		TSimplePrior** simplePriors;
		TCombinedPrior** combinedPriors;
		TPrior* curPrior;
		long _Idum;
		ofstream* logFile;

		TPriorVector(my_string fileName, ofstream* gotLogFile);
		~TPriorVector(){
		   mapSimplePrior.clear();
		   mapCombinedPrior.clear();
		}

		void getNewValues();
		void resetOldValues();
		void saveOldValues();
		void readPriorsAndRules(my_string fileName);
		void updateCombinedParameters();
		void writeHeader(ofstream& ofs);
		void writeHeaderSimplePriors(ofstream& ofs);
		void writeHeaderCombinedPriors(ofstream& ofs);
		void writeParameters(ofstream& ofs);
		void writeParametersSimplePriors(ofstream& ofs);
		void writeParametersCombinedPriors(ofstream& ofs);
		void getNewValuesMcmcUpdateOnePriorOnly();
		void getNewValuesMcmc(const int& priorNumber);
		void getNewValuesMcmc();
		void getNewValuesMcmc(TSimplePrior* thisSimplePrior);
		bool getNewValuesPMC(double* newParams);
		double getPriorDensity();
		double getOldPriorDensity();
		double getPriorDensity(double* values);
		double getValue(const my_string& Name);
		bool isPrior(const my_string& name);
		double calcEquation(my_string equat);
		bool writeCurValueToFileFromName(const my_string& name, ofstream& file);
	    bool writeCurValueWithHyperprior(const my_string& name, ofstream& file);
		TPrior* getPriorFromName(const my_string& name);
		TSimplePrior* getSimplePriorFromName(const my_string& name);
		TCombinedPrior* getCombinedPriorFromName(const my_string& name);
		int getNumberOfSimplePriorFromName(const my_string& name);
		void setSimplePriorValue(const int& priorNumber, const double& value);
};





#endif
