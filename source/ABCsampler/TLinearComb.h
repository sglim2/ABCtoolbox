//---------------------------------------------------------------------------

#ifndef TLinearCombH
#define TLinearCombH
//---------------------------------------------------------------------------
#include <strstream>
#include <vector>
#include "my_cstring.h"

class TLinearComb{
	public:
	   vector<my_string> statNamesVector;
	   vector<my_string>::iterator curStatName, endStatName;
	   double* mean;
	   double* sd;
	   double* obs;
	   vector<double>::iterator curDouble, endDouble;
	   vector<my_string> namesInInputFiles;
	   vector<my_string>::iterator curNameInInputFiles, endNameInInputFiles;
	   my_string LinearCombFileName;
	   double** pca;
	   int* statIsUsed; // if>0 it points to the row in the linear comb file
	   int numLinearComb;
	   int numParams;
	   int numStats;
	   double* linearCombVariances;
	   double* obsPCA;
	   double* simPCA;
	   double* oldSimPCA;


	   TLinearComb(){};
	   TLinearComb(my_string gotLinearCombFileName, vector<my_string> NamesInInputFiles);
	   int checkName(my_string name);
	   my_string getFirstMissingStat();
	   int getNumberOfFirstUsedStat();
	   void writePCA(ofstream* output, vector<double> tempInput);
	   void writeDistance(ofstream* output, vector<double> tempInput);
	   virtual void readLinearCombFile();
	   virtual void calcObsLineraComb(vector<double> obsData);
	   double getPCA(int num, double* stats);
	   double getPCA(int num, double** stats);
	   virtual void calcSimDataPCA(double** simData);
	   virtual void calcSimDataPCA(double* simData);
	   void writeHeader(ofstream& ofs);
	   void writeSimPCA(ofstream& ofs);
	   void saveOldValues();
	   void resetOldValues();
};
#endif
