//---------------------------------------------------------------------------
#ifndef TDataVectorH
#define TDataVectorH
#include "global.h"
#include "TParameters.h"
#include "TExecuteProgram.h"
#include "TData.h"
#include "TLinearComb.h"
#include "TInputFileVector.h"
#include <vector>
#include <strstream>
#include <sys/stat.h>
#include <dirent.h>
//---------------------------------------------------------------------------
class TDataVector{
   private:
		double calculateDistanceLinearComb(TLinearComb* pcaObject);

   public:
		vector<TData> vecDataObjects;
		vector<TData>::iterator curDataObject, endDataObject;
		int numDataObjects;
		ofstream* logFile;
		int numObsData;
		bool stdLinearCombForDist;
		my_string arpDirectory;
		my_string sumStatsTempFileName;
		double** pointersToObsData;
		double** pointersToSimData;
		double** pointersToVariances;
		my_string* obsDataNamesWithPrefix;
		TInputFileVector* inputFiles;
		bool useSumstatProgram;
		bool launchScriptAfterSS;
		bool launchScriptBeforeSS;

		TDataVector(TParameters* GotParameters, TInputFileVector* gotInputFiles, ofstream* gotlogfile);

		double calculateDistance();
		double calculateDistance(TLinearComb* pcaObject);
		double calculateDistance(double* pointerToData);
		double calculateDistance(double* pointerToData, TLinearComb* pcaObject);
		my_string locateFileRecursively(my_string file, my_string dir);
		void initializeDataObjects();
		void summarizeDataObjects();
		void calculateSumStats(int simnum);
		my_string getObsFileNames();
		void saveOldValues();
		void resetOldValues();
		void writeHeader(ofstream& ofs);
		void writeHeaderOnlyOneDataObject(int thisObs, ofstream& ofs);
		void writeData(ofstream& ofs);
		void writeDataOnlyOneDataObject(int thisObsofstream, ofstream& ofs);
		bool matchColumns(vector<my_string> & colNames, int numColNames, int** pointerToArrayToFill);
		vector<my_string> getNamesVector();
		vector<double> getObsDataVector();
		void fillObsDataArray(double* & obsData);
};

#endif
