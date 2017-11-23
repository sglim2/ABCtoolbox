#ifndef TDataH
#define TDataH

#include "global.h"
#include "TParameters.h"
#include "TExecuteProgram.h"
#include "TPriorVector.h"
#include <strstream>
#include <vector>

//---------------------------------------------------------------------------
class TData{
   public:
	  int numObsData;
	  double* obsData;
	  double* oldData;
	  double* obsDataVariance; //assed during calibration step
	  my_string* obsDataName;
	  my_string obsFileName;
	  double** simDataPointers;
	  double* simDataInput;
	  my_string* simDataName;
	  int numSimData;
	  ofstream* logFile;
	  my_string simDataFileName;
	  bool doPCA;
	  TExecuteProgram* sumStatProgram;
	  bool checkParametersPassedToSumStatProgram;
	  int numValuesToPassToSumStatProgram;
	  my_string* valuesToPassToSumStatProgram;

	  TExecuteProgram* scriptBeforeSSCalc;
	  bool checkParametersPassedToScriptBeforeSSCalc;
	  int numValuesToPassToScriptBeforeSSCalc;
	  my_string* valuesToPassToScriptBeforeSSCalc;

	  TExecuteProgram* scriptAfterSSCalc;
	  bool checkParametersPassedToScriptAfterSSCalc;
	  int numValuesToPassToScriptAfterSSCalc;
	  my_string* valuesToPassToScriptAfterSSCalc;

	  TData();
	  TData(my_string obsname, my_string simDataName, ofstream* gotlogfile);
	  TData(my_string obsname, ofstream* gotlogfile);
	  ~TData(){};

	  void readObsfile();
	  void initializeScriptBeforeSSCalc(my_string gotScript, my_string gotParameter, my_string sumStatsTempFileName);
	  void initializeScriptAfterSSCalc(my_string gotScript, my_string gotParameter, my_string sumStatsTempFileName);
	  void initializeSumstatProgramm(my_string gotSumstatprogram, my_string gotParameter, my_string sumStatsTempFileName);
	  int launchScriptBeforeSSCalc(int simnum);
	  int launchScriptAfterSSCalc(int simnum);
	  void calculateSumStats(int simnum);
	  void readInitialSimData(my_string simFilename);
	  void readSimData(my_string simFilename);
	  void readSimDataLine(strstream & line);
	  void writeData(ofstream& ofs);
	  void writeHeader(ofstream& ofs);
	  void writeObsData(ofstream& ofs);
	  int getObsDataIndexFromName(my_string name);
	  void resetOldValues();
	  void saveOldValues();

};



#endif
