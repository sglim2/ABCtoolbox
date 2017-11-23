//---------------------------------------------------------------------------

#ifndef TInputFileVectorH
#define TInputFileVectorH

#include <vector>
#include "TPriorVector.h"
#include "global.h"
#include "TInputFile.h"
#include "TParameters.h"
#include "TExecuteProgram.h"
//---------------------------------------------------------------------------
class TInputFileVector{
	public:
	  vector<TInputFile> vecInputFileObjects;
	  vector<TInputFile>::iterator curInputFileObject, endInputFileObject;
	  ofstream* logFile;
	  TPriorVector* priors;
	  my_string estFile;
	  bool launchScriptAfterSimulation;
	  bool launchScriptBeforeSimulation;

	  TInputFileVector(TParameters* gotParameters, ofstream* gotLogFile);
	  void writeNewInputFiles(int simnum=-1);
	  void createNewInputFiles(int simnum=-1);
	  void createNewInputFilesMcmc(int simnum=-1);
	  void createNewInputFilesMcmcUpdateOnePriorOnly(int simnum=-1);
	  bool creatNewInputFilesPMC(double* newParams, int simnum);
	  void getNewValueMcmc(const int& priorNumber);
	  void performSimulations(int simnum=-1);
	  vector<my_string> getVectorOfInputFileNames();
	  double getPriorDensity();
	  double getOldPriorDensity();
	  double getPriorDensity(double* values);
	  int getNumberOfSimplePriorFromName(my_string name);
	  void setPriorValue(int priorNumber, double value);
      void updateCombinedParameters();
};





#endif
