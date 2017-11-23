//---------------------------------------------------------------------------

#ifndef TInputFileH
#define TInputFileH
#include "global.h"
#include "TPriorVector.h"
#include "TExecuteProgram.h"
//---------------------------------------------------------------------------
class TInputFile{
   public:
	  vector<my_string> simulationProgramInputFilenames;
	  vector<my_string> newSimulationProgramInputFilenames;

	  bool checkParametersPassedToSimulationProgram;
	  bool checkParametersPassedToScriptBeforeSimulations;
	  bool checkParametersPassedToScriptAfterSimulations;
	  TPriorVector* priors;
	  TExecuteProgram* simulationProgram;
	  TExecuteProgram* scriptBeforeSimulations;
	  TExecuteProgram* scriptAfterSimulations;
	  my_string* valuesToPassToSimulationProgramm;
	  my_string* valuesToPassToScriptBeforeSimulations;
	  my_string* valuesToPassToScriptAfterSimulations;
	  int numValuesToPassToSimulationProgramm;
	  int numValuesToPassToScriptBeforeSimulations;
	  int numValuesToPassToScriptAfterSimulations;

	  TInputFile(my_string siminputname, TPriorVector* gotPriors);
	  ~TInputFile(){};

	  void writeNamesToLogfile(ofstream* logfile);
	  void writeNewNamesToLogfile(ofstream* logfile);
	  void createNewInputFile();
	  void initializeSimulationProgramm(my_string simulationprogram, my_string gotParameter);
	  int performSimulation(int simnum=-1);
	  void initializeScriptBeforeSimulation(my_string gotScript, my_string gotParameter);
	  int launchScriptBeforeSimulations(int simnum=-1);
	  void initializeScriptAfterSimulation(my_string gotScript, my_string gotParameter);
	  int launchScriptAfterSimulations(int simnum=-1);


};

#endif
