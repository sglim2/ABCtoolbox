#ifndef TsimH
#define TsimH
//#include "global.h"
#include "TInputFileVector.h"
#include "TDataVector.h"
#include "TParameters.h"
#include "TExecuteProgram.h"
#include "TLinearComb.h"
#include "TLinearCombBoxCox.h"
#include "TOutput.h"
#include <vector>

class TSim{

	public:
	  my_string outName;
	  bool append;
	  bool writeSeparateOutputFiles;
	  bool addDistanceToOutputfile;
	  vector<TInputFile> vecDataObjects;
	  vector<TData>::iterator curDataObject, endDataObject;
	  TInputFileVector* inputFiles;
	  ofstream est;
	  my_string exeDir;
	  my_string arpFile;
	  my_string arpDirectory;
	  my_string sumStatsTempFileName;
	  ofstream* logFile;
	  TParameters* gotParameters;
	  TDataVector* myData;
	  bool doLinearComb, doBoxCox;
	  my_string linearCombFileName;
	  TLinearComb* myLinearComb;

	  virtual int runSimulations();
	  void createLinearCombObject();
	  double calcDistance();
      TSim(){};
      TSim(const TSim & sim);
		TSim(TParameters* gotParameters, my_string gotexedir, ofstream* gotLogFile);
      virtual ~TSim(){};
};
#endif
