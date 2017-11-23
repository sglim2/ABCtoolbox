//---------------------------------------------------------------------------

#ifndef TOutputH
#define TOutputH
#include "global.h"
#include "TDataVector.h"
#include "TInputFileVector.h"
#include "TLinearComb.h"
#include <vector>
//---------------------------------------------------------------------------

class TOutput{
	public:
		int sampling; //every ith simulation is recorded
		my_string filename;
		ofstream myFile;
		bool writeDistance;
		TDataVector* data;
		TInputFileVector* inputFiles;

		TOutput(){};
		TOutput(int gotSampling, my_string outputFilename, bool gotWriteDistance, TDataVector* Data, TInputFileVector* InputFiles);
		~TOutput(){myFile.close();};
		virtual void writeHeader();
        void writeHeader(TLinearComb* pcaObject);
		virtual void writeSimulations(int run, double distance=-1);
		void writeSimulations(int run, TLinearComb* pcaObject, double distance);
};
class TOutputOneObs:public TOutput{
	public:
	   int obsFileNumber;

	   TOutputOneObs(int gotSampling, int gotObsFileNumber, my_string outputFilename, bool gotWriteDistance, TDataVector* Data, TInputFileVector* InputFiles);
	   ~TOutputOneObs(){myFile.close();};
	   virtual void writeHeader();
	   virtual void writeSimulations(int run, double distance=-1);
};

class TOutputVector{
	public:
		vector<TOutput*> vecOutputFiles;
		vector<TOutput*>::iterator curOutput, endOutput;
		TDataVector* data;
		TInputFileVector* inputFiles;
		bool separateOutputFiles;

		TOutputVector(my_string mcmcsampling, my_string outName, TDataVector* Data, TInputFileVector* InputFiles, bool gotWriteDistance=false, bool gotSeparetOutputFiles=false);
		~TOutputVector(){
			for(unsigned int i=0; i< vecOutputFiles.size();++i){
				delete vecOutputFiles[i];
			}
			vecOutputFiles.clear();
		};                                    
		void writeHeader();
		void writeHeader(TLinearComb* pcaObject);
		void writeSimulations(int run);
		void writeSimulations(int run, double distance);
		void writeSimulations(int run, TLinearComb* pcaObject);
		void writeSimulations(int run, TLinearComb* pcaObject, double distance);             
};






#endif
