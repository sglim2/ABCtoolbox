#ifndef TStandardSimH
#define TStandardSimH
//---------------------------------------------------------------------------

//#include "global.h"
#include "TSim.h"
#include "TParameters.h"


//---------------------------------------------------------------------------
class TStandardSim:public TSim{
   public:
	   int nbSims;
		double distance;
		int repeatsPerParameterVector;
		TOutputVector* myOutputFiles;
        virtual int runSimulations();

      TStandardSim(TParameters* gotParameters, my_string gotexedir, ofstream* gotLogFile);
      ~TStandardSim(){}

};
#endif
