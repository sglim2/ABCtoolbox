//---------------------------------------------------------------------------


#include "TStandardSim.h"
//#include <sys/wait.h>
//---------------------------------------------------------------------------

TStandardSim::TStandardSim(TParameters* GotParameters, my_string gotexedir, ofstream* gotLogFile):TSim(GotParameters, gotexedir, gotLogFile){
	nbSims=gotParameters->getdoubleParameter("nbSims");
	*logFile << "Number of simulations to perform: " <<  nbSims << endl;
	repeatsPerParameterVector=gotParameters->getdoubleParameter("runsPerParameterVector",0);
	if(!repeatsPerParameterVector) repeatsPerParameterVector=1;
	*logFile << "Number of repeats per parameter vector is set to: " <<  repeatsPerParameterVector << endl;


};

//---------------------------------------------------------------------------
int TStandardSim::runSimulations(){
		myOutputFiles=new TOutputVector("1", outName, myData, inputFiles, false, writeSeparateOutputFiles);
		if(gotParameters->getParameter("writeHeader", 0)=="1") myOutputFiles->writeHeader();
		*logFile << endl << "Simulations" << endl << "**********" << endl;
		*logFile << "Loop is started....." << endl;

		double distance;
		for (int s=0; s<nbSims; ++s){
				if(s%50==0) *logFile << s << "/" << nbSims << endl;
				// create new input file
				inputFiles->createNewInputFiles(s);
				for(int i=0; i<repeatsPerParameterVector; ++i){
					// perform simulations
					inputFiles->performSimulations(s);
					//calculate SumStats
					myData->calculateSumStats(s);
					//write Simulation
					if(addDistanceToOutputfile){
						distance=calcDistance();
						myOutputFiles->writeSimulations(s, distance);
					} else myOutputFiles->writeSimulations(s);
				}
		}
		// close the files
		delete myOutputFiles;
      return true;
};
