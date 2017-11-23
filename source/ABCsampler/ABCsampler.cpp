//---------------------------------------------------------------------------
#pragma hdrstop
#include "TMcmcSim.h"
#include "TStandardSim.h"
#include "TPmcSim.h"

#ifdef _GCC_
	#include <unistd.h>
	#include <sys/types.h>
	#include <sys/wait.h>
	#include <errno.h>
	#include <time.h>
#else
	#include <process.h>
#endif


//---------------------------------------------------------------------------
#pragma argsused
int main(int argc, char* argv[]){
	ofstream logFile;
	clock_t starttime=clock();
	logFile.open("ABCsampler.log");
	logFile << "ABCsampler Logfile" << endl << "*******************" << endl;
	my_string buf;
   try{
		if (argc<2) throw TException("Missing first argument, the name of an input file!", _FATAL_ERROR);
       //read Inputfile
		TParameters* myParameters;
		logFile << "Reading inputfile '" << argv[1] << "' .....";
		if(argc>2){
			my_string* params=new my_string[argc-2];
			for(int i=2;i<argc;++i) params[i-2]=argv[i];
			myParameters = new TParameters(argv[1], argc-2, params);
		} else myParameters = new TParameters(argv[1]);
		logFile << " done!" << endl;

		//initialize random generator
		long seed=0;
		while(seed==0) seed=myParameters->getdoubleParameter("addToSeed", 0)+get_randomSeedFromCurrentTime();
		if(seed>0) seed=-seed;
		while(seed>161803398) seed=seed/10;
		initialize_ran3(seed);
		logFile << "Random generator initialized with seed " << seed << endl;

		//set directory where this program is launched
		my_string exeDir=argv[0];
		exeDir.extract_path();

		//check which type of simulation to perform
		TSim* mySim;
		bool typeExists=false;
		my_string type=myParameters->getParameter("samplerType");
		if(type=="standard"){
			logFile << "Performing a Standard Simulation" << endl;
			mySim=new TStandardSim(myParameters, exeDir, &logFile);
			typeExists=true;
		}
		if(type=="MCMC"){
			logFile << "Performing an MCMC Simulation" << endl;
			mySim=new TMcmcSim(myParameters, exeDir, &logFile);
			typeExists=true;
		}
		if(type=="PMC"){
			logFile << "Performing an PMC Simulation" << endl;
			mySim=new TPmcSim(myParameters, exeDir, &logFile);
			typeExists=true;
			}
		if(!typeExists) throw TException("Unknown samplerType '"+type+"'!", _FATAL_ERROR);

		//Launch simulation
		mySim->runSimulations();
		my_string unusedParams=myParameters->getListOfUnusedParameters();
		if(unusedParams!=""){
			cout << "\nThe following parameters were not used: " << unusedParams << "!" << endl;
			logFile << "\nThe following parameters were not used: " << unusedParams << "!" << endl;
		}
		delete myParameters;
		delete mySim;
   }
	catch (TException & error){
	  logFile << "\nERROR: " << error.getMessage() << endl;
	  cout <<"\nERROR: " << error.getMessage() << endl;
	}
	catch (...){
      logFile <<"\nERROR: unhandeled error!" << "\nProgram is finished!"  << endl;
	  cout <<"\nERROR: unhandeled error!" << "\nProgram is finished!"  << endl;
	}
	logFile << "End of Program!!" << endl;
	logFile.close();
	cout << "\nProgram terminated in " << (float) (clock()-starttime)/CLOCKS_PER_SEC/60<< "min!\n";

   return 0;
}


//---------------------------------------------------------------------------

