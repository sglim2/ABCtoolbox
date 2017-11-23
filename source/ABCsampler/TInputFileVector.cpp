//---------------------------------------------------------------------------

#pragma hdrstop

#include "TInputFileVector.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
TInputFileVector::TInputFileVector(TParameters* gotParameters, ofstream* gotLogFile){
	logFile=gotLogFile;

	//est file=file with prior definitions
	estFile=gotParameters->getParameter("estName");
	*logFile << "Name of file with prior definitions: " <<  estFile << endl;
	priors = new TPriorVector(estFile, logFile);

	//create input files
	my_string siminputname=gotParameters->getParameter("simInputName");
	while(!siminputname.empty()){
	   vecInputFileObjects.push_back(TInputFile(siminputname.extract_sub_str(*";"), priors));
	   siminputname.remove(0,1);
	}
	endInputFileObject=vecInputFileObjects.end();
	//write names to logfile
	*logFile << "Names of the simulation input files:\n";
	curInputFileObject=vecInputFileObjects.begin();
	for(;curInputFileObject!=endInputFileObject; ++curInputFileObject){
	   curInputFileObject->writeNamesToLogfile(logFile);
	}
	*logFile << "Corresponding names of the temporal simulation input files:\n";
	curInputFileObject=vecInputFileObjects.begin();
	for(;curInputFileObject!=endInputFileObject; ++curInputFileObject){
	   curInputFileObject->writeNewNamesToLogfile(logFile);
	}

	//read simulation programs...
	vector<my_string> simProgList;
	my_string temp=gotParameters->getParameter("simulationProgram");
	while(!temp.empty()){
		simProgList.push_back(temp.extract_sub_str(*";"));
		temp.remove(0,1);
	}

	//read parameters for the simulation program.
	//these may be strings, priors or the tag FILENAME, which prints the same of the inputfile (*-temp.par)
    //the TExecuteProram will be a variable of TInputFile
	//first split with ";" to get for each inputfile
    vector<my_string> simParamList;
    temp=gotParameters->getParameter("simParam");
	while(!temp.empty()){
		simParamList.push_back(temp.extract_sub_str(*";"));
		temp.remove(0,1);
	}
	//check number of entries and initialize the simulation program in each TInputfile. three options:
	//- one program, one argument list
	//- one program, several argument lists
	//- several programs, several argument lists
	curInputFileObject=vecInputFileObjects.begin();
    if(simParamList.size()==1 && simProgList.size()==1){
    	*logFile << "Using simulation program: '" <<  simProgList[0] << "'" << endl;
    	for(;curInputFileObject!=endInputFileObject; ++curInputFileObject){
    	   curInputFileObject->initializeSimulationProgramm(simProgList[0], simParamList[0]);
    	}
    } else {
    	if(simParamList.size()==vecInputFileObjects.size()){
    		if(simProgList.size()==1){
    			*logFile << "Using simulation program: '" <<  simProgList[0] << "'" << endl;
				for(int i=0;i<simParamList.size(); ++i, ++curInputFileObject){
					curInputFileObject->initializeSimulationProgramm(simProgList[0], simParamList[i]);
				}
    		} else {
    			if(simProgList.size()==simParamList.size()){
    				*logFile << "Using the following simulation programs:" <<  endl;
    				for(int i=0;i<simParamList.size(); ++i, ++curInputFileObject){
    					*logFile << "   - '" << simProgList[i] << "'" << endl;
    					curInputFileObject->initializeSimulationProgramm(simProgList[i], simParamList[i]);
    				}
    			} else throw TException("The number of simulation programs and the number of argument lists for the simulation program do not match!", _FATAL_ERROR);
    		}
    	} else throw TException("Unequal number of input files and argument lists for the simulation program!", _FATAL_ERROR);
    }


	//script/program before and after each simulation
	#ifdef _GCC_
    //script before
    my_string scriptbeforesimulation=gotParameters->getParameter("launchBeforeSim", 0);
    if(scriptbeforesimulation!=""){
    	//check what arguments have to be passed
    	my_string temp=gotParameters->getParameter("launchBeforeSimParam", false);
    	//split for different inputfiles
    	vector<my_string> scriptParamList;
    	while(!temp.empty()){
    		scriptParamList.push_back(temp.extract_sub_str(*";"));
    		temp.remove(0,1);
    	}
    	if(scriptParamList.size()==0){
    		//initialize simulation program in each TInputfile
    	    curInputFileObject=vecInputFileObjects.begin();
    	    for(;curInputFileObject!=endInputFileObject; ++curInputFileObject){
    	    	curInputFileObject->initializeScriptBeforeSimulation(scriptbeforesimulation, "");
    	    }
    	} else {
    		if(scriptParamList.size()==1){
    			//initialize simulation program in each TInputfile
    			curInputFileObject=vecInputFileObjects.begin();
    			for(;curInputFileObject!=endInputFileObject; ++curInputFileObject){
    				curInputFileObject->initializeScriptBeforeSimulation(scriptbeforesimulation, scriptParamList[0]);
    			}
    		} else {
    			if(scriptParamList.size()==vecInputFileObjects.size()){
    				curInputFileObject=vecInputFileObjects.begin();
    				for(int i=0;i<scriptParamList.size(); ++i, ++curInputFileObject){
    					curInputFileObject->initializeScriptBeforeSimulation(scriptbeforesimulation, scriptParamList[i]);
    				}
    			} else throw TException("Unequal number of input files and parameter specifications for the script/program to launch before the simulations!", _FATAL_ERROR);
    		}
    	}
    	launchScriptBeforeSimulation=true;
        *logFile << "Before each simulation, the following script is launched: " <<  scriptbeforesimulation << endl;
    } else launchScriptBeforeSimulation=false;

   	//script after
	my_string scriptaftersimulation=gotParameters->getParameter("launchAfterSim", 0);
	if(scriptaftersimulation!=""){
		//check what arguments have to be passed
		my_string temp=gotParameters->getParameter("launchAfterSimParam", false);
		//split for different inputfiles
		vector<my_string> scriptParamList;
		while(!temp.empty()){
			scriptParamList.push_back(temp.extract_sub_str(*";"));
			temp.remove(0,1);
		}
		if(scriptParamList.size()==0){
			//initialize simulation program in each TInputfile
		    curInputFileObject=vecInputFileObjects.begin();
		    for(;curInputFileObject!=endInputFileObject; ++curInputFileObject){
		       curInputFileObject->initializeScriptAfterSimulation(scriptaftersimulation, "");
		    }
		} else {
			if(scriptParamList.size()==1){
				//initialize simulation program in each TInputfile
				curInputFileObject=vecInputFileObjects.begin();
				for(;curInputFileObject!=endInputFileObject; ++curInputFileObject){
					curInputFileObject->initializeScriptAfterSimulation(scriptaftersimulation, scriptParamList[0]);
				}
			} else {
				if(scriptParamList.size()==vecInputFileObjects.size()){
					curInputFileObject=vecInputFileObjects.begin();
					for(int i=0;i<scriptParamList.size(); ++i, ++curInputFileObject){
						curInputFileObject->initializeScriptAfterSimulation(scriptaftersimulation, scriptParamList[i]);
					}
				} else throw TException("Unequal number of input files and parameter specifications for the script to launch after the simulations!", _FATAL_ERROR);
			}
		}
		launchScriptAfterSimulation=true;
	   *logFile << "After each simulation, the following script is launched: " <<  scriptaftersimulation << endl;
	} else launchScriptAfterSimulation=false;
	#else
	   launchScriptAfterSimulation=false;
	   launchScriptBeforeSimulation=false;
    #endif


}
//---------------------------------------------------------------------------
void TInputFileVector::writeNewInputFiles(int simnum){
   //write files
   curInputFileObject=vecInputFileObjects.begin();
   for(;curInputFileObject!=endInputFileObject; ++curInputFileObject){
	   //launch script before simulations
	   if(launchScriptBeforeSimulation) curInputFileObject->launchScriptBeforeSimulations(simnum);
	   curInputFileObject->createNewInputFile();
   }
}
//---------------------------------------------------------------------------
void TInputFileVector::getNewValueMcmc(const int& priorNumber){
   priors->getNewValuesMcmc(priorNumber);
}
//---------------------------------------------------------------------------
void TInputFileVector::createNewInputFiles(int simnum){
   priors->getNewValues();
   writeNewInputFiles(simnum);
}
void TInputFileVector::createNewInputFilesMcmc(int simnum){
   priors->getNewValuesMcmc();
   writeNewInputFiles(simnum);
}
void TInputFileVector::createNewInputFilesMcmcUpdateOnePriorOnly(int simnum){
   priors->getNewValuesMcmcUpdateOnePriorOnly();
   writeNewInputFiles(simnum);
}
bool TInputFileVector::creatNewInputFilesPMC(double* newParams, int simnum){
	if(priors->getNewValuesPMC(newParams)){
		writeNewInputFiles(-1);
		return true;
	}
	return false;
}
//---------------------------------------------------------------------------
void TInputFileVector::performSimulations(int simnum){
   curInputFileObject=vecInputFileObjects.begin();
   for(;curInputFileObject!=endInputFileObject; ++curInputFileObject){
	  if(!curInputFileObject->performSimulation()){
		 *logFile << "Simulation program returned an error! The following parameters were used:" << endl << "Parameter names:";
		 priors->writeHeader(*logFile);
		 *logFile << endl << "Values:";
		 priors->writeParameters(*logFile);
		 *logFile << endl;
		 //throw TException("Simulation program returned an error!", _FATAL_ERROR);
	  }
      if(launchScriptAfterSimulation) curInputFileObject->launchScriptAfterSimulations(simnum);
   }
}
//---------------------------------------------------------------------------
vector<my_string> TInputFileVector::getVectorOfInputFileNames(){
	//return a vector with all inputFileNames that are used to execute
	//the simulation program (-> temporary file names!)
	//if several files are given per Inputfile only the first one is returned....
	vector<my_string> filenames;
	curInputFileObject=vecInputFileObjects.begin();
	for(;curInputFileObject!=endInputFileObject; ++curInputFileObject){
		filenames.push_back(curInputFileObject->newSimulationProgramInputFilenames[0]);
	}
	return filenames;
}
//---------------------------------------------------------------------------
double TInputFileVector::getPriorDensity(){
   return priors->getPriorDensity();
}
double TInputFileVector::getOldPriorDensity(){
   return priors->getOldPriorDensity();
}
double TInputFileVector::getPriorDensity(double* values){
	return priors->getPriorDensity(values);
}
//---------------------------------------------------------------------------
int TInputFileVector::getNumberOfSimplePriorFromName(my_string name){
   return priors->getNumberOfSimplePriorFromName(name);
}
//---------------------------------------------------------------------------
void TInputFileVector::setPriorValue(int priorNumber, double value){
   priors->setSimplePriorValue(priorNumber, value);
}
//---------------------------------------------------------------------------
void TInputFileVector::updateCombinedParameters(){
   priors->updateCombinedParameters();
}
