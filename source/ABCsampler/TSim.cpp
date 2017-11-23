//---------------------------------------------------------------------------


#include "TSim.h"
#ifdef _GCC_
   #include <unistd.h>
   #include <sys/types.h>
   #include <sys/wait.h>
   #include <errno.h>
#else
	#include <process.h>
  //	#include <sys/wait.h>
#endif
//---------------------------------------------------------------------------
TSim::TSim(TParameters* GotParameters, my_string gotexedir, ofstream* gotLogFile){
  logFile=gotLogFile;
  gotParameters=GotParameters;
  exeDir=gotexedir;

  //read parameters from inputfile and write logfile
  *logFile << endl << "Filenames" << endl << "**********" << endl;
  outName=gotParameters->getParameter("outName");
  *logFile << "Prefix for output files: " <<  outName << endl;
  //create input Files
  inputFiles = new TInputFileVector(gotParameters, gotLogFile);

  //create first input file and perform first simulation
  *logFile << endl << "Initialization" << endl << "**********" << endl;
  *logFile << "create first input file and initialize ...." << endl;
  inputFiles->createNewInputFiles();
  inputFiles->priors->saveOldValues();
  inputFiles->performSimulations();
  //create Data Object
  myData= new TDataVector(gotParameters, inputFiles, gotLogFile);
  myData->initializeDataObjects();

  //linear tranformation to statistics applied?
  linearCombFileName=gotParameters->getParameter("linearCombName",0);
  if(!linearCombFileName.empty()){
	  if(writeSeparateOutputFiles) throw TException("Linear transformation can not be used if several output files are written!", _FATAL_ERROR);
	  doLinearComb=true;
	  *logFile << "Linear combinations are computed as defined in the file: " <<  linearCombFileName << endl;
	  doBoxCox=gotParameters->getdoubleParameter("doBoxCox",0);
	  if(doBoxCox) *logFile << "All statistics are Box-Cox tranformed before linear combinations are computed, as defined in the File: " <<  linearCombFileName << endl;
	  //initialize pcaObject
	  createLinearCombObject();
  } else doLinearComb=false;

  //output files
	if(gotParameters->getdoubleParameter("addDistanceToOutputfile", 0)){
		addDistanceToOutputfile=true;
		*logFile << "The Distance is added to the output files." << endl;
	 } else addDistanceToOutputfile=false;
	if(gotParameters->getdoubleParameter("separateOutputFiles", 0)){
		writeSeparateOutputFiles=true;
		*logFile << "A separate Outputfile is written for every obs file." << endl;
	} else {
		writeSeparateOutputFiles=false;
		*logFile << "All Statistics are written in one single Outputfile." << endl;
	}

};
//---------------------------------------------------------------------------
void TSim::createLinearCombObject(){
	//create vector with the names from the input file, therefore the names of the simulated stats in the right order
	vector<my_string> statNames=myData->getNamesVector();
	vector<double> obsData=myData->getObsDataVector();
	//double* obsData;
	//myData->fillObsDataArray(obsData);
	if(doBoxCox){
		myLinearComb=new TLinearCombBoxCox(linearCombFileName, statNames);
	}
	else myLinearComb=new TLinearComb(linearCombFileName, statNames);
	myLinearComb->calcObsLineraComb(obsData);
}
//---------------------------------------------------------------------------
double TSim::calcDistance(){
   if(doLinearComb) return myData->calculateDistance(myLinearComb);
   else return myData->calculateDistance();
}
//---------------------------------------------------------------------------
int TSim::runSimulations(){
   return 0;
};
//---------------------------------------------------------------------------



