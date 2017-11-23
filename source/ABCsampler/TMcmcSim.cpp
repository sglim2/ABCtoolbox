//---------------------------------------------------------------------------

#pragma hdrstop
#include "TMcmcSim.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)

TMcmcSim::TMcmcSim(TParameters* GotParameters, my_string gotexedir, ofstream* gotLogFile):TSim(GotParameters, gotexedir, gotLogFile){
	nbSims=gotParameters->getdoubleParameter("nbSims");
	*logFile << "Number of simulations to perform: " <<  nbSims << endl;

	calibrationFile=gotParameters->getParameter("calName", 0);
	thresholdProportion=gotParameters->getdoubleParameter("tolerance");
	*logFile << "Tolerance is set to: " <<  thresholdProportion << endl;
	mcmc_range_proportion=gotParameters->getdoubleParameter("rangeProp");
	*logFile << "Rangeproportion is set to: " <<  mcmc_range_proportion << endl;

    burninLength=gotParameters->getdoubleParameter("startupLength",0);
	if(!burninLength) burninLength=100;
	*logFile << "Startup Length is set to: " <<  burninLength << endl;
	burninTrials=gotParameters->getdoubleParameter("startupAttempts",0);
	if(!burninTrials) burninTrials=20;
	*logFile << "Max Startup Attempts is set to: " <<  burninTrials << endl;

	stopIfBurninFailed=gotParameters->getdoubleParameter("stopIfStartupFailed",0);

	repeatsPerParameterVector=gotParameters->getdoubleParameter("runsPerParameterVector",0);
	if(!repeatsPerParameterVector) repeatsPerParameterVector=1;
	*logFile << "Number of repeats per Parameter Vector is set to: " <<  repeatsPerParameterVector << endl;

	_Idum=1L;
	/* --> while implemented, there is no theoretical argument for using the distancedensity nor the distance ration....
	if(gotParameters->getParameter("mcmctype", 0)!=""){
		my_string mcmctype=gotParameters->getParameter("mcmctype", 0);
		if(mcmctype.extract_sub_str(*";").toInt()==1) includePriorDens=true; else includePriorDens=false;
		mcmctype.remove(0,1);
		if(mcmctype.extract_sub_str(*";").toInt()==1) includeDistDens=true; else includeDistDens=false;
		mcmctype.remove(0,1);
		if(mcmctype.extract_sub_str(*";").toInt()==1) includeDistRatio=true; else includeDistRatio=false;
	 } else { */
		 includePriorDens=true;
		 includeDistDens=false;
		 includeDistRatio=false;
	 //}
};
//---------------------------------------------------------------------------
void TMcmcSim::performCalibratingSimulations(){
	//This function performs a given number fo simulations and calculates for
	//all summary statistics mean and variance (used to normalize the distances).
	myCaliData=new TSimDatabase((int)gotParameters->getdoubleParameter("numCaliSims"), myData, inputFiles, logFile);
	*logFile << "Calibration is started, performing " << myCaliData->num_sims << " simulations." << endl;
	//Run simulations to fill the arrays
	*logFile << "Start Calibration Loop ......" << endl;
	int nbSimulationsPerformed=0;
	int step_for_logfile=myCaliData->num_sims/50;
	for(int i=0; i< myCaliData->num_sims ; ++i){
		inputFiles->createNewInputFiles(-1);
		if(i%step_for_logfile==0) *logFile << i << "/" << myCaliData->num_sims << endl;
		inputFiles->performSimulations(-1);
		myData->calculateSumStats(-1);
		myCaliData->addSimToDB(myData, inputFiles);
		++nbSimulationsPerformed;
	}
	*logFile << " done!" << endl;
}
//---------------------------------------------------------------------------
void TMcmcSim::initialCalibration(){
	  if(!calibrationFile.empty()){
		  *logFile << "Calibrate the system with already existing simulations ...." << endl;
		  myCaliData=new TSimDatabase(calibrationFile, (int)gotParameters->getdoubleParameter("numCaliSims"), myData, inputFiles, logFile);
		  //myCaliData->writeToFile("calibration_file_2.txt", true);
	  } else {
		  *logFile << "Calibrate the system by performing simulations ...." << endl;
		  performCalibratingSimulations();
		  myCaliData->writeToFile("calibration_file.txt");
	  }
	  myCaliData->calculateMeanVariances();
	  if(doLinearComb) myCaliData->calculateDistances(myLinearComb);
	  else myCaliData->calculateDistances();
	  //if(includeDistDens) myCaliData->fillNormalizedDistances();
	  threshold=myCaliData->getThreshold(thresholdProportion);
	  myCaliData->fillArrayOfSimsBelowThreshold();
	  myCaliData->setMcmcRanges(mcmc_range_proportion);
	  //set starting point
	  if(gotParameters->getParameter("setStartingPoint", 0)=="best")
		 myCaliData->setPriorStartingConditionsFromBestSimulation();
	  else myCaliData->setPriorStartingConditionsAtRandom();
	  *logFile << "Calibration done!" << endl;
	  *logFile << "Threshold is set to: " << threshold << endl;
}
//---------------------------------------------------------------------------
void TMcmcSim::resetCounters(){
	//counters
	nbSimulationsAccepted=0;
	nbSimulationsWhereHAboveOne=0;
	nbSimulationsLargerThreshold=0;
}
//---------------------------------------------------------------------------
void TMcmcSim::performSimulationSeveralRepeats(int numRepeatsPerParameterVector, int s){
   //run simulation program. Run it runPerParameterVector times and report average distance
   distance=0;
   for(int i=0; i<numRepeatsPerParameterVector; ++i){
	  inputFiles->performSimulations();
	  myData->calculateSumStats(s);
	  distance+=calcDistance();
   }
   distance=distance/numRepeatsPerParameterVector;
}
//---------------------------------------------------------------------------
void TMcmcSim::performSimulation(int s){
	  inputFiles->performSimulations();
	  myData->calculateSumStats(s);
	  distance=calcDistance();
}
//---------------------------------------------------------------------------
void TMcmcSim::checkAcceptance(){
   if(distance < threshold){
	  if(calculcateH()){
		 ++nbSimulationsAccepted;
		 inputFiles->priors->saveOldValues();
		 myData->saveOldValues();
		 if(doLinearComb) myLinearComb->saveOldValues();
		 oldDistance=distance;
		 oldDistanceDensity=distanceDensity;
	  } else {
		 inputFiles->priors->resetOldValues();
		 myData->resetOldValues();
		 if(doLinearComb) myLinearComb->resetOldValues();
		 distance=oldDistance;
		 distanceDensity=oldDistanceDensity;
	  }
   } else {
	  inputFiles->priors->resetOldValues();
	  myData->resetOldValues();
	  if(doLinearComb) myLinearComb->resetOldValues();
	  ++nbSimulationsLargerThreshold;
	  distance=oldDistance;
	  distanceDensity=oldDistanceDensity;
   }
}
//---------------------------------------------------------------------------
void TMcmcSim::performBurnin(){
   *logFile << "Startup of " <<  burninLength << " is started ...";
   nbSimulationsAccepted=0;
   for(int s=1; s<=burninLength; ++s, ++nbSimulationsPerformed){
	  // create new input file
	  inputFiles->createNewInputFilesMcmc();
	  //perform simulation
	  if(repeatsPerParameterVector>1) performSimulationSeveralRepeats(repeatsPerParameterVector);
	  else performSimulation();
	  //check if we accept new parameters
	  checkAcceptance();
	  if(s==burninLength){
		 if(nbSimulationsAccepted>1){
			*logFile << "Burnin over!" << endl;
			resetCounters();
		 } else {
			if(burninAttemps>=burninTrials){
			   *logFile << "Burnin failed!" << endl;
			   if(!gotParameters->getdoubleParameter("dontStopIfBurninFailed", 0)) s=0;
			} else {
			   ++burninAttemps;
			   *logFile << endl << burninAttemps << ". attempt!...";
			   myCaliData->setPriorStartingConditionsAtRandom();
			   s=0;
			}
		 }
	  }
   }
}
//---------------------------------------------------------------------------
int TMcmcSim::runSimulations(){
	 *logFile << endl << "Simulations" << endl << "**********" << endl;
	 *logFile << "An MCMC chain is launched ..... !" << endl;
	 //reset Counters
	 resetCounters();
	 nbSimulationsPerformed=0;

	 //open the files for output
	 my_string mcmcSampling=gotParameters->getParameter("mcmcsampling", 0);
	 if(mcmcSampling=="") mcmcSampling="1";
	 cout << outName << endl;
	 myOutputFiles=new TOutputVector(mcmcSampling, outName, myData, inputFiles, addDistanceToOutputfile, writeSeparateOutputFiles);
	 if(doLinearComb) myOutputFiles->writeHeader(myLinearComb);
	 else myOutputFiles->writeHeader();

	 //Do calibrating simulations
	 initialCalibration();

	 //First Simulation after Calibration: get distance etc. for priors chosen
     performSimulation();
	 oldDistance=distance;
	 //if(includeDistDens) oldDistanceDensity=myCaliData->getDistanceDensity(oldDistance);

	 *logFile << "Distance of first Simulation: " << oldDistance << endl;
	 // Loop of Chain, start only if number of sims > 0
	 if(nbSims>0){
		//start burnin
		burninAttemps=0;
		performBurnin();
		//start Loop
		nbSims=nbSims-nbSimulationsPerformed;
		int step_for_logfile=nbSims/50;
		*logFile << "Loop of " << nbSims << " simulations is started ..." << endl;
		for (int s=1; s<=nbSims; ++s, ++nbSimulationsPerformed){
		   if(s%step_for_logfile==0) *logFile << s << "/" << nbSims << endl;
		   // create new input file
		   inputFiles->createNewInputFilesMcmc();
		   //perform simulation
		   if(repeatsPerParameterVector>1) performSimulationSeveralRepeats(repeatsPerParameterVector);
		   else performSimulation(s);
		   //check if we accept new parameters
		   checkAcceptance();
		   if(doLinearComb) myOutputFiles->writeSimulations(s, myLinearComb, distance);
		   else myOutputFiles->writeSimulations(s, distance);
		}
		*logFile << " done!" << endl;
		//write counters to logFile
		*logFile << "Number of bunrin attemps needed: " << burninAttemps+1 << endl;
		*logFile << "Number of performed jumps: " << nbSimulationsPerformed << endl;
		*logFile << "Number of simulations in chain: " << nbSims << endl;
		*logFile << "Number of accepted simulations: " << nbSimulationsAccepted << endl;
		*logFile << "Number of simulations where h>1: " << nbSimulationsWhereHAboveOne << endl;
		*logFile << "Number of simulations where distance larger than threshold: " << nbSimulationsLargerThreshold << endl;
		//output of chain producing new calibration, if performed

	 } //end if nbsims >0
	 delete myOutputFiles;
	 return true;
};

//------------------------------------------------------------------------------
//calculate the acceptance probability of new priors and check, if we accept them
bool TMcmcSim::calculcateH(){
    //as we use a symmetric transition kernel, this part is reduced to the quotient
	double h=inputFiles->getPriorDensity()/inputFiles->getOldPriorDensity();
	hFile << h << endl;
	if(h>1){
		++nbSimulationsWhereHAboveOne;
		return true;
	}
	if(ran3(&_Idum)< h)	return true;
	else return false;
}

//function below not used anymore --> no theoretical argument to use distance densities!
/*
bool TMcmcSim::calculcateH(){
   //as we use a symmetric transition kernel, this part is reduced to the quotient
	//of Pi(new)/Pi(old) times the quotient of the distance densities
	//Pi(oldDistance)/Pi(newDistance) (see NSF Grant ) or related ratios
	double h=1;
	//if we inlcude the prior density (see input File)
	if(includePriorDens){
		double Pi_old=inputFiles->getOldPriorDensity();
		double Pi_new=inputFiles->getPriorDensity();
		h=h* (Pi_new/Pi_old);
	}
	//if we inlcude the distance density (see input file)
	if(includeDistDens){
		double h_density;
		distanceDensity=myCaliData->getDistanceDensity(distance);
		if(oldDistanceDensity==0) h_density=1; else h_density= oldDistanceDensity/distanceDensity;
		h =  h * h_density;
	}
	//if we inlcude the distance ratio (see input file)
	if(includeDistRatio){
		double h_distratio;
		if(distance==0) h_distratio=1; else h_distratio= oldDistance/distance;
		h=h*h_distratio;
	}
   //chek now h
	if(h>1){
		++nbSimulationsWhereHAboveOne;
		return true;
	}
	if(ran3(&_Idum)< h)	return true;
	else return false;
}
*/
//------------------------------------------------------------------------------

