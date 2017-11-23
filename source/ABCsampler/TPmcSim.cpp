//---------------------------------------------------------------------------

#pragma hdrstop
#include "TPmcSim.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)

TPmcSim::TPmcSim(TParameters* GotParameters, my_string gotexedir, ofstream* gotLogFile):TSim(GotParameters, gotexedir, gotLogFile){
	calibrationFile=gotParameters->getParameter("calName", 0);

	numInterations=gotParameters->getdoubleParameter("numInterations");
	*logFile << "The number of iterations is set to: " <<  numInterations << endl;
	particleSize=gotParameters->getdoubleParameter("sampleSize");
	*logFile << "The sample size is set to: " <<  particleSize << endl;
	lastParticleSize=(int)gotParameters->getdoubleParameter("lastSampleSize", false);
	if(lastParticleSize<=0) lastParticleSize=particleSize;
	*logFile << "The sample size of the last iteration is set to: " <<  lastParticleSize << endl;
	tolerance=gotParameters->getdoubleParameter("tolerance");
	numRetained=tolerance*particleSize;
	if(numRetained<1) throw TException("Given the defined threshold, no simulations is retained!", _FATAL_ERROR);
	if(numRetained>particleSize) throw TException("Given the defined threshold, all simulations are retained!", _FATAL_ERROR);
	*logFile << "Tolerance is set to: " <<  tolerance << ". " << numRetained << " simulations will be used to generate the next population." << endl;
	pmc_range_proportion=gotParameters->getdoubleParameter("rangeProp", false);
	if(pmc_range_proportion<=0) pmc_range_proportion=2;
	*logFile << "Range proportion is set to: " <<  pmc_range_proportion << endl;
	numCaliSims=(int)gotParameters->getdoubleParameter("numCaliSims");
	if(numCaliSims<numRetained) throw TException("Less calibration simulations than simulations to retain in the first iteration!", _FATAL_ERROR);
	repeatsPerParameterVector=gotParameters->getdoubleParameter("runsPerParameterVector",0);
	if(!repeatsPerParameterVector) repeatsPerParameterVector=1;
	*logFile << "Number of repeats per parameter vector is set to: " <<  repeatsPerParameterVector << endl;

	_Idum=1L;
	//prepare the two sim databases
	mySimData=new TSimDatabase*[3];
	weights=new double*[3];
	weights[0]=new double[numRetained];
	weights[1]=new double[numRetained];
	weights[2]=new double[numRetained];
};
//---------------------------------------------------------------------------
void TPmcSim::performCalibratingSimulations(){
	//This function performs a given number of simulations and calculates for
	//all summary statistics mean and variance (used to normalize the distances).
	mySimData[0]=new TSimDatabase(numCaliSims, myData, inputFiles, logFile);

	*logFile << "Calibration is started, performing " << mySimData[0]->num_sims << " simulations." << endl;
	//Run simulations to fill the arrays
	*logFile << "Start Calibration Loop ......" << endl;
	for(int i=0; i< mySimData[0]->num_sims ; ++i){
		inputFiles->createNewInputFiles(-1);
		inputFiles->performSimulations(-1);
		myData->calculateSumStats(-1);
		mySimData[0]->addSimToDB(myData, inputFiles);
	}
	*logFile << " done!" << endl;
}
//---------------------------------------------------------------------------
void TPmcSim::initialCalibration(){
	  if(!calibrationFile.empty()){
		  *logFile << "Calibrate the system with already existing simulations ...." << endl;
		  mySimData[0]=new TSimDatabase(calibrationFile, numCaliSims, myData, inputFiles, logFile);
	  } else {
		  *logFile << "Calibrate the system by performing simulations ...." << endl;
		  performCalibratingSimulations();
		  mySimData[0]->writeToFile("calibration_file.txt");
	  }

	  //calc things for calibration set
	  mySimData[0]->calculateMeanVariances();
	  if(doLinearComb) mySimData[0]->calculateDistances(myLinearComb);
	  else mySimData[0]->calculateDistances();
	  *logFile << "Calibration done!" << endl;
}
//---------------------------------------------------------------------------
int TPmcSim::runSimulations(){
	 *logFile << endl << "Simulations" << endl << "**********" << endl;
	 *logFile << "An ABC-PMC algorithm is launched ..... !" << endl;

	 initialCalibration();

	 //set up the second database
	 mySimData[1]=new TSimDatabase(particleSize, myData, inputFiles, logFile);
	 mySimData[2]=new TSimDatabase(particleSize, myData, inputFiles, logFile);
	 nextDb=1;
	 currentDb=0;
	 oldDb=2;

	 ColumnVector vij(inputFiles->priors->numSimplePrior);
	 ColumnVector z(inputFiles->priors->numSimplePrior);
	 double* cumulativeWeights=new double[numRetained];
	 double* newParams=new double[inputFiles->priors->numSimplePrior];

	 for(int iteration=0; iteration<numInterations; ++iteration){
		 //make rejection
		 mySimData[currentDb]->setThreshold(numRetained);
		 mySimData[currentDb]->fillArrayOfSimsBelowThreshold();
		 *logFile << "The threshold in iteration " << iteration << " is set to: " << mySimData[currentDb]->threshold << endl;

		 //standardize the retained parameters
		 if(iteration==0) mySimData[currentDb]->calculateMinMaxofParameters(paramMin, paramMax);
		 mySimData[currentDb]->standardizeRetainedParameters(paramMin, paramMax);

		 //set weights
		 if(iteration==0){
			 for(int i=0;i<numRetained;++i) weights[nextDb][i]=(double)1/(double) numRetained;
		 } else {
			 double weightSum=0;
			 for(int i=0;i<numRetained;++i){
				 double sum=0;
				 for(int j=0;j<numRetained;++j){
					 //calculate vij, which is the column vector of the differences between the particle i and j
					 for(int n=0; n<inputFiles->priors->numSimplePrior; ++n){
					 	vij.element(n)=mySimData[currentDb]->standardizedParameters[i][n]-mySimData[oldDb]->standardizedParameters[j][n];
					 }
					 Matrix temp=-0.5*vij.t()*Sigma_inv*vij;
					 sum+=weights[currentDb][j]*exp(temp.element(0,0));
				 }
				 weights[nextDb][i]=mySimData[currentDb]->getPriorDensityOneRetainedSim(i)/sum;
				 weightSum+=weights[nextDb][i];
			 }
			 for(int i=0;i<numRetained;++i) weights[nextDb][i]=weights[nextDb][i]/weightSum;
		 }
		 cumulativeWeights[0]=weights[nextDb][0];
		 for(int i=1;i<numRetained;++i) cumulativeWeights[i]=cumulativeWeights[i-1]+weights[nextDb][i];

		 //calculate Sigma
		 mySimData[currentDb]->calculateWeightedSigmaOfRetainedParameters(Sigma, weights[nextDb]);
		 //*logFile << Sigma << endl;
		 Sigma=pmc_range_proportion*Sigma;
		 Sigma_inv=Sigma.i();

		 //DiagonalMatrix D;
		 //Jacobi(Sigma, D);
		 //*logFile << D << endl << endl;
		 try{
			 A=Cholesky(Sigma);
		 } catch (...){
		 	throw TException("Problems solving the Cholesky decomposition of Sigma!", _FATAL_ERROR);
		 }
		 //generate new population
		 if(iteration==numInterations-1) particleSize=lastParticleSize;
		 mySimData[nextDb]->empty(particleSize);
		 myOutputFiles=new TOutputVector("1", outName+"_iteration_"+iteration, myData, inputFiles, false, writeSeparateOutputFiles);
		 if(doLinearComb) myOutputFiles->writeHeader(myLinearComb);
		 else myOutputFiles->writeHeader();
		 for(int j=0; j<particleSize;++j){
			 double distance=0;
			 for(int r=0; r<repeatsPerParameterVector;++r){
				 bool parametersAccepted=false;
				 while(!parametersAccepted){
					 double rand=UniformRandom(0, 1);
					 //find sim according to weights....
					 int thisSim=0;
					 while(cumulativeWeights[thisSim]<rand) ++thisSim;
					 ColumnVector mu(inputFiles->priors->numSimplePrior);
					 for(int i=0;i<inputFiles->priors->numSimplePrior;++i) mu.element(i)=mySimData[currentDb]->standardizedParameters[thisSim][i];
					 //get new multivariate normal values
					 for(int i=0;i<inputFiles->priors->numSimplePrior;++i) z.element(i)=NormalRandom(0,1);
					 //get new parameters
					 z=mu+A*z;
					 //now destandardize them....
					 for(int i=0;i<inputFiles->priors->numSimplePrior;++i) newParams[i]=(z.element(i)*(paramMax[i]-paramMin[i]))+paramMin[i];
					 parametersAccepted=inputFiles->getPriorDensity(newParams);
					 if(parametersAccepted) parametersAccepted=inputFiles->creatNewInputFilesPMC(newParams, j);
				 }
				 inputFiles->performSimulations(j);
				 myData->calculateSumStats(j);
				 distance+=calcDistance();
			 }
			 distance=distance/repeatsPerParameterVector;
			 if(addDistanceToOutputfile) myOutputFiles->writeSimulations(j, distance);
			 else myOutputFiles->writeSimulations(j);
			 mySimData[nextDb]->addSimToDB(myData, inputFiles, distance);
		 }
		 delete myOutputFiles;
		 int temp=currentDb;
		 currentDb=nextDb;
		 nextDb=oldDb;
		 oldDb=temp;
	 }
	 return true;
};
