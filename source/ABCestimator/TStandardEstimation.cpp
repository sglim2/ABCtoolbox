//---------------------------------------------------------------------------

#pragma hdrstop

#include "TStandardEstimation.h"
#pragma package(smart_init)
//---------------------------------------------------------------------------

TStandardEstimation::TStandardEstimation(TParameters* parameters):TEstimation(parameters){
	trueParamName=parameters->getParameter("trueParamName", false);
	if(!trueParamName.empty()){
		trueParamsAvailable=true;
		myObsData->readTrueFile(trueParamName);
		//are the parameters with true values estimated?
		trueParamNumInSimData=new int[myObsData->numTrueParams];
		for(int i=0; i<myObsData->numTrueParams; ++i){
			trueParamNumInSimData[i]=mySimData->getParamNumberFromName(myObsData->trueNameVector[i]);
			if(trueParamNumInSimData[i]<0) throw TException("The parameter '"+myObsData->trueNameVector[i]+"' is not estimated, but requested from the file with true parameter values!", _FATAL_ERROR);
		}
	} else trueParamsAvailable=false;
}
//---------------------------------------------------------------------------
void TStandardEstimation::standardize(){
	TEstimation::standardize();
    if(standardizeStats && distFromDifferentStats){
    	cout << "Standardizing statistics used to calculate the distances..." << endl;
    	mySimDataForDistForDist->standardizeStatistics();
    	myObsDataForDist->standardizeObservedValues(mySimDataForDistForDist->statMeans, mySimDataForDistForDist->statSDs);
    	cout << " done!" << endl;
    }
}
//---------------------------------------------------------------------------
void TStandardEstimation::performEstimations(){
   cout << endl << "ABC-GLM Estimations:" << endl << "------------" << endl;
   //standardize
   standardize();
   preparePosteriorDensityPoints();

   //open file for quantiles
   ofstream quantileFile, accuracyFileRMISE, accuracyFileREMODE, accuracyFileREMEAN, accuracyFileREMEDIAN;
   if(trueParamsAvailable){
	   my_string filename=outputPrefix + "quantilesOfTrueParameters.txt";
	   quantileFile.open(filename.c_str());
	   filename=outputPrefix + "accuracyMeasuresRRMISE.txt";
	   accuracyFileRMISE.open(filename.c_str());
	   filename=outputPrefix + "accuracyMeasuresREMODE.txt";
	   accuracyFileREMODE.open(filename.c_str());
	   filename=outputPrefix + "accuracyMeasuresREMEAN.txt";
	   accuracyFileREMEAN.open(filename.c_str());
	   filename=outputPrefix + "accuracyMeasuresREMEDIAN.txt";
	   accuracyFileREMEDIAN.open(filename.c_str());
	   //write header
	   quantileFile << "observedSet";
	   accuracyFileRMISE << "observedSet";
	   accuracyFileREMODE << "observedSet";
	   accuracyFileREMEAN << "observedSet";
	   accuracyFileREMEDIAN << "observedSet";
	   myObsData->writeTrueParamsheader(&quantileFile);
	   myObsData->writeTrueParamsheader(&accuracyFileRMISE);
	   myObsData->writeTrueParamsheader(&accuracyFileREMODE);
	   myObsData->writeTrueParamsheader(&accuracyFileREMEAN);
	   myObsData->writeTrueParamsheader(&accuracyFileREMEDIAN);
	   quantileFile << endl;
	   accuracyFileRMISE << endl;
	   accuracyFileREMODE << endl;
	   accuracyFileREMEAN << endl;
	   accuracyFileREMEDIAN << endl;
   }

   //loop for all observed data sets
   for(int i=0;i<myObsData->numObsDataSets;++i){
	  cout << "performing estimations for observed data set " << i << ":" << endl;
	  //make rejection
	  makeRejection(i);
	  //check if a statistic is monomorphic
	  my_string monomorphicStat=mySimData->checkIfstatsArePolymorphic();
	  if(monomorphicStat!="") cout << "   - the statistic '" << monomorphicStat << "' is monomorphic among the retained simulations -> Skipping this dataset!" << endl;
	  else {
		  //prepare Dirac peaks -> also used for writing the prior!
		  prepareDiracPeaks();
		  //prepare prior
		  preparePrior();
		  //perform linear Regression
		  myObsData->getObsValuesIntoColumnVector(i, obsValues);
		  performLinearRegression();
		  //calculate posterior
		  calculatePosterior();
		  //calculate marginal density
		  printMarginalDensity(i);
		  //calculate smoothed used sims
		  calculateSmoothedRetainedMatrix();
		  //calculate common surface
		  cout << "   - calculate common surface of prior and posterior ...";
		  propCommonSurfacePriorPosterior << "Obs" << i;
		  for(int p=1; p<=mySimData->numParams;++p) propCommonSurfacePriorPosterior << "\t" << getCommonSurfaceRetainedPosterior(p);
		  propCommonSurfacePriorPosterior << endl;
		  cout << "done!" << endl;
		  //write posterior
		  writePosteriorFile("_Obs"+(my_string)i);
		  writePosteriorCharacteristics("_Obs"+(my_string)i);
		  //write prior
		  if(writeSmoothedsimsUsed) writeSmoothedSimsUsed("_Obs"+(my_string)i);
		  //if requested, calculate quantiles of true parameter within cumulative posterior distribution
		  if(trueParamsAvailable){
			  quantileFile << i;
			  accuracyFileRMISE << i;
			  accuracyFileREMODE << i;
			  accuracyFileREMEAN << i;
			  accuracyFileREMEDIAN << i;
			  for(int tp=0; tp<myObsData->numTrueParams; ++tp){
				  quantileFile << "\t" << getquantileOfTrueParam(trueParamNumInSimData[tp]+1, myObsData->getTrueParam(i, tp));
				  accuracyFileRMISE << "\t" << getRMISEOfTrueParam(trueParamNumInSimData[tp]+1, myObsData->getTrueParam(i, tp));
				  accuracyFileREMODE << "\t" << getREMODEOfTrueParam(trueParamNumInSimData[tp]+1, myObsData->getTrueParam(i, tp));
				  accuracyFileREMEAN << "\t" << getREMEANOfTrueParam(trueParamNumInSimData[tp]+1, myObsData->getTrueParam(i, tp));
				  accuracyFileREMEDIAN << "\t" << getREMEDIANOfTrueParam(trueParamNumInSimData[tp]+1, myObsData->getTrueParam(i, tp));
			  }
			  quantileFile << endl;
			  accuracyFileRMISE << endl;
			  accuracyFileREMODE << endl;
			  accuracyFileREMEAN << endl;
			  accuracyFileREMEDIAN << endl;
		  }
	  }
	  cout << "   - finished estimations for data set " << i << "!" << endl;
   }
   if(trueParamsAvailable){
	   quantileFile.close();
	   accuracyFileRMISE.close();
	   accuracyFileREMODE.close();
	   accuracyFileREMEAN.close();
	   accuracyFileREMEDIAN.close();
   }
   mdFile.close();
   propCommonSurfacePriorPosterior.close();
}
//---------------------------------------------------------------------------
void TStandardEstimation::makeRejection(int thisObsDataSet){
	//prepare distance array
	float* distances=new float[mySimData->numReadSims];
	//first we have to calculate the distances
	if(distFromDifferentStats){
		mySimDataForDistForDist->calculateDistances(myObsDataForDist->getObservedValues(thisObsDataSet), distances);
		if(thisObsDataSet==myObsData->numObsDataSets-1) delete mySimDataForDistForDist;
	} else {
		mySimData->calculateDistances(myObsData->getObservedValues(thisObsDataSet), distances);
	}
	//keep only retained -> check how they are defined
	if(thresholdDefined){
		if(retainedDefined) mySimData->fillStatAndParamMatricesFromThreshold(distances, threshold, numToRetain);
		else mySimData->fillStatAndParamMatricesFromThreshold(distances, threshold);
	} else {
		if(retainedDefined) mySimData->fillStatAndParamMatricesFromNumRetained(distances, numToRetain);
		else mySimData->fillStatAndParamMatricesKeepAll();
	}
	//check if at least 10 simulations has been retained
	if(mySimData->numUsedSims<10) throw TException("Given the current specified tolerance / numRetained, less than 10 (only "+ (my_string) mySimData->numUsedSims + ")simulation are retained!", _FATAL_ERROR);
	//write retained, if requested
	if(writeRetainedSimulations) mySimData->writeRetained(outputPrefix, distances, thisObsDataSet);
	delete[] distances;
}
//---------------------------------------------------------------------------
void TStandardEstimation::printMarginalDensity(int thisObsDataSet){
   double f_M=0.0;
   Matrix D = SigmaS+CHat*SigmaTheta*CHat.t();
   Matrix D_inv=D.i();
   ColumnVector theta_j, m_j;

   if(calcObsPValue) calculateMarginalDensitiesOfRetained();

   //loop over all dirac peaks
   for(int d=1; d<=mySimData->numUsedSims;++d){
	  theta_j=mySimData->paramMatrix.submatrix(d,d,2,mySimData->numParams+1).as_column();
	  m_j=c_zero+CHat*theta_j;
	  Matrix temp=-0.5*(obsValues-m_j).t()*D_inv*(obsValues-m_j);
	  f_M+=exp(temp(1,1));
   }
   double determinant=(2*3.14159265358979323846*D).determinant();
   if(determinant==0) cout << "\nWARNING: Problems calculating the marginal density: determinant is zero!";
  f_M=f_M/(mySimData->numUsedSims*sqrt(determinant));
  f_M=f_M*((double)mySimData->numUsedSims/(double)mySimData->numReadSims);
  if(calcObsPValue){
	  float p=0;
	  for(int i=0;i<calcObsPValue;++i) if(fmFromRetained[i]<=f_M) ++p;
	  cout << "   - marginal density: " << f_M << "\t(p-value " << p/(float)calcObsPValue << ")" << endl;
      mdFile << thisObsDataSet << "\t" << f_M << "\t" << p/(float)calcObsPValue << endl;
  } else {
	  cout << "   - marginal density: " << f_M << endl;
	  mdFile << thisObsDataSet << "\t" << f_M << endl;
  }
}
//---------------------------------------------------------------------------
float TStandardEstimation::getquantileOfTrueParam(int param, float trueValue){
	float totalArea=getArea(posteriorMatrix.column(param), theta);
	float area=0;
	if(theatParameterScale(param, 1)>trueValue) return 0;
	if(theatParameterScale(param, posteriorDensityPoints)<trueValue) return 1;
	for(int k=2;k<=posteriorDensityPoints;++k){
		if(theatParameterScale(param, k)<trueValue) area+=step*(posteriorMatrix(k,param)+posteriorMatrix(k-1,param));
		else{
			float prop=(trueValue-theatParameterScale(param, k-1))/(theatParameterScale(param, k)-theatParameterScale(param, k-1));
			float thisHeight=posteriorMatrix(k-1,param)+prop*(posteriorMatrix(k,param)-posteriorMatrix(k-1,param));
			float thisWidth=step*prop;
			area += thisWidth*(posteriorMatrix(k-1,param)+thisHeight);
			return area/(2*totalArea);
		}
	}
	return 1;
}
//---------------------------------------------------------------------------
float TStandardEstimation::getRMISEOfTrueParam(int param, float trueValue){
	float totalArea=getArea(posteriorMatrix.column(param), theta);
	float rmise=0;
	for(int k=1;k<=posteriorDensityPoints;++k){
		rmise+=step*posteriorMatrix(k,param)*(theatParameterScale(param, k)-trueValue)*(theatParameterScale(param, k)-trueValue);
	}
	return sqrt(rmise/totalArea);
}
//---------------------------------------------------------------------------
float TStandardEstimation::getREMODEOfTrueParam(int param, float trueValue){
	return fabs(getPositionofMode(param)-trueValue);
}
//---------------------------------------------------------------------------
float TStandardEstimation::getREMEANOfTrueParam(int param, float trueValue){
	return fabs(getPosteriorMean(param)-trueValue);
}
//---------------------------------------------------------------------------
float TStandardEstimation::getREMEDIANOfTrueParam(int param, float trueValue){
	return fabs(getPosteriorquantile(param, 0.5)-trueValue);
}
//---------------------------------------------------------------------------






