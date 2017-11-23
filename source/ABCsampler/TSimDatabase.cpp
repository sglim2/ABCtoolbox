//------------------------------------------------------------------------------


#include "TSimDatabase.h"
//------------------------------------------------------------------------------
TSimDatabase::TSimDatabase(int gotNum_cali_sims, TDataVector* gotDataPointer, TInputFileVector* gotinputFiles, ofstream* gotLogFile){
	num_sims=gotNum_cali_sims;
	myDataObject=gotDataPointer;
	inputFiles=gotinputFiles;
	logFile=gotLogFile;
	//prepare arrays to store the distances
	distances=new double[num_sims];
	//prepare arrays to store data and priors
	data=new double*[num_sims];
	priors=new double*[num_sims];
	for(int i=0; i<num_sims; ++i){
		data[i]=new double[myDataObject->numObsData];
		priors[i]=new double[inputFiles->priors->numSimplePrior];
	}
	num_sims_in_DB=0;
	thresholdSet=false;
	distancesCalculated=false;
	matrixWithSimsBelowThresholdFilled=false;
}
//------------------------------------------------------------------------------
void TSimDatabase::empty(int newsize){
	if(newsize!=num_sims){
		for(int i=0; i<num_sims; ++i){
			delete[] data[i];
			delete[] priors[i];
		}
		delete[] data;
		delete[] priors;
		num_sims=newsize;
		data=new double*[num_sims];
		priors=new double*[num_sims];
		for(int i=0; i<num_sims; ++i){
			data[i]=new double[myDataObject->numObsData];
			priors[i]=new double[inputFiles->priors->numSimplePrior];
		}
		delete[] distances;
		distances=new double[num_sims];
		if(matrixWithSimsBelowThresholdFilled) delete[] simsBelowThreshold;

	}
	num_sims_in_DB=0;
	thresholdSet=false;
	matrixWithSimsBelowThresholdFilled=false;
	distancesCalculated=false;
}
//------------------------------------------------------------------------------
TSimDatabase::TSimDatabase(my_string fileName, int gotNum_cali_sims, TDataVector* gotDataPointer, TInputFileVector* gotinputFiles, ofstream* gotLogFile){
	num_sims=gotNum_cali_sims;
	myDataObject=gotDataPointer;
	inputFiles=gotinputFiles;
	logFile=gotLogFile;
	*logFile << "Reading calibrationfile '" << fileName << "' ......";
	ifstream is (fileName.c_str()); // opening the file for reading
	if(!is)	throw TException("Calibration file '" + fileName + "' could not be opened!", _FATAL_ERROR);

	my_string buf;
	my_string my_value;
	strstream isLine;  // create the new stream

	//get number of simulations to read from parameter...
	distances=new double[num_sims];
	distances_normalized=new double[num_sims];

	//prepare arrays to store data and priors
	data_mean=new double[myDataObject->numObsData];
	data_var=new double[myDataObject->numObsData];
	data=new double*[num_sims];
	priors=new double*[num_sims];
	for(int i=0; i<num_sims; ++i){
		data[i]=new double[myDataObject->numObsData];
		priors[i]=new double[inputFiles->priors->numSimplePrior];
	}

	//read header line with names
	buf.read_line(is);

	//make a vector of all names
	vector<my_string> caliColumnNameVector;
	isLine.clear();
	isLine << buf;
	while(isLine){
		my_value.read_to_delim(isLine);
		if(my_value=="") break;
		caliColumnNameVector.push_back(my_value);
	}

	//create a vector of pointers pointing to the cali arrays storing the values for each cali File column
	int numCaliColumns=(int)caliColumnNameVector.size();
	double* caliColumnsOneLine= new double[numCaliColumns];
	int* caliColumnIndexPrior = new int[inputFiles->priors->numSimplePrior];
	int* caliColumnIndexData= new int[myDataObject->numObsData];
	//for the priors
	for(int k=0; k<inputFiles->priors->numSimplePrior;++k){
		for(int i=0; i<(numCaliColumns+1); ++i){
			if(i==numCaliColumns)
				throw TException("The parameter column '" + inputFiles->priors->simplePriors[k]->name + "' is missing in the calibration File!", _FATAL_ERROR);
			if(inputFiles->priors->simplePriors[k]->name==caliColumnNameVector[i]){
				caliColumnIndexPrior[k]=i;
				break;
			}
		}
	}
	//for the data
	myDataObject->matchColumns(caliColumnNameVector,numCaliColumns,&caliColumnIndexData);
	//read lines with data
	for(int i=0; i< num_sims ; ++i){
		if(is.good() && !is.eof()){
			buf.read_line(is);
			isLine.clear();
			isLine << buf;
			//read line into an array
			my_value.read_to_delim(isLine);
			caliColumnsOneLine[0]=my_value.toDouble();
			for(int k=1; k<numCaliColumns;++k){
				isLine >> ws;
				my_value.read_to_delim(isLine);
				caliColumnsOneLine[k]=my_value.toDouble();
			}
			//now write the values into the specific arrays
			for(int k=0; k<inputFiles->priors->numSimplePrior;++k){
				priors[i][k]=caliColumnsOneLine[caliColumnIndexPrior[k]];
			}
			for(int k=0; k<myDataObject->numObsData;++k){
				data[i][k] = caliColumnsOneLine[caliColumnIndexData[k]];
			}
		} else
			throw TException("File with simulations '" + fileName + "' is too short!", _FATAL_ERROR);
	}
	num_sims_in_DB=num_sims;
	*logFile << "done!" << endl;
}
//------------------------------------------------------------------------------
void TSimDatabase::addSimToDB(TDataVector* gotDataPointer, TInputFileVector* gotinputFiles, double distance){
	if(num_sims_in_DB==num_sims) throw TException("Database of simulations already full!", _FATAL_ERROR); //NEW ERROR
	for(int k=0; k<gotDataPointer->numObsData;++k){
		data[num_sims_in_DB][k]=*gotDataPointer->pointersToSimData[k];
	}
	for(int k=0; k<gotinputFiles->priors->numSimplePrior;++k){
		priors[num_sims_in_DB][k]=gotinputFiles->priors->simplePriors[k]->curValue;
	}
	if(distance>=0) distances[num_sims_in_DB]=distance;
	++num_sims_in_DB;
}
//------------------------------------------------------------------------------
void TSimDatabase::calculateMeanVariances(){
	data_mean=new double[myDataObject->numObsData];
	data_var=new double[myDataObject->numObsData];
	for(int k=0; k<myDataObject->numObsData;++k){
		data_mean[k]=0;
		data_var[k]=0;
		for(int i=0; i< num_sims ; ++i){
			data_mean[k]+=data[i][k];
		}
	}
	//calculate mean of simulated data
	for(int k=0; k<myDataObject->numObsData;++k){
		  data_mean[k]=data_mean[k]/num_sims;
	}
	//calculate variance of simulated data
	for(int k=0; k<myDataObject->numObsData;++k){
		data_var[k]=0;
		for(int i=0; i< num_sims ; ++i){
			data_var[k]=data_var[k]+(data_mean[k]-data[i][k])*(data_mean[k]-data[i][k]);
		}
		data_var[k]=data_var[k]/(num_sims-1);
		*myDataObject->pointersToVariances[k]=data_var[k];
	}

}
//---------------------------------------------------------------------------
void TSimDatabase::calculateMinMaxofParameters(double* & min, double* & max){
	min=new double[inputFiles->priors->numSimplePrior];
	max=new double[inputFiles->priors->numSimplePrior];
	for(int m=0; m<inputFiles->priors->numSimplePrior; ++m){
		min[m]=priors[0][m];
		max[m]=priors[0][m];
	}
	for(int i=1; i< num_sims ; ++i){
		for(int m=0; m<inputFiles->priors->numSimplePrior; ++m){
			if(min[m]>priors[i][m]) min[m]=priors[i][m];
			if(max[m]<priors[i][m]) max[m]=priors[i][m];
		}
	}
}
//---------------------------------------------------------------------------
void TSimDatabase::standardizeRetainedParameters(double* min, double* max){
	//standardizes the parameters between 0 and 1
	standardizedParameters=new double*[num_sims_below_threshold];
	for(int i=0; i<num_sims_below_threshold; ++i){
		standardizedParameters[i]=new double[inputFiles->priors->numSimplePrior];
	}
	for(int i=0; i< num_sims_below_threshold ; ++i){
		for(int m=0; m<inputFiles->priors->numSimplePrior; ++m){
			standardizedParameters[i][m]=(priors[simsBelowThreshold[i]][m]-min[m])/(max[m]-min[m]);
		}
	}
}
//---------------------------------------------------------------------------
void TSimDatabase::calculateWeightedSigmaOfRetainedParameters(SymmetricMatrix & Sigma, double* weights){
	if(!matrixWithSimsBelowThresholdFilled) throw TException("SimDatabase: Attempt to use the array of simulations below threshold before filling it!", _FATAL_ERROR);
	//get weighted mean of the different parameters...
	double means[inputFiles->priors->numSimplePrior];
	for(int m=0; m<inputFiles->priors->numSimplePrior; ++m) means[m]=0;
	for(int i=0; i<num_sims_below_threshold ; ++i){
		for(int m=0; m<inputFiles->priors->numSimplePrior; ++m){
			means[m]+=weights[i]*standardizedParameters[i][m];
		}
	}
	Sigma.ReSize(inputFiles->priors->numSimplePrior);
	Sigma=0;
	for(int i=0; i<num_sims_below_threshold ; ++i){
		for(int m=0; m<inputFiles->priors->numSimplePrior; ++m){
			for(int n=m; n<inputFiles->priors->numSimplePrior; ++n){
				Sigma.element(m, n)=Sigma.element(m, n)+weights[i]*(standardizedParameters[i][m]-means[m])*(standardizedParameters[i][n]-means[n]);
			}
		}
	}
	//Sigma=Sigma/(num_sims_below_threshold-1);
}
//---------------------------------------------------------------------------
double TSimDatabase::getPriorDensityOneRetainedSim(int sim){
	return inputFiles->getPriorDensity(priors[simsBelowThreshold[sim]]);
}
//---------------------------------------------------------------------------
void TSimDatabase::calculateDistances(){
	//calculate distances
	smallestDistCaliSim=0;
	distances[0]=myDataObject->calculateDistance(data[0]);
	for(int i=1; i< num_sims ; ++i){
		distances[i]=myDataObject->calculateDistance(data[i]);
		if(distances[i]<distances[smallestDistCaliSim]) smallestDistCaliSim=i;
	}
	distancesCalculated=true;
}
void TSimDatabase::calculateDistances(TLinearComb* linearComb){
	//do if linear combinations will be standardized for distance calculation...
	if(myDataObject->stdLinearCombForDist){
		//Save values to standardize
		double** cali_linearComb;
		double* cali_data_mean_linearCombination;
		double* cali_data_var_linearCombination;
		cali_linearComb=new double*[linearComb->numLinearComb];
		cali_data_mean_linearCombination=new double[linearComb->numLinearComb];
		cali_data_var_linearCombination=new double[linearComb->numLinearComb];

		for(int i=0; i<linearComb->numLinearComb;++i){
			cali_linearComb[i]=new double[num_sims];
			cali_data_var_linearCombination[i]=0;
			cali_data_mean_linearCombination[i]=0;
		}
		//calculate pca's
		for(int j=0; j< num_sims ; ++j){
			linearComb->calcSimDataPCA(data[j]);
			for(int i=0; i<linearComb->numLinearComb;++i){
				cali_linearComb[i][j]=linearComb->simPCA[i];
				cali_data_mean_linearCombination[i]+=linearComb->simPCA[i];
			}
		}
		//calculate mean and variance
		for(int i=0; i<linearComb->numLinearComb;++i){
			cali_data_mean_linearCombination[i]=cali_data_mean_linearCombination[i]/num_sims;
		}
		linearComb->linearCombVariances=new double[linearComb->numLinearComb];
		for(int i=0; i<linearComb->numLinearComb;++i){
			for(int j=0; j< num_sims ; ++j){
				cali_data_var_linearCombination[i]+=(cali_linearComb[i][j]-cali_data_mean_linearCombination[i])*(cali_linearComb[i][j]-cali_data_mean_linearCombination[i]);
			}
			linearComb->linearCombVariances[i]=cali_data_var_linearCombination[i]/(num_sims-1);
		}
	}

	//calculate distances
	smallestDistCaliSim=0;
	distances[0]=myDataObject->calculateDistance(data[0], linearComb);
	for(int j=0; j< num_sims ; ++j){
	   distances[j]=myDataObject->calculateDistance(data[j], linearComb);
	   if(distances[j]<distances[smallestDistCaliSim]) smallestDistCaliSim=j;
    }
	distancesCalculated=true;
}
//---------------------------------------------------------------------------
void TSimDatabase::fillNormalizedDistances(){
   //only needed to perform distance density estimates...
	delete distances_normalized;
	distances_normalized = new double[num_sims];
   //normalize distance
   distance_mean=0;
   for(int i=0; i<num_sims;++i) distance_mean+=distances[i];
   distance_mean=distance_mean/num_sims;
   distance_std=0;
   for(int i=0; i<num_sims;++i) distance_std+=(distance_mean-distances[i])*(distance_mean-distances[i]);
   distance_std=distance_std/(num_sims-1);
   distance_std=sqrt(distance_std);
   for(int i=0; i<num_sims;++i) distances_normalized[i]=(distances[i]-distance_mean)/distance_std;
}
//------------------------------------------------------------------------------
double TSimDatabase::getThreshold(double thresholdProportion){
	setThreshold((int)num_sims * thresholdProportion);
	return threshold;
}
void TSimDatabase::setThreshold(int numToRetain){
	num_sims_below_threshold=numToRetain;
	distances_ordered=new double*[num_sims];
	double* curDist=distances;
	double** curCaliDist=distances_ordered;
	for(int i=0; i< num_sims ;++curDist, ++curCaliDist, ++i){
		*curCaliDist=curDist;
	}
	quicksortP(distances_ordered, 0, num_sims-1);
	threshold=*distances_ordered[numToRetain];
	delete[] distances_ordered;
	thresholdSet=true;
}
//------------------------------------------------------------------------------
void TSimDatabase::writeToFile(my_string fileName){
	ofstream cf;
	cf.open(fileName.c_str());
	cf << inputFiles->priors->simplePriors[0]->name;
	for(int k=1; k<inputFiles->priors->numSimplePrior;++k){
		cf << "\t" << inputFiles->priors->simplePriors[k]->name;
	}
	myDataObject->writeHeader(cf);
	cf << endl;
	for(int i=0; i< num_sims ; ++i){
		cf << priors[i][0];
		for(int k=1; k<inputFiles->priors->numSimplePrior;++k){
			cf << "\t" << priors[i][k];
		}
		for(int k=0; k<myDataObject->numObsData;++k){
			cf << "\t" << data[i][k];
		}
		cf << endl;
	}
	cf.close();
}
//------------------------------------------------------------------------------
void TSimDatabase::writeDistanceFile(my_string fileName){
   *logFile << "distance File written...: " << fileName << endl;
	ofstream cf;
	cf.open(fileName.c_str());
	cf << num_sims << endl;
	for(int i=0; i< num_sims ; ++i){
		cf << distances[i] << endl;
	}
	cf.close();
}
//------------------------------------------------------------------------------
void TSimDatabase::fillArrayOfSimsBelowThreshold(){
   if(!thresholdSet) throw TException("SimDatabase: Attempt to fill array of simulations below threshold before the threshold was computed!", _FATAL_ERROR);
   simsBelowThreshold=new int[num_sims_below_threshold];
   int p=0;
   for(int i=0; i<num_sims ; ++i){
	  if(distances[i]<=threshold && p<num_sims_below_threshold){
		 simsBelowThreshold[p]=i;
		 ++p;
	  }
   }
   //fix num_sims_below_threshold as in the case of an iterative mcmc some sims may be replicate several times
   num_sims_below_threshold=p;
   matrixWithSimsBelowThresholdFilled=true;
}
//------------------------------------------------------------------------------
double TSimDatabase::getDistanceDensity(double dist){
//get the density at the given distance "dist"
	return epanechnikovDensityEstimate(distances_normalized, num_sims, (dist-distance_mean)/distance_std);
}
//------------------------------------------------------------------------------
void TSimDatabase::setMcmcRanges(float mcmc_range_proportion){
   fillArrayOfSimsBelowThreshold();
   for(int k=0; k<inputFiles->priors->numSimplePrior;++k)
	  setMcmcRangesOnePrior(k, mcmc_range_proportion);
}
void TSimDatabase::setMcmcRangesOnePrior(int prior, float mcmc_range_proportion){
	  double mean=0;
	  for(int i=0; i<num_sims_below_threshold ; ++i)
			mean+=priors[simsBelowThreshold[i]][prior];
	  mean=mean/num_sims_below_threshold;
	  inputFiles->priors->simplePriors[prior]->mcmcStep=0;
	  for(int i=0; i<num_sims_below_threshold; ++i)
		 inputFiles->priors->simplePriors[prior]->mcmcStep+=(priors[simsBelowThreshold[i]][prior]-mean)*(priors[simsBelowThreshold[i]][prior]-mean);
	  inputFiles->priors->simplePriors[prior]->mcmcStep=sqrt(inputFiles->priors->simplePriors[prior]->mcmcStep/(num_sims_below_threshold-1))*mcmc_range_proportion;
}
//------------------------------------------------------------------------------
//Functions to set the priors to a certain point in the parameter space.
void TSimDatabase::setPriorStartingConditionsFromBestSimulation(){
   for(int k=0; k < inputFiles->priors->numSimplePrior; ++k){
	  inputFiles->priors->simplePriors[k]->curValue=priors[smallestDistCaliSim][k];
	  inputFiles->priors->simplePriors[k]->oldValue=priors[smallestDistCaliSim][k];
   }
   *logFile << "Starting point is set from the best Cali-Simulation" << endl;
   // calc all combined parameters
   inputFiles->priors->updateCombinedParameters();
}
void TSimDatabase::setPriorStartingConditionsAtRandom(){
   int i=floor(UniformRandom(0, 1)*(num_sims_below_threshold-0.000001));
   for(int k=0; k < inputFiles->priors->numSimplePrior; ++k){
	  inputFiles->priors->simplePriors[k]->curValue=priors[simsBelowThreshold[i]][k];
	  inputFiles->priors->simplePriors[k]->oldValue=priors[simsBelowThreshold[i]][k];
   }
   *logFile << "Starting point is set randomly at simulation " << i << endl;
   // calc all combined parameters
   inputFiles->priors->updateCombinedParameters();
}

