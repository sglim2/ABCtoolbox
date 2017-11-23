//---------------------------------------------------------------------------

#pragma hdrstop

#include "TSimData.h"
#pragma package(smart_init)
//---------------------------------------------------------------------------

TSimData::TSimData(my_string simfilename, my_string params, int maxSimsToRead, TObsData* obsData){
   //this constructer reads the file
   simFileName=simfilename;
   paramString=params;
   pointerToObsData=obsData;
   statMeansCalculated=false;
   statSDsCalculated=false;
   statsAreStandardized=false;
   statMinMaxCalculated=false;
   paramsStandardized=false;
   //read the file
   readSimFile(maxSimsToRead);
   numUsedSims=0;
}
//---------------------------------------------------------------------------
void TSimData::readSimFile(int maxSimsToRead){
   cout << "- Reading simulation file '" << simFileName << "' ...";
   //open file
   ifstream is (simFileName.c_str()); // opening the file for reading
   if(!is) throw TException("File with simulations '" + simFileName + "' could not be opened!", _FATAL_ERROR);
   my_string buf;
   strstream isLine;  // create the new stream
   buf.read_line(is); // read the line
   isLine << buf;     //allocate the line to a new stream for easy reading
   vector<my_string> simNameVector;
   //read header line with name
   while(isLine){
	  buf.read_to_delim(isLine);
	  if(buf=="") break;
	  simNameVector.push_back(buf);
   }

   //go through the names and find corresponding columns for the stats
   int* matchSimToObsStats;
   bool* simColIsStat;
   simColIsStat=new bool[simNameVector.size()];
   matchSimToObsStats=new int[simNameVector.size()];
   for(unsigned int j=0; j<simNameVector.size(); ++j) simColIsStat[j]=false;
   int colnum;
   int i;
   int numStatsMatched=0;
   for(i=0; i<simNameVector.size();++i){
	  colnum=pointerToObsData->getStatColNumberFromName(simNameVector[i]);
	  if(colnum>=0){
		 simColIsStat[i]=true;
		 matchSimToObsStats[i]=colnum;
		 ++numStatsMatched;
	  }
   }
   if(numStatsMatched!=pointerToObsData->numStats) throw TException("Some statistic(s) is/are missing in the simulation file!", _FATAL_ERROR);;

   //now match parameters
   numParams=0;
   bool* simColIsParam;
   simColIsParam=new bool[simNameVector.size()];
   for(unsigned int j=0; j<simNameVector.size(); ++j) simColIsParam[j]=false;
   //read parameter string... -> get which params to estimate!
   my_string temp, Temp;
   i=0;
   while(!paramString.empty()){
	   temp=paramString.extract_sub_str(*",");
	   paramString.remove(0,1);
	   //if number, put back.

	   if(!temp.contains(*"-") && temp.toInt()>0){
		   simColIsParam[temp.toInt()-1]=true;
		   ++numParams;
	   }
	   //If sequence fill sequence...
	   else {
		   Temp=temp.extract_sub_str(*"-");
		   temp.remove(0,1);
		   if(Temp.toInt()>0 && temp.toInt()>Temp.toInt()){
			   for(unsigned int j=Temp.toInt()-1; j<temp.toInt(); ++j) simColIsParam[j]=true;
			   numParams=numParams+(temp.toInt()-Temp.toInt())+1;
		   }
		   else throw TException("Problem reading the parameter string!", _FATAL_ERROR);
	   }
   }

   //check if some params are also stats....
   for(unsigned int j=0; j<simNameVector.size(); ++j){
	   if(simColIsStat[j] && simColIsParam[j]) throw TException("Column '"+simNameVector[j]+"' is requested as parameter and statistic!", _FATAL_ERROR);
   }

   //store parameter names
   paramNames=new my_string[numParams];
   i=0;
   for(unsigned int j=0; j<simNameVector.size(); ++j){
	   if(simColIsParam[j]){
		   paramNames[i]=simNameVector[j];
		   ++i;
	   }
   }

   //prepare arrays to store sims and stats
   simParams=new float*[maxSimsToRead];
   simStats=new float*[maxSimsToRead];
   for(int j=0; j<maxSimsToRead; ++j){
	  simParams[j]=new float[numParams];
	  simStats[j]=new float[pointerToObsData->numStats];
   }

   //now read params and stats into big arrays
   numReadSims=0;
   while(is.good() && !is.eof() && numReadSims<maxSimsToRead){
	  isLine.clear();
	  buf.read_line(is); // read the line
	  isLine << buf;     //allocate the line to a new stream for easy reading
	  if(!buf.empty()){
		  i=0;
		  int p=0;
		  while(isLine){
			  if(i==simNameVector.size()) throw TException("Too many columns on line " + (my_string) (numReadSims+1) + " of the simulation file '" + simFileName + "'!", _FATAL_ERROR);
			  buf.read_to_delim(isLine);
			  if(buf.empty()) break;
			  else {
				  if(simColIsParam[i]){
					  simParams[numReadSims][p]=buf.toDouble();
					  ++p;
				  }
				  else {
					  if(simColIsStat[i]){
						  simStats[numReadSims][matchSimToObsStats[i]]=buf.toDouble();
					  }
				  }
				  ++i;
			  }
		  }
		  if(i<simNameVector.size()) throw TException("Too few columns on line " + (my_string) (numReadSims+1) + " of the simulation file '" + simFileName + "'!", _FATAL_ERROR);
		  ++numReadSims;
	  }
   }

   //if less than maxSimsToRead sims have been read: delete some arrays!
   if(numReadSims<maxSimsToRead){
      for(int j=maxSimsToRead;j>numReadSims;--j){
         delete simStats[j-1];
		 delete simParams[j-1];
	  }
   }

   //delete variables
   delete[] matchSimToObsStats;
   delete[] simColIsParam;

   //done!
   is.close();
   cout << " done!" << endl;
   cout << "   -> " << numReadSims << " simulations with " << numParams << " parameters and " << pointerToObsData->numStats << " statistics each." << endl;
}
//---------------------------------------------------------------------------
int TSimData::getParamNumberFromName(my_string name){
	for(int i=0; i<numParams;++i){
		if(paramNames[i]==name) return i;
	}
	return -1;
}
//---------------------------------------------------------------------------
void TSimData::calculateStatisticsMeans(){
	//prepare Mean array
	statMeans=new float[pointerToObsData->numStats];
	for(int j=0; j<pointerToObsData->numStats; ++j){
		statMeans[j]=0;
	}
   //calculate means
   for(int i=0; i<numReadSims;++i){
	  for(int j=0; j<pointerToObsData->numStats; ++j){
	     statMeans[j]+=simStats[i][j];
	  }
   }
   for(int j=0; j<pointerToObsData->numStats; ++j){
	  statMeans[j]=statMeans[j]/numReadSims;
	  //statMeans[j]=0;
   }
   statMeansCalculated=true;
}
//---------------------------------------------------------------------------
void TSimData::calculateStatisticsSDs(){
	if(!statMeansCalculated) calculateStatisticsMeans();
	//prepare SD array
	statSDs=new float[pointerToObsData->numStats];
	for(int j=0; j<pointerToObsData->numStats; ++j){
		statSDs[j]=0;
	}
	//calculate SDs
	for(int i=0; i<numReadSims;++i){
	  for(int j=0; j<pointerToObsData->numStats; ++j){
		 statSDs[j]+=(simStats[i][j]-statMeans[j])*(simStats[i][j]-statMeans[j]);
	  }
	}
	for(int j=0; j<pointerToObsData->numStats; ++j){
	  statSDs[j]=sqrt(statSDs[j]/(numReadSims-1));
	  //test if its is nan or inf
	  if(statSDs[j] != statSDs[j] || statSDs[j] == std::numeric_limits<float>::infinity())
		  throw TException("Problem calculating the standard deviation of the simulated statistic '"+pointerToObsData->obsNameVector[j]+"'!", _FATAL_ERROR);
	}
	statSDsCalculated=true;
}
//---------------------------------------------------------------------------
void TSimData::standardizeStatistics(){
	if(!statSDsCalculated) calculateStatisticsSDs();
	//standardize the stats
	for(int i=0; i<numReadSims;++i){
		for(int j=0; j<pointerToObsData->numStats; ++j){
			simStats[i][j]=(simStats[i][j]-statMeans[j])/statSDs[j];
		}
	}
	statsAreStandardized=true;
}
//---------------------------------------------------------------------------
void TSimData::standardizeParameters(){
   paramMaxima=new float[numParams];
   paramMinima=new float[numParams];
   for(int j=0; j<numParams; ++j){
      paramMaxima[j]=simParams[0][j];
      paramMinima[j]=simParams[0][j];
   }
   //get max value
   for(int i=1; i<numReadSims;++i){
	  for(int j=0; j<numParams; ++j){
		 if(simParams[i][j]>paramMaxima[j]) paramMaxima[j]=simParams[i][j];
		 if(simParams[i][j]<paramMinima[j]) paramMinima[j]=simParams[i][j];
		 if(simParams[i][j]==0){
			 float a=5;
		 }
	  }
   }
   //now standardize the params by substracting the minimum and dividing them by the maximum
   for(int i=0; i<numReadSims;++i){
	  for(int j=0; j<numParams; ++j){
		 simParams[i][j]=(simParams[i][j]-paramMinima[j])/(paramMaxima[j]-paramMinima[j]);
	  }
   }
   paramsStandardized=true;
}
//---------------------------------------------------------------------------
my_string TSimData::checkIfstatsArePolymorphic(){
	for(int i=1; i<=pointerToObsData->numStats;++i){
		double min=statMatrix.row(i).minimum();
		double max=statMatrix.row(i).maximum();
		if(min==max)
			return pointerToObsData->obsNameVector[i-1];
	}
	return "";
}
//---------------------------------------------------------------------------
void TSimData::getMinMaxofStats(){
	//compute min and max
	statMaxima=new float[pointerToObsData->numStats];
	statMinima=new float[pointerToObsData->numStats];
	for(int j=0; j<pointerToObsData->numStats; ++j){
		statMinima[j]=simStats[0][j];
		statMaxima[j]=simStats[0][j];
	}
	for(int i=1; i<numReadSims;++i){
		for(int j=0; j<pointerToObsData->numStats; ++j){
			if(statMinima[j]>simStats[i][j]) statMinima[j]=simStats[i][j];
			if(statMaxima[j]<simStats[i][j]) statMaxima[j]=simStats[i][j];
		}
	}
	statMinMaxCalculated=true;
}
//---------------------------------------------------------------------------
void TSimData::checkIfstatsAreWithinRange(float* obsValueArray){
	if(!statMinMaxCalculated) getMinMaxofStats();
	for(int i=0; i<pointerToObsData->numStats;++i){
		if(obsValueArray[i]>statMaxima[i] || obsValueArray[i]<statMinima[i])
			cout << "\nWARNING: the observed value of the statistic '"+pointerToObsData->obsNameVector[i]+"' is outside the simulated range of this statistic!";
	}
}

//---------------------------------------------------------------------------
void TSimData::calculateDistances(float* obsValueArray, float* distances){
   cout << "   - calculating distances ...";
   checkIfstatsAreWithinRange(obsValueArray);
   float distance;
   for(int i=0; i<numReadSims;++i){
	  distance=0;
	  for(int j=0; j<pointerToObsData->numStats; ++j){
		 distance+=(simStats[i][j]-obsValueArray[j])*(simStats[i][j]-obsValueArray[j]);
	  }
	  //distances[i]=sqrt(distance); --> does not matter since only the order is important....
	  distances[i]=distance;
   }
   cout << " done!" << endl;
}
//---------------------------------------------------------------------------
int TSimData::fillStatAndParamMatricesKeepAll(){
	cout << "   - use all " << numReadSims << " simulations for the estimation (no rejection)...";
	//create matrices
	paramMatrix=Matrix(numReadSims, numParams+1);
	statMatrix=Matrix(pointerToObsData->numStats, numReadSims);
	//fill matrices
	for(int i=1; i<=numReadSims;++i){
		paramMatrix.submatrix(i,i,2,numParams+1) << simParams[i-1];
		statMatrix.column(i) << simStats[i-1];
	}
   paramMatrix.column(1) = 1;
   numUsedSims=numReadSims;
   cout << " done!" << endl;
   return numUsedSims;
}
//---------------------------------------------------------------------------
int TSimData::fillStatAndParamMatricesFromNumRetained(float* distances, int numToRetain){
   if(numToRetain<numReadSims){
	   cout << "   - retain the best " << numToRetain << " simulations ...";
	  //sort the distances
	  float** distances_ordered=new float*[numReadSims];
	  float* curDist=distances;
	  float** curCaliDist=distances_ordered;
	  for(int i=0; i< numReadSims ;++curDist, ++curCaliDist,++i){
		*curCaliDist=curDist;
	  }
	  quicksortP(distances_ordered, 0, numReadSims-1);
	  //get threshold
	  threshold=*distances_ordered[numToRetain];
	  delete[] distances_ordered;
	  //create matrices

	  paramMatrix=Matrix(numToRetain, numParams+1);
	  statMatrix=Matrix(pointerToObsData->numStats, numToRetain);
	  //fill matrices
	  int j=1;
	  for(int i=0; i<numReadSims;++i){
		  if(distances[i]<=threshold && j<=numToRetain){
			 paramMatrix.submatrix(j,j,2,numParams+1) << simParams[i];
			 statMatrix.column(j) << simStats[i];
			 ++j;
		 }
	  }

	  cout << " done!" << endl;
	  paramMatrix.column(1) = 1;
	  numUsedSims=numToRetain;
	  return numUsedSims;
   } else {
	   cout << " given the numretained specified, all simulations are retained!" << endl;
	   return fillStatAndParamMatricesKeepAll();
   }
}
//---------------------------------------------------------------------------
int TSimData::fillStatAndParamMatricesFromThreshold(float* distances, float gotThreshold, int numToRetain){
   if(numToRetain<0) cout << "   - retain the best simulations with tolerance " << gotThreshold << " ...";
   else cout << "   - retain the best simulations with tolerance " << gotThreshold << " and numRetained " << numToRetain << " ...";
   //save threshold
   threshold=gotThreshold;

   //get number of sims to retain
   if(numToRetain<0){
	   for(int i=0; i<numReadSims;++i) if(distances[i]<=threshold) ++numToRetain;
	   numUsedSims=numToRetain;
   } else {
	   numUsedSims=0;
	   for(int i=0; i<numReadSims;++i){
	       if(distances[i]<=threshold && numUsedSims<numToRetain) ++numUsedSims;
	   }
   }
   //create matrices
   paramMatrix=Matrix(numUsedSims, numParams+1);
   statMatrix=Matrix(pointerToObsData->numStats, numUsedSims);
   //fill matrices
   int j=1;
   for(int i=0; i<numReadSims;++i){
	   if(distances[i]<=threshold && j<=numToRetain){
		   paramMatrix.submatrix(j,j,2,numParams+1) << simParams[i];
		   statMatrix.column(j) << simStats[i];
		   ++j;
	   }
   }
   paramMatrix.column(1) = 1;
   cout << " done, " << numUsedSims << " simulations retained!" << endl;
   return numUsedSims;
}
//---------------------------------------------------------------------------
void TSimData::writeRetained(my_string outputPrefix, float* distances, int dataSet){
   cout << "   - writing retained ...";
   my_string filename=outputPrefix+"BestSimsParamStats";
   if(dataSet>=0) filename+="_Obs"+(my_string)dataSet;
   filename+=".txt";
   ofstream output;
   output.open(filename.c_str());
   //write Header
   output << "Sim_num\tDist";
   for(int i=0;i<numParams;++i) output << "\t" << paramNames[i];
   for(int i=0;i<pointerToObsData->numStats;++i) output << "\t" << pointerToObsData->obsNameVector[i];
   output << endl;
   //write stats -> destandardize if necessary / param -> destandardize!!
   int j=1;
   for(int i=0; i<numReadSims;++i){
	  if(distances[i]<=threshold && j<=numUsedSims){
		 output << (i+1) << "\t" << distances[i];
		 //for(int k=0;k<numParams;++k) output << "\t" << simParams[i][k]*paramMaxima[k];
		 for(int k=0;k<numParams;++k) output << "\t" << simParams[i][k]*(paramMaxima[k]-paramMinima[k])+paramMinima[k];
		 if(statsAreStandardized) for(int k=0;k<pointerToObsData->numStats;++k) output << "\t" << (simStats[i][k]*statSDs[k])+statMeans[k];
		 else for(int k=0;k<pointerToObsData->numStats;++k) output << "\t" << simStats[i][k];
		 output << endl;
		 ++j;
	  }
   }
   output.close();
   cout << " done!" << endl;
}
//---------------------------------------------------------------------------
//functions to sort an array of pointers
float TSimData::partitionP(float** a, int left, int right, int pivotIndex) {
	int storeIndex=left;
	//swap
	float* temp=a[right];
	a[right]=a[pivotIndex];
	a[pivotIndex]=temp;
	for (int i=left; i < right; i++){
		if(*a[i] <= *a[right]){
			//swap
			temp=a[storeIndex];
			a[storeIndex]=a[i];
			a[i]=temp;
			++storeIndex;
		}
	}
	//Pivot back to final place
	temp=a[right];
	a[right]=a[storeIndex];
	a[storeIndex]=temp;
	return storeIndex;
}

void TSimData::quicksortP(float** a, int left, int right){
	if(right > left){
		int pivotIndex=(left+right)/2;
		pivotIndex=partitionP(a, left, right, pivotIndex);
		quicksortP(a, left, pivotIndex-1);
		quicksortP(a, pivotIndex+1, right);
	}
}
//---------------------------------------------------------------------------

