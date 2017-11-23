//---------------------------------------------------------------------------
#include "TObsData.h"
//---------------------------------------------------------------------------
TObsDataSet::TObsDataSet(int gotNumStats){
   internalCounter=0;
   numStats=gotNumStats;
   myObsValues=new float[numStats];
   myStandardizedObsValues=NULL;
};
void TObsDataSet::add(float value){
   if(internalCounter==numStats) throw TException("Too many values added to an observed data set!", _FATAL_ERROR);
   myObsValues[internalCounter]=value;
   ++internalCounter;
};
void TObsDataSet::setNumTrueParams(int value){
	internalTrueParamCounter=0;
	numTrueParams=value;
	myTrueParams=new float[numTrueParams];
};
void TObsDataSet::addTrueParams(float value){
   if(internalTrueParamCounter==numTrueParams) throw TException("Too many true parameter values added to an observed data set!", _FATAL_ERROR);
   myTrueParams[internalTrueParamCounter]=value;
   ++internalTrueParamCounter;
};
void TObsDataSet::standardize(float* means, float* sds){
   myStandardizedObsValues=new float[numStats];
   for(int i=0; i<numStats;++i){
	  myStandardizedObsValues[i]=(myObsValues[i]-means[i])/sds[i];
   }
}
void TObsDataSet::getStandardizedObsValuesIntoColumnVector(ColumnVector& colVector){
   colVector.resize(numStats);
   for(int i=0; i<numStats; ++i){
	   colVector(i+1)= myStandardizedObsValues[i];
   }
}
void TObsDataSet::getObsValuesIntoColumnVector(ColumnVector& colVector){
   colVector.resize(numStats);
   for(int i=0; i<numStats; ++i){
	   colVector(i+1)= myObsValues[i];
   }
}

//---------------------------------------------------------------------------
TObsData::TObsData(my_string obsfilename){
   //this constructer reads the files
   obsFileName=obsfilename;
   hasBeenStandardized=false;
   readObsFile();
}
//---------------------------------------------------------------------------
void TObsData::readObsFile(){
   cout << "- Reading observed data file '" << obsFileName << "' ...";
   ifstream is (obsFileName.c_str()); // opening the file for reading
   if(!is) throw TException("File with observed statistics '" + obsFileName + "' could not be opened!", _FATAL_ERROR);
   my_string buf;
   strstream isLine;  // create the new stream
   buf.read_line(is); // read the line
   isLine << buf;     //allocate the line to a new stream for easy reading
   //read header line with name
   while(isLine){
	  buf.read_to_delim(isLine);
	  if(buf.empty()) break;
	  obsNameVector.push_back(buf);
   }
   //read values into an array
   endObsNameVector=obsNameVector.end();
   numStats=obsNameVector.size();
   obsValuesMatrix=Matrix(numStats,1); //prepare for later..
   //every line is a new dataset
   int dataSetCounter=0;
   while(is.good() && !is.eof()){
	  isLine.clear();
	  buf.read_line(is); // read the line
	  if(!buf.empty()){
		 isLine << buf;     //allocate the line to a new stream for easy reading
		 obsDataSets.push_back(new TObsDataSet(numStats));
		 int i=0;
		 while(isLine){
			buf.read_to_delim(isLine);
			if(buf.empty( )) break;
			if(i==numStats) throw TException("More values than names in the .obs-file '" + obsFileName + "'!", _FATAL_ERROR);
			obsDataSets.at(dataSetCounter)->add(buf.toDouble());
			++i;
		 }
		 if(i<(numStats-1)) throw TException("More names than values in the .obs-file '" + obsFileName + "'!", _FATAL_ERROR);
		 ++dataSetCounter;
	  }
   }

   curObsNameVector=obsNameVector.begin();
   numObsDataSets=obsDataSets.size();
   is.close();
   cout << " done!" << endl;
   cout << "   -> " << numObsDataSets << " data sets with " << numStats << " statistics each." << endl;
}
//---------------------------------------------------------------------------
void TObsData::readTrueFile(my_string truefilename){
   cout << "- Reading file with true parameters'" << truefilename << "' ...";
   ifstream is (truefilename.c_str()); // opening the file for reading
   if(!is) throw TException("File with true parameters '" + truefilename + "' could not be opened!", _FATAL_ERROR);
   my_string buf;
   strstream isLine;  // create the new stream
   buf.read_line(is); // read the line
   isLine << buf;     //allocate the line to a new stream for easy reading
   //read header line with name
   while(isLine){
	  buf.read_to_delim(isLine);
	  if(buf.empty()) break;
	  trueNameVector.push_back(buf);
   }
   //read values into an array
   numTrueParams=trueNameVector.size();
   //every line is a new dataset and corresponds to the observed data sets...
   int dataSetCounter=0;
   while(is.good() && !is.eof()){
	   buf.read_line(is); // read the line
	   if(!buf.empty()){
		   if(dataSetCounter==numObsDataSets) throw TException("Too many lines in the file with true parameters '" + truefilename + "', more than observed data set!", _FATAL_ERROR);
		   isLine.clear();
		   isLine << buf;     //allocate the line to a new stream for easy reading
		   obsDataSets.at(dataSetCounter)->setNumTrueParams(numTrueParams);
		   int i=0;
		   while(isLine){
			   buf.read_to_delim(isLine);
			   if(buf.empty( )) break;
			   if(i==numTrueParams) throw TException("More values than names in the file with true parameters'" + truefilename + "'!", _FATAL_ERROR);
			   obsDataSets.at(dataSetCounter)->addTrueParams(buf.toDouble());
			   ++i;
		   }
		   if(i<(numTrueParams-1)) throw TException("More names than values in the file with true parameters'" + truefilename + "'!", _FATAL_ERROR);
		   ++dataSetCounter;
	   }
   }
   if(dataSetCounter<numObsDataSets) throw TException("Less lines in the file with true parameters'" + truefilename + "' than observed data set!", _FATAL_ERROR);
   is.close();
   cout << " done!" << endl;
   cout << "   -> " << numObsDataSets << " sets with " << numTrueParams << " true parameters each." << endl;
}
//---------------------------------------------------------------------------
int TObsData::getStatColNumberFromName(my_string name){
   int i=0;
   curObsNameVector=obsNameVector.begin();
   for(;curObsNameVector!=endObsNameVector;++curObsNameVector, ++i){
	  if(*curObsNameVector==name) return i;
   }
   return -1;
}
//---------------------------------------------------------------------------
void TObsData::fillObsDataMatrix(int dataSet){
   if(dataSet<numObsDataSets) obsValuesMatrix << obsDataSets[dataSet]->myStandardizedObsValues;
   else {
	  throw TException("There are less data sets in the obsfile than requested!", _FATAL_ERROR);
   }
}
//---------------------------------------------------------------------------
void TObsData::standardizeObservedValues(float* means, float*sds){
	for(int i=0; i<numObsDataSets;++i){
		obsDataSets[i]->standardize(means, sds);
	  }
	  hasBeenStandardized=true;
}
//---------------------------------------------------------------------------
float* TObsData::getObservedValues(int num){
	if(hasBeenStandardized) return obsDataSets[num]->myStandardizedObsValues;
	else return obsDataSets[num]->myObsValues;
}
//---------------------------------------------------------------------------
void TObsData::getObsValuesIntoColumnVector(int num, ColumnVector& colVector){
   if(hasBeenStandardized) obsDataSets[num]->getStandardizedObsValuesIntoColumnVector(colVector);
   else obsDataSets[num]->getObsValuesIntoColumnVector(colVector);
}
//---------------------------------------------------------------------------
void TObsData::writeTrueParamsheader(ofstream* file){
	for(int i=0; i<numTrueParams; ++i) *file << "\t" << trueNameVector[i];
}
//---------------------------------------------------------------------------
float TObsData::getTrueParam(int num, int trueparam){
	return obsDataSets[num]->myTrueParams[trueparam];
}




