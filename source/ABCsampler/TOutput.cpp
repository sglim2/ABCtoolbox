//---------------------------------------------------------------------------

#pragma hdrstop

#include "TOutput.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
TOutput::TOutput(int gotSampling, my_string outputFilename, bool gotWriteDistance, TDataVector* Data, TInputFileVector* InputFiles){
	data=Data;
	inputFiles=InputFiles;
	sampling=gotSampling;
	writeDistance=gotWriteDistance;
	filename=outputFilename+"_sampling"+sampling+".txt";
	myFile.open(filename.c_str());
}
TOutputOneObs::TOutputOneObs(int gotSampling, int gotObsFileNumber, my_string outputFilename, bool gotWriteDistance, TDataVector* Data, TInputFileVector* InputFiles){
	data=Data;
	inputFiles=InputFiles;
	sampling=gotSampling;
	obsFileNumber=gotObsFileNumber;
	writeDistance=gotWriteDistance;
	filename=outputFilename+"_Obs"+obsFileNumber+"_sampling"+sampling+".txt";
	myFile.open(filename.c_str());
}
//------------------------------------------------------------------------------
void TOutput::writeHeader(){
   myFile << "Sim";
   inputFiles->priors->writeHeader(myFile);
   data->writeHeader(myFile);
   if(writeDistance) myFile << "\t" << "Distance";
   myFile << endl;
}
void TOutput::writeHeader(TLinearComb* pcaObject){
   myFile << "Sim";
   inputFiles->priors->writeHeader(myFile);
   data->writeHeader(myFile);
   pcaObject->writeHeader(myFile);
   if(writeDistance) myFile << "\t" << "Distance";
   myFile << endl;
}
void TOutputOneObs::writeHeader(){
   myFile << "Sim";
   inputFiles->priors->writeHeader(myFile);
   data->writeHeaderOnlyOneDataObject(obsFileNumber, myFile);
   if(writeDistance) myFile << "\t" << "Distance";
   myFile << endl;
}
//------------------------------------------------------------------------------
void TOutput::writeSimulations(int run, double distance){
   if(run%sampling==0){
	  myFile << run;
	  inputFiles->priors->writeParameters(myFile);
	  data->writeData(myFile);
	  if(writeDistance) myFile << "\t" << distance;
	  myFile << endl;
   }
}
void TOutput::writeSimulations(int run, TLinearComb* pcaObject, double distance){
   if(run%sampling==0){
	  myFile << run;
	  inputFiles->priors->writeParameters(myFile);
	  data->writeData(myFile);
	  pcaObject->writeSimPCA(myFile);
	  if(writeDistance) myFile << "\t" << distance;
	  myFile << endl;
   }
}
void TOutputOneObs::writeSimulations(int run, double distance){
   if(run%sampling==0){
	  myFile << run;
	  inputFiles->priors->writeParameters(myFile);
	  data->writeDataOnlyOneDataObject(obsFileNumber, myFile);
	  if(writeDistance) myFile << "\t" << distance;
	  myFile << endl;
   }
}
//------------------------------------------------------------------------------
TOutputVector::TOutputVector(my_string mcmcsampling, my_string outName, TDataVector* Data, TInputFileVector* InputFile, bool gotWriteDistance, bool gotSeparateOutputFiles){
	separateOutputFiles=gotSeparateOutputFiles;
	if(separateOutputFiles){
	   while(!mcmcsampling.empty()){
		  int sampling=mcmcsampling.extract_sub_str(*";").toInt();
		  for(int i=0; i<Data->numDataObjects;++i){
			 TOutput* newObject = new TOutputOneObs(sampling, i, outName, gotWriteDistance, Data, InputFile);
			 //this is done because the Vector wants a TOutput object and produces an error if the TOutputOneObs is directly created on the next line
			 vecOutputFiles.push_back(newObject);
			 //vecOutputFiles.push_back(new TOutput(mcmcsampling.extract_sub_str(*";").toInt(), outName, gotWriteDistance, Data, InputFile));
		  }
		  mcmcsampling.remove(0,1);
	   }
	} else {
	   while(!mcmcsampling.empty()){
	      vecOutputFiles.push_back(new TOutput(mcmcsampling.extract_sub_str(*";").toInt(), outName, gotWriteDistance, Data, InputFile));
		  mcmcsampling.remove(0,1);
	   }
    }
	endOutput=vecOutputFiles.end();
}
//------------------------------------------------------------------------------
void TOutputVector::writeHeader(){
	for(unsigned int i=0; i< vecOutputFiles.size();++i){
		vecOutputFiles[i]->writeHeader();
	}
}
void TOutputVector::writeHeader(TLinearComb* pcaObject){
	for(unsigned int i=0; i< vecOutputFiles.size();++i){
		vecOutputFiles[i]->writeHeader(pcaObject);
	}
}
//------------------------------------------------------------------------------
void TOutputVector::writeSimulations(int run){
   writeSimulations(run, 0.0);
}
void TOutputVector::writeSimulations(int run, double distance){
	//curOutput=vecOutputFiles.begin();
	//for(; curOutput!=endOutput; ++curOutput){
	for(unsigned int i=0; i< vecOutputFiles.size();++i){
	   vecOutputFiles[i]->writeSimulations(run, distance);
	}
}
void TOutputVector::writeSimulations(int run, TLinearComb* pcaObject){
   writeSimulations(run, pcaObject, 0.0);
}
void TOutputVector::writeSimulations(int run, TLinearComb* pcaObject, double distance){
	//curOutput=vecOutputFiles.begin();
	//for(; curOutput!=endOutput; ++curOutput){
	for(unsigned int i=0; i< vecOutputFiles.size();++i){
	   vecOutputFiles[i]->writeSimulations(run, pcaObject, distance);
	}
}
