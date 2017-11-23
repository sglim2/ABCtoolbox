//---------------------------------------------------------------------------


#pragma hdrstop
#include "TLinearCombBoxCox.h"
#include "TException.h"
//---------------------------------------------------------------------------

#pragma package(smart_init)

TLinearCombBoxCox::TLinearCombBoxCox(my_string gotLinearCombFileName, vector<my_string> NamesInInputFiles){
	LinearCombFileName=gotLinearCombFileName;
	namesInInputFiles=NamesInInputFiles;
	readLinearCombFile();
}
//------------------------------------------------------------------------------
void TLinearCombBoxCox::readLinearCombFile(){
	numLinearComb=0;
	my_string name, buf, filename;
	strstream isLine;                   // create the new stream
	vector<double> tempPCAVector;
	vector<double>::iterator curTempPCA, endTempPCA;
	vector<double*> pcaVector;
	vector<double*>::iterator curPCA, endPCA;
	vector<double> statMeanVector;
	vector<double> statSdVector;
	vector<double> statMaxVector;
	vector<double> statMinVector;
	vector<double> lambdaVector;
	vector<double> statGMVector;

	//open PCA_file
	ifstream pcaFile (LinearCombFileName.c_str()); // opening the file for reading
	if(!pcaFile)
	   throw TException("The Linear-Combination-File '" + LinearCombFileName + "' can not be read!");

	//read LinearCombFile
	while(pcaFile.good() && !pcaFile.eof()){
		isLine.clear();                   // create the new stream
		buf.read_line(pcaFile);
		buf=buf.extract_before_doubleSlash();   // remove the commentated part
		if(!buf.empty()){
			//allocate the line to a new stream for easy reading
			isLine << buf;
			//read the stat name
			name.read_to_delim(isLine);
			if(!name.empty()){
				statNamesVector.push_back(name);
			}
			//read max
			name.read_to_delim(isLine);
			if(!name.empty()){
				statMaxVector.push_back(name.toDouble());
			}
			//read min
			name.read_to_delim(isLine);
			if(!name.empty()){
				statMinVector.push_back(name.toDouble());
			}
			//read lambda
			name.read_to_delim(isLine);
			if(!name.empty()){
				lambdaVector.push_back(name.toDouble());
			}
			//read GM
			name.read_to_delim(isLine);
			if(!name.empty()){
				statGMVector.push_back(name.toDouble());
			}
			//read Mean after Boxcox
			name.read_to_delim(isLine);
			if(!name.empty()){
				statMeanVector.push_back(name.toDouble());
			}
			//read the sd after Box cox
			name.read_to_delim(isLine);
			if(!name.empty()){
			   statSdVector.push_back(name.toDouble());
			}
			//read PLS into temp vector
			tempPCAVector.clear();
			while(!name.empty()){
			   name.read_to_delim(isLine);
			   if(name.empty()) break;
			   tempPCAVector.push_back(name.toDouble());
			}
			//create PLS array
			if((int)tempPCAVector.size()==numLinearComb || numLinearComb==0){
			   numLinearComb=tempPCAVector.size();
			   double* tempArray=new double[numLinearComb];
			   curTempPCA=tempPCAVector.begin();
			   endTempPCA=tempPCAVector.end();
			   int i=0;
			   for(;curTempPCA!=endTempPCA; ++curTempPCA){
				  tempArray[i]=*curTempPCA;
				  ++i;
			   }
			   pcaVector.push_back(tempArray);
			} else
			   throw TException("Unequal number of columns among the lines in the Linear-Combination-File '" + LinearCombFileName + "'!", _FATAL_ERROR);
		}
	}
	pcaFile.close();

	curNameInInputFiles=namesInInputFiles.begin();
	endNameInInputFiles=namesInInputFiles.end();

	//write pls into an array.
	//The order of the array is similar to the order in the inputFile
	endPCA=pcaVector.end();
	pca=new double*[pcaVector.size()];
	max=new double[pcaVector.size()];
	min=new double[pcaVector.size()];
	lambda=new double[pcaVector.size()];
	gm=new double[pcaVector.size()];
	mean=new double[pcaVector.size()];
	sd=new double[pcaVector.size()];
	numStats=namesInInputFiles.size();
	statIsUsed=new int[numStats];
	int check;
	int k=0;
	for(int j=0;curNameInInputFiles!=endNameInInputFiles; ++curNameInInputFiles, ++j){
		curPCA=pcaVector.begin();
		check=false;
		for(int i=0;curPCA!=endPCA; ++curPCA, ++i){
			if(statNamesVector[i]==*curNameInInputFiles){
				pca[k]=new double[numLinearComb];
				pca[k]=*curPCA;
				max[k]=statMaxVector[i];
				min[k]=statMinVector[i];
				lambda[k]=lambdaVector[i];
				gm[k]=statGMVector[i];
				mean[k]=statMeanVector[i];
				sd[k]=statSdVector[i];
				check=true;
				statIsUsed[j]=k;
				++k;
			}
		}
		if(!check){
			statIsUsed[j]=-1;
		}
	}
	//check stats
	my_string missingStat=getFirstMissingStat();
	if(!missingStat.empty())
	   throw TException("The summary statistic '"+ missingStat +"' s missing in the simulated File, but required by the Linear-Combination-File '"+ LinearCombFileName + "'!", _FATAL_ERROR);
	simPCA=new double[numLinearComb];
	oldSimPCA=new double[numLinearComb];
}
//------------------------------------------------------------------------------
void TLinearCombBoxCox::calcObsLineraComb(vector<double> obsData){
	//order has to be the same as the names passed to the constructor!!!
	obs=new double[numStats];
	for(int i=0; i<numStats; ++i){
		if(statIsUsed[i]>=0){
			obs[statIsUsed[i]]=getBoxCox(obsData[i], i);;
		}
	}
	//calculate pca for obsData
	obsPCA=new double[numLinearComb];
	for(int i=0; i<numLinearComb; ++i){
		obsPCA[i]=getPCA(i, obs);
	}
}
//------------------------------------------------------------------------------
double TLinearCombBoxCox::getBoxCox(double simData, int stat){
	simData=1+(simData-min[stat])/(max[stat]-min[stat]);
	return (pow(simData, lambda[stat])-1)/(lambda[stat]*pow(gm[stat], (lambda[stat]-1)));
}

void TLinearCombBoxCox::calcSimDataPCA(double** simData){
	//first transform stats via boxcox
	int stat=0;
	for(int i=0; i<numStats; ++i){
		if(statIsUsed[i]>=0){
			*simData[i]=getBoxCox(*simData[i], stat);
			++stat;
		}
	}
	for(int i=0; i<numLinearComb; ++i){
		simPCA[i]=getPCA(i, simData);
	}
}
void TLinearCombBoxCox::calcSimDataPCA(double* simData){
	//first transform stats via boxcox
	int stat=0;
	for(int i=0; i<numStats; ++i){
		if(statIsUsed[i]>=0){
			simData[i]=getBoxCox(simData[i], stat);
			++stat;
		}
	}
	for(int i=0; i<numLinearComb; ++i){
		simPCA[i]=getPCA(i, simData);
	}
}
//------------------------------------------------------------------------------



