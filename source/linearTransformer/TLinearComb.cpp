//---------------------------------------------------------------------------


#pragma hdrstop

#include "TLinearComb.h"
#include "TException.h"
//---------------------------------------------------------------------------

#pragma package(smart_init)

TLinearComb::TLinearComb(my_string gotLinearCombFileName, vector<my_string> NamesInInputFiles){
		LinearCombFileName=gotLinearCombFileName;
		namesInInputFiles=NamesInInputFiles;
		readLinearCombFile();
}
//------------------------------------------------------------------------------
void TLinearComb::readLinearCombFile(){
		numLinearComb=0;
		my_string name, buf, filename;
		strstream isLine;                   // create the new stream
		vector<double> tempPCAVector;
		vector<double>::iterator curTempPCA, endTempPCA;
		vector<double*> pcaVector;
		vector<double*>::iterator curPCA, endPCA;
		vector<double> statMeanVector;
		vector<double> statSdVector;

		//open LinearCombFile
		ifstream linearCombFile (LinearCombFileName.c_str()); // opening the file for reading
		if(!linearCombFile)
		   throw TException("The Linear-Combination-File '" + LinearCombFileName + "' can not be read!");

		//read PCA File
		while(linearCombFile.good() && !linearCombFile.eof()){
			isLine.clear();                   // create the new stream
			buf.read_line(linearCombFile);
			buf=buf.extract_before_doubleSlash();   // remove the commented part
			if(!buf.empty()){
				//allocate the line to a new stream for easy reading
				isLine << buf;
				//read the stat name
				name.read_to_delim(isLine);
				if(!name.empty()){
					statNamesVector.push_back(name);
				}
				//read the mean
				name.read_to_delim(isLine);
				if(!name.empty()){
				   statMeanVector.push_back(name.toDouble());
				}
				//read the sd
				name.read_to_delim(isLine);
				if(!name.empty()){
				   statSdVector.push_back(name.toDouble());
				}
				//read pca into temp vector
				tempPCAVector.clear();
				while(!name.empty()){
				   name.read_to_delim(isLine);
				   if(name.empty()) break;
				   tempPCAVector.push_back(name.toDouble());
				}
				//create pca array
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
		linearCombFile.close();

		//write pca into an array.
		//The order of the array is similar to the order in the inputFile
		//all stats have an array which contains zeros if the stat ist not specified in the PCA_file
		endPCA=pcaVector.end();
		pca=new double*[pcaVector.size()];
		mean=new double[pcaVector.size()];
		sd=new double[pcaVector.size()];
		numStats=namesInInputFiles.size();
		statIsUsed=new int[numStats];
		int check;
		int k=0;
		endNameInInputFiles=namesInInputFiles.end();
		curNameInInputFiles=namesInInputFiles.begin();
		for(int j=0;curNameInInputFiles!=endNameInInputFiles; ++curNameInInputFiles, ++j){
			curPCA=pcaVector.begin();
			check=false;
			for(int i=0;curPCA!=endPCA; ++curPCA, ++i){
				if(statNamesVector[i]==*curNameInInputFiles){
					pca[k]=new double [numLinearComb];
					pca[k]=*curPCA;
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
void TLinearComb::calcObsLineraComb(vector<double> obsData){
	//order has to be the same as the names passed to the constructor!!!
	obs=new double[numStats];
	for(int i=0; i<numStats; ++i){
		if(statIsUsed[i]>=0) obs[statIsUsed[i]]=obsData[i];
	}
	//calculate pca for obsData
	obsPCA=new double[numLinearComb];
	for(int i=0; i<numLinearComb; ++i) obsPCA[i]=getPCA(i, obs);
}
//------------------------------------------------------------------------------
int TLinearComb::checkName(my_string name){
	curStatName=statNamesVector.begin();
	endStatName=statNamesVector.end();
	for(; curStatName!=endStatName; ++curStatName){
		if(*curStatName==name) return true;
	}
	return false;
}
//------------------------------------------------------------------------------
int TLinearComb::getNumberOfFirstUsedStat(){
	for(int i=0; i<numStats; ++i){
		if(statIsUsed[i]>=0) return i;
	}
	return numStats+1;
}
//------------------------------------------------------------------------------
my_string TLinearComb::getFirstMissingStat(){
   curStatName=statNamesVector.begin();
   endStatName=statNamesVector.end();
   vector<my_string>::iterator curInputFileName, endInputFileName;
   endNameInInputFiles=namesInInputFiles.end();
   int check;
   for(;curStatName!=endStatName; ++curStatName){
	  curNameInInputFiles=namesInInputFiles.begin();
	  check=false;
	  for(;curNameInInputFiles!=endNameInInputFiles; ++curNameInInputFiles){
		if(*curNameInInputFiles==*curStatName){
		   check=true;
		   break;
		}
	  }
	  if(!check) return *curStatName;
   }
   return "";
}
//------------------------------------------------------------------------------
void TLinearComb::calcSimDataPCA(double** simData){
   for(int i=0; i<numLinearComb; ++i) simPCA[i]=getPCA(i, simData);
}
void TLinearComb::calcSimDataPCA(double* simData){
   for(int i=0; i<numLinearComb; ++i) simPCA[i]=getPCA(i, simData);
}
//------------------------------------------------------------------------------
void TLinearComb::writeHeader(ofstream& ofs){
   for(int i=0; i<numLinearComb; ++i){
	  ofs << "\t" << "LinearCombination_" << i;
   }
}
//------------------------------------------------------------------------------
void TLinearComb::writeSimPCA(ofstream& ofs){
   for(int i=0; i<numLinearComb; ++i){
	  ofs << "\t" << simPCA[i];
   }
}
//------------------------------------------------------------------------------
//Functions to calculate the linear combinations
//!!!!change both or none!!!
double TLinearComb::getPCA(int num, double* stats){
	 double pcaCalc=0;
	 int stat=0;
	 for(int i=0; i<numStats; ++i){
		if(statIsUsed[i]>=0){
			pcaCalc+= pca[stat][num]*(stats[i]-mean[stat])/sd[stat];
			++stat;
		}
	 }
	 return pcaCalc;
}
double TLinearComb::getPCA(int num, double** stats){
	double pcaCalc=0;
	int stat=0;
	for(int i=0; i<numStats; ++i){
		if(statIsUsed[i]>=0){
			pcaCalc+= pca[stat][num]*(*stats[i]-mean[stat])/sd[stat];
			++stat;
		}
	}
	return pcaCalc;
}
//------------------------------------------------------------------------------
void TLinearComb::saveOldValues(){
   for(int i=0; i<numLinearComb; ++i){
	  oldSimPCA[i]=simPCA[i];
   }
}
void TLinearComb::resetOldValues(){
   for(int i=0; i<numLinearComb; ++i){
	  simPCA[i]=oldSimPCA[i];
   }
}



