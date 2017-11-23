//---------------------------------------------------------------------------

#pragma hdrstop

#include "TDataVector.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
//------------------------------------------------------------------------------
TDataVector::TDataVector(TParameters* gotParameters, TInputFileVector* gotInputFiles, ofstream* gotlogfile){
   logFile=gotlogfile;
   inputFiles=gotInputFiles;
   sumStatsTempFileName=gotParameters->getParameter("sumStatFile", false);
   //standardization of linear combinations?
   stdLinearCombForDist=gotParameters->getdoubleParameter("stdLinearCombForDist",0);
   if(stdLinearCombForDist) *logFile << "Linear combinations are standardized for the distance calculation" << endl;

   if(sumStatsTempFileName=="") sumStatsTempFileName="summary_stats_temp.txt";
   *logFile << "The following file is read to get the summary statistics: '" <<  sumStatsTempFileName << "'" << endl;

   //read obsnames
   my_string obsname=gotParameters->getParameter("obsName");
   vector<my_string> vecObsNames;
   while(!obsname.empty()){
	   vecObsNames.push_back(obsname.extract_sub_str(*";"));
	   obsname.remove(0,1);
   }

   //read names of files with simulated data
	  my_string simDataName=gotParameters->getParameter("simDataName", false);
	  vector<my_string> vecSimDataName;
	  vector<my_string>::iterator curSimDataName;
	  bool simDataNameGiven=false;
	  if(!simDataName.empty()){
		 while(!simDataName.empty()){
			vecSimDataName.push_back(simDataName.extract_sub_str(*";"));
			simDataName.remove(0,1);
		 }
		 if(vecSimDataName.size()==1){
			 //check for tag...
			 simDataName=vecSimDataName[0];
			 int pos=simDataName.find("SIMINPUTNAME");
			 if(pos>=0){
				 my_string before=simDataName.extract_before(pos);
				 my_string after=simDataName.extract_after(pos+11);
				 //construct from siminputfiles
				 vecSimDataName=inputFiles->getVectorOfInputFileNames();
				 curSimDataName=vecSimDataName.begin();
				 for(;curSimDataName!=vecSimDataName.end();++curSimDataName){
					 curSimDataName->remove_file_extension();
					 *curSimDataName=before+*curSimDataName + after;
				 }
			 }
		 }
		 simDataNameGiven=true;
		 if(vecObsNames.size()!=vecSimDataName.size())
	 		 throw TException("The number of .obs files and files with simulated data are unequal!", _FATAL_ERROR);
	  }

	  //create vector for TData objects
	  vector<my_string>::iterator curObsName;
	  curObsName=vecObsNames.begin();
	  if(simDataNameGiven){
		  curSimDataName=vecSimDataName.begin();
		  for(;curObsName!=vecObsNames.end() && curSimDataName!=vecSimDataName.end(); ++curObsName, ++curSimDataName){
			  vecDataObjects.push_back(TData(*curObsName, *curSimDataName, gotlogfile));
		  }
	  } else {
		  for(;curObsName!=vecObsNames.end(); ++curObsName){
		  	  vecDataObjects.push_back(TData(*curObsName, gotlogfile));
		  }
	  }
	  endDataObject=vecDataObjects.end();
	  numDataObjects=vecDataObjects.size();

	  //read sumstat programs
	  int i=0;


	  //since a simulation was already done, try to find the simulated data files "SIMDATANAME"
	  if(simDataNameGiven){
		  my_string foundThere;
		  curDataObject=vecDataObjects.begin();
		  for(;curDataObject!=endDataObject; ++curDataObject){
			  foundThere=locateFileRecursively(curDataObject->simDataFileName, ".");
			  if(foundThere=="") throw TException("The file containing the simulated data '" + curDataObject->simDataFileName + "' has not been found!", _FATAL_ERROR);
			  *logFile << "File '" << curDataObject->simDataFileName << "' is located here: '" << foundThere << "':" << endl;
			  curDataObject->simDataFileName=foundThere;
		  }
	  }

	  //initialize sumstat program
	  my_string sumstatprogram=gotParameters->getParameter("sumStatProgram", false);
	  if(sumstatprogram==""){
		  useSumstatProgram=false;
	  	  *logFile << "sumStatProgram not defined -> simulation program or script writes the file: '" <<  sumStatsTempFileName << "' with the statistics." << endl;
	  } else {
		  useSumstatProgram=true;
		  vector<my_string> vecSumstatprogram;
		  while(!sumstatprogram.empty()){
			  vecSumstatprogram.push_back(sumstatprogram.extract_sub_str(*";"));
			  sumstatprogram.remove(0,1);
		  }
		  //read param
		  my_string parameterSumStatProgram=gotParameters->getParameter("sumStatParam");
		  vector<my_string> sumStatParam;
		  while(!parameterSumStatProgram.empty()){
			  sumStatParam.push_back(parameterSumStatProgram.extract_sub_str(*";"));
		 	  parameterSumStatProgram.remove(0,1);
		  }
		  curDataObject=vecDataObjects.begin();
		  if(vecSumstatprogram.size()==1 && sumStatParam.size()==1){
			  *logFile << "Using SumStat program: '" <<  vecSumstatprogram[0] << "'" << endl;
			  for(;curDataObject!=endDataObject; ++curDataObject){
			  	 curDataObject->initializeSumstatProgramm(vecSumstatprogram[0], sumStatParam[0], sumStatsTempFileName);
			  }
		  } else {
			  i=0;
			  if(vecSumstatprogram.size()==1 && sumStatParam.size()==vecSimDataName.size()){
				  *logFile << "Using SumStat program: '" <<  vecSumstatprogram[0] << "'" << endl;
				  for(;curDataObject!=endDataObject; ++curDataObject, ++i){
					  curDataObject->initializeSumstatProgramm(vecSumstatprogram[0], sumStatParam[i], sumStatsTempFileName);
				  }
			  } else {
				  if(vecSumstatprogram.size()==vecSimDataName.size() && sumStatParam.size()==vecSimDataName.size()){
					  *logFile << "Using SumStat programs:" << endl;
					  for(;curDataObject!=endDataObject; ++curDataObject, ++i){
						  *logFile << "   - '" << vecSumstatprogram[i]<< "'" << endl;
				  		  curDataObject->initializeSumstatProgramm(vecSumstatprogram[i], sumStatParam[i], sumStatsTempFileName);
					  }
				  } else throw TException("Wrong number of programs to calculate summary statistics or defined argument lists for these programs!", _FATAL_ERROR);
			  }
		  }
	  }
	  summarizeDataObjects();

	//scripts before and after sumstat calculation....
	#ifdef _GCC_
    //script before
    my_string script=gotParameters->getParameter("launchBeforeSS", 0);
    if(script!=""){
    	//check what arguments have to be passed
    	my_string temp=gotParameters->getParameter("launchBeforeSSParam", false);
    	//split argument lists
    	vector<my_string> scriptParamList;
    	while(!temp.empty()){
    		scriptParamList.push_back(temp.extract_sub_str(*";"));
    		temp.remove(0,1);
    	}
    	//initialize in each TData
    	curDataObject=vecDataObjects.begin();
    	if(scriptParamList.size()==0){
    		for(;curDataObject!=endDataObject; ++curDataObject){
    			curDataObject->initializeScriptBeforeSSCalc(script, "", sumStatsTempFileName);
    		}
    	} else {
    		if(scriptParamList.size()==1){
        		for(;curDataObject!=endDataObject; ++curDataObject){
        			curDataObject->initializeScriptBeforeSSCalc(script, scriptParamList[0], sumStatsTempFileName);
        		}
    		} else {
    			if(scriptParamList.size()==vecDataObjects.size()){
    				i=0;
            		for(;curDataObject!=endDataObject; ++curDataObject, ++i){
            			curDataObject->initializeScriptBeforeSSCalc(script, scriptParamList[i], sumStatsTempFileName);
            		}
    			} else throw TException("Unequal number of files with simulated data sets and parameter specifications for the script/program to launch before the calculation of summary statistics!", _FATAL_ERROR);
    		}
    	}
    	launchScriptBeforeSS=true;
        *logFile << "Before each launch of the program calculating summary statistics, the following script is launched: " <<  script << endl;
    } else launchScriptBeforeSS=false;

   	//script after
    script=gotParameters->getParameter("launchAfterSS", 0);
    if(script!=""){
    	//check what arguments have to be passed
    	my_string temp=gotParameters->getParameter("launchAfterSSParam", false);
    	//split argument lists
    	vector<my_string> scriptParamList;
    	while(!temp.empty()){
    		scriptParamList.push_back(temp.extract_sub_str(*";"));
    		temp.remove(0,1);
    	}
    	//initialize in each TData
    	curDataObject=vecDataObjects.begin();
    	if(scriptParamList.size()==0){
    		for(;curDataObject!=endDataObject; ++curDataObject){
    			curDataObject->initializeScriptAfterSSCalc(script, "", sumStatsTempFileName);
    		}
    	} else {
    		if(scriptParamList.size()==1){
        		for(;curDataObject!=endDataObject; ++curDataObject){
        			curDataObject->initializeScriptAfterSSCalc(script, scriptParamList[0], sumStatsTempFileName);
        		}
    		} else {
    			if(scriptParamList.size()==vecDataObjects.size()){
    				i=0;
            		for(;curDataObject!=endDataObject; ++curDataObject, ++i){
            			curDataObject->initializeScriptAfterSSCalc(script, scriptParamList[i], sumStatsTempFileName);
            		}
    			} else throw TException("Unequal number of files with simulated data sets and parameter specifications for the script/program to launch after the calculation of summary statistics!", _FATAL_ERROR);
    		}
    	}
    	launchScriptAfterSS=true;
        *logFile << "After each launch of the program calculating summary statistics, the following script is launched: " <<  script << endl;
    } else launchScriptAfterSS=false;
	#else
		launchScriptAfterSS=false;
		launchScriptBeforeSS=false;
    #endif

}
//------------------------------------------------------------------------------
void TDataVector::summarizeDataObjects(){
	//summarize some features of the TData objects
	  //numObsData
	  numObsData=0;
	  curDataObject=vecDataObjects.begin();
	  for(;curDataObject!=endDataObject; ++curDataObject){
		 numObsData+=curDataObject->numObsData;
	  }

	  //create an array with pointers to all obs Data
	  pointersToObsData=new double*[numObsData];
	  curDataObject=vecDataObjects.begin();
	  int i=0;
	  for(;curDataObject!=endDataObject; ++curDataObject){
		 for(int j=0; j<curDataObject->numObsData; ++j){
			pointersToObsData[i]=&curDataObject->obsData[j];
			++i;
		 }
	  }
	  //create an array with pointers to all ObsVariances Data
	  pointersToVariances=new double*[numObsData];
	  curDataObject=vecDataObjects.begin();
	  i=0;
	  for(;curDataObject!=endDataObject; ++curDataObject){
		 for(int j=0; j<curDataObject->numObsData; ++j){
			pointersToVariances[i]=&curDataObject->obsDataVariance[j];
			++i;
		 }
	  }
	  //create an Array with all ObsDataNames (with Prefix added if more than one object)
	  obsDataNamesWithPrefix=new my_string[numObsData];
	  curDataObject=vecDataObjects.begin();
	  i=0;
	  int numObs=0;
	  if(vecDataObjects.size()>1){
		  for(;curDataObject!=endDataObject; ++curDataObject, ++numObs){
			  for(int j=0; j<curDataObject->numObsData; ++j){
				  my_string a= "Obs" + (my_string) numObs;
				  obsDataNamesWithPrefix[i]=a + "_" + curDataObject->obsDataName[j];
				  ++i;
			  }
		  }
	  } else {
		  for(int j=0; j<curDataObject->numObsData; ++j){
			  obsDataNamesWithPrefix[i]=curDataObject->obsDataName[j];
			  ++i;
		  }
	  }
}
//------------------------------------------------------------------------------
void TDataVector::calculateSumStats(int simnum){
   curDataObject=vecDataObjects.begin();
   for(;curDataObject!=endDataObject; ++curDataObject){
	   if(launchScriptBeforeSS) curDataObject->launchScriptBeforeSSCalc(simnum);
	   if(useSumstatProgram) curDataObject->calculateSumStats(simnum);
	   if(launchScriptAfterSS) curDataObject->launchScriptAfterSSCalc(simnum);
	   curDataObject->readSimData(sumStatsTempFileName);
   }
}
//------------------------------------------------------------------------------
my_string TDataVector::locateFileRecursively(my_string file, my_string dir){
	DIR *pdir;
	struct dirent *pent;
	my_string name;
	my_string ret;
	pdir=opendir(dir.c_str());
	if(pdir==NULL){
		closedir(pdir);
		return "";
	}
	while ((pent=readdir(pdir))){
		name=(my_string) pent->d_name;
		//if file found
		if(name==file){
			closedir(pdir);
			return file;
		}
		if(name.find(".")<0){
			ret=locateFileRecursively(file, dir+"/"+name);
			if(ret!=""){
				closedir(pdir);
				return name+"/"+ret;
			}
		}
	}
	closedir(pdir);
	return "";
}
//------------------------------------------------------------------------------
void TDataVector::initializeDataObjects(){
   curDataObject=vecDataObjects.begin();
   pointersToSimData=new double*[numObsData];
   int i=0;
   for(;curDataObject!=endDataObject; ++curDataObject){
	   if(launchScriptBeforeSS) curDataObject->launchScriptBeforeSSCalc(-1);
	   if(useSumstatProgram) curDataObject->calculateSumStats(-1);
	   if(launchScriptAfterSS) curDataObject->launchScriptAfterSSCalc(-1);
	   curDataObject->readInitialSimData(sumStatsTempFileName);
	   curDataObject->saveOldValues();
	   for(int j=0; j<curDataObject->numObsData; ++j){
		   pointersToSimData[i]=curDataObject->simDataPointers[j];
		   ++i;
	   }
   }
}
//------------------------------------------------------------------------------
//Functions to calculate the distances
//!!!!if one is changed, change all others as well!!!!!!
double TDataVector::calculateDistance(){
   double dis=0;
   for(int i=0; i<numObsData; ++i){
      dis+=(*pointersToObsData[i]-*pointersToSimData[i])*(*pointersToObsData[i]-*pointersToSimData[i]) / *pointersToVariances[i];
   }
   return sqrt(dis);
}
double TDataVector::calculateDistance(TLinearComb* pcaObject){
   //first calculate pca's
   pcaObject->calcSimDataPCA(pointersToSimData);
   //then calculate distance
   return calculateDistanceLinearComb(pcaObject);
}
double TDataVector::calculateDistance(double* data){
   double dis=0;
   for(int i=0; i<numObsData; ++i){
	  dis+=(*pointersToObsData[i]-data[i])*(*pointersToObsData[i]-data[i]) / *pointersToVariances[i];
   }
   return sqrt(dis);
}
double TDataVector::calculateDistance(double* data, TLinearComb* pcaObject){
   //first calculate pca's
   pcaObject->calcSimDataPCA(data);
   //then calculate distance
   return calculateDistanceLinearComb(pcaObject);
}
double TDataVector::calculateDistanceLinearComb(TLinearComb* pcaObject){
   double dis=0;
   if(stdLinearCombForDist){
	   for(int i=0; i<pcaObject->numLinearComb; ++i){
	      dis+=(pcaObject->obsPCA[i]-pcaObject->simPCA[i])*(pcaObject->obsPCA[i]-pcaObject->simPCA[i])/pcaObject->linearCombVariances[i];
	   }
   } else {
	   for(int i=0; i<pcaObject->numLinearComb; ++i){
	      dis+=(pcaObject->obsPCA[i]-pcaObject->simPCA[i])*(pcaObject->obsPCA[i]-pcaObject->simPCA[i]);
	   }
   }
   return sqrt(dis);
}
//------------------------------------------------------------------------------
my_string TDataVector::getObsFileNames(){
   curDataObject=vecDataObjects.begin();
   my_string allObsFileNames=curDataObject->obsFileName;
   ++curDataObject;
   for(;curDataObject!=endDataObject; ++curDataObject){
	  allObsFileNames+=", "+curDataObject->obsFileName;
   }
   return allObsFileNames;
}
//------------------------------------------------------------------------------
void TDataVector::saveOldValues(){
   curDataObject=vecDataObjects.begin();
   for(;curDataObject!=endDataObject; ++curDataObject){
	  curDataObject->saveOldValues();
   }

}
//------------------------------------------------------------------------------
void TDataVector::resetOldValues(){
   curDataObject=vecDataObjects.begin();
   for(;curDataObject!=endDataObject; ++curDataObject){
	  curDataObject->resetOldValues();
   }

}
//------------------------------------------------------------------------------
void TDataVector::writeHeader(ofstream& ofs){
	for(int i=0; i<numObsData; ++i){
		ofs << "\t" << obsDataNamesWithPrefix[i];
	}
}
void TDataVector::writeHeaderOnlyOneDataObject(int thisObs, ofstream& ofs){
   vecDataObjects[thisObs].writeHeader(ofs);
}
//------------------------------------------------------------------------------
void TDataVector::writeData(ofstream& ofs){
   for(int i=0; i<numObsData; ++i){
	  ofs << "\t" << *pointersToSimData[i];
   }
}
void TDataVector::writeDataOnlyOneDataObject(int thisObs, ofstream& ofs){
	vecDataObjects[thisObs].writeData(ofs);
}
//------------------------------------------------------------------------------
bool TDataVector::matchColumns(vector<my_string> & colNames, int numColNames, int** pointerToArrayToFill){
   curDataObject=vecDataObjects.begin();
   int* arrayToFill=*pointerToArrayToFill;
   for(int k=0; k<numObsData; ++k){
	  for(int j=0; j<(numColNames+1); ++j){
		 if(j==numColNames)
			throw TException("The data column '" + curDataObject->obsDataName[k] + "' is missing in the calibration File!!!", _FATAL_ERROR);
		 if(obsDataNamesWithPrefix[k]==colNames[j]){
			arrayToFill[k]=j;
			break;
		 }
	  }
   }
   return true;
}
//------------------------------------------------------------------------------
vector<my_string> TDataVector::getNamesVector(){
   vector<my_string> statNames;
   for(int i=0; i < numObsData; ++i)
	   statNames.push_back(obsDataNamesWithPrefix[i]);
   return statNames;
}
//------------------------------------------------------------------------------
vector<double> TDataVector::getObsDataVector(){
   vector<double> obsData;
   for(int i=0; i < numObsData; ++i)
	   obsData.push_back(*pointersToObsData[i]);
   return obsData;
}
//------------------------------------------------------------------------------
void TDataVector::fillObsDataArray(double* & obsData){
   obsData=new double[numObsData];
   for(int i=0; i < numObsData; ++i)
	   obsData[i]=*pointersToObsData[i];
}
//------------------------------------------------------------------------------
