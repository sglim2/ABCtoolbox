//---------------------------------------------------------------------------

#ifndef TSimDataH
#define TSimDataH
#include <strstream>
#include <fstream>
#include <iostream>
#include "my_cstring.h"
#include "newmat.h"
#include "newmatap.h"
#include "TParameters.h"
#include <vector>
#include "TObsData.h"
#include "math.h"
#include <limits>
//---------------------------------------------------------------------------
class TSimData{
   public:
	  my_string simFileName;
	  int numParams;
	  my_string paramString;
	  int numReadSims;
	  int numUsedSims;
	  my_string* paramNames;
	  float** simParams;
	  float** simStats;
	  float* statMeans;
	  float* statSDs;
	  float* statMaxima;
	  float* statMinima;
	  bool statMinMaxCalculated;
	  float* paramMaxima;
	  float* paramMinima;
	  TObsData* pointerToObsData;
	  Matrix statMatrix, paramMatrix;
	  float threshold;
	  bool statMeansCalculated, statSDsCalculated, statsAreStandardized, paramsStandardized;


	  TSimData(my_string simfilename, my_string sparams, int maxSimsToRead, TObsData* obsData);
	  ~TSimData(){
		  delete [] paramNames;
		  for(int j=0; j<numReadSims; ++j){
			 delete[] simParams[j];
			 delete[]simStats[j];
		  }
		  delete[] simParams;
		  delete[] simStats;
		  if(statMeansCalculated) delete[] statMeans;
		  if(statSDsCalculated) delete[] statSDs;
		  if(statMinMaxCalculated){
			  delete[] statMaxima;
			  delete[] statMinima;
		  }
		  if(paramsStandardized){
			  delete[] paramMaxima;
			  delete[] paramMinima;
		  }
	  }
	  void readSimFile(int maxSimsToRead);
	  void calculateDistances(float* obsValueArray, float* distances);
	  int getParamNumberFromName(my_string name);
	  void calculateStatisticsMeans();
	  void calculateStatisticsSDs();
	  void standardizeStatistics();
	  void standardizeParameters();
	  my_string checkIfstatsArePolymorphic();
	  void getMinMaxofStats();
	  void checkIfstatsAreWithinRange(float* obsValueArray);
	  int fillStatAndParamMatricesKeepAll();
	  int fillStatAndParamMatricesFromNumRetained(float* distances, int numToRetain);
	  int fillStatAndParamMatricesFromThreshold(float* distances, float gotThreshold, int numToRetain=-1);
	  void writeRetained(my_string outputPrefix, float* distances, int dataSet=-1);
	  float partitionP(float** a, int left, int right, int pivotIndex);
	  void quicksortP(float** a, int left, int right);
	 // float getParamMin(int param);
	 // float getParamMax(int param);
};



#endif
