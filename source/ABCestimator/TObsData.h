//---------------------------------------------------------------------------

#ifndef TObsDataH
#define TObsDataH
//---------------------------------------------------------------------------
#include <strstream>
#include <fstream>
#include <iostream.h>
#include "my_cstring.h"
#include "newmat.h"
#include "newmatio.h"
#include "TParameters.h"
#include <vector>

//---------------------------------------------------------------------------
class TObsDataSet{
	public:
	   float* myObsValues;
	   float* myStandardizedObsValues;
	   float* myTrueParams;
	   int internalCounter;
	   int internalTrueParamCounter;
	   int numObsDataSets;
	   int numStats;
	   int numTrueParams;

	   TObsDataSet(int gotNumStats);
	   TObsDataSet(){
		   myObsValues=NULL;
		   myStandardizedObsValues=NULL;
		   internalCounter=0;
		   numObsDataSets=0;
		   numStats=0;
	   }

	   ~TObsDataSet(){
		  // delete [] myObsValues;
		  // delete [] myStandardizedObsValues;
	   }
	   void add(float value);
	   void setNumTrueParams(int value);
	   void addTrueParams(float value);
	   void standardize(float* means, float*sds);
	   void getStandardizedObsValuesIntoColumnVector(ColumnVector& colVector);
	   void getObsValuesIntoColumnVector(ColumnVector& colVector);
};

class TObsData{
   public:
	  my_string obsFileName;
	  vector<my_string> obsNameVector;
	  vector<my_string> trueNameVector;
	  vector<my_string>::iterator curObsNameVector, endObsNameVector;
	  int numStats;
	  int numTrueParams;
	  int numObsDataSets;
	  vector<TObsDataSet*> obsDataSets;

	  Matrix obsValuesMatrix;
	  bool hasBeenStandardized;

      TObsData(my_string obsfilename);
      ~TObsData(){
    	  //obsNameVector.clear();
    	  //obsDataSets.clear();
      }
	  void readObsFile();
	  void readTrueFile(my_string truefilename);
	  int getStatColNumberFromName(my_string name);
	  void fillObsDataMatrix(int dataSet);
	  void standardizeObservedValues(float* means, float*sds);
	  float* getObservedValues(int num);
	  void getObsValuesIntoColumnVector(int num, ColumnVector& colVector);
	  void writeTrueParamsheader(ofstream* file);
	  float getTrueParam(int num, int trueparam);
};

#endif
