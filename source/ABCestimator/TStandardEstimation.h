//---------------------------------------------------------------------------

#ifndef TStandardEstimationH
#define TStandardEstimationH
#include "TEstimation.h"
//---------------------------------------------------------------------------
class TStandardEstimation:public TEstimation {
   public:
	  my_string trueParamName;
	  bool trueParamsAvailable;
	  int* trueParamNumInSimData;

	  TStandardEstimation(TParameters* parameters);
	  void standardize();
	  void performEstimations();
	  void makeRejection(int thisObsDataSet);
	  void printMarginalDensity(int thisObsDataSet);
	  float getquantileOfTrueParam(int param, float trueValue);
	  float getRMISEOfTrueParam(int param, float trueValue);
	  float getREMODEOfTrueParam(int param, float trueValue);
	  float getREMEANOfTrueParam(int param, float trueValue);
	  float getREMEDIANOfTrueParam(int param, float trueValue);




};


#endif
