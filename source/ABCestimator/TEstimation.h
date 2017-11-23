//---------------------------------------------------------------------------

#ifndef TEstimationH
#define TEstimationH
#include "TObsData.h"
#include "TSimData.h"
#include <stdlib.h>
#include <time.h>

//---------------------------------------------------------------------------
class TEstimation {
   public:
	  TObsData* myObsData;
	  TSimData* mySimData;
	  int numToRetain;
	  float threshold;
	  bool retainedDefined;
	  bool thresholdDefined;
	  bool writeRetainedSimulations;
	  ofstream mdFile, propCommonSurfacePriorPosterior;
	  bool distFromDifferentStats;
	  int calcObsPValue;
	  TObsData* myObsDataForDist;
	  TSimData* mySimDataForDistForDist;

	  float diracPeakWidth;
	  float step; //step between posterior density points
	  int posteriorDensityPoints;
	  ColumnVector obsValues;
	  Matrix C;
	  ColumnVector c_zero;
	  Matrix CHat;
	  Matrix RHat;
	  Matrix SigmaS;
	  Matrix SigmaSInv;
	  Matrix SigmaTheta;
	  Matrix SigmaThetaInv;
	  Matrix posteriorMatrix;
	  Matrix retainedMatrix;
	  Matrix priorMatrix;
	  ColumnVector theta;
	  int numParams;
	  my_string outputPrefix;
	  bool standardizeStats;
	  bool writeSmoothedsimsUsed;
	  float* posteriorHDIMin;
	  float* posteriorHDIMax;
	  float* fmFromRetained;

	  TEstimation(TParameters* parameters);
	  ~TEstimation(){
		  delete myObsData;
		  delete mySimData;
		  if(distFromDifferentStats){
			  delete myObsDataForDist;
			  delete mySimDataForDistForDist;
		  }
	  }
	  virtual void standardize();
	  void preparePosteriorDensityPoints();
	  virtual void performEstimations();
	  void preparePrior();
	  void calculatePosterior();
	  virtual Matrix calculate_T_j();
	  virtual ColumnVector calculate_v_j(ColumnVector theta_j);
	  void performLinearRegression();
	  void prepareDiracPeaks();
	  float theatParameterScale(int param, int k);
	  void writePosteriorFile(my_string filenameTag="");
	  float getPositionofMode(int param);
	  float getPosteriorMean(int param);
	  float getPosteriorquantile(int param, float q);
	  void writePosteriorCharacteristics(my_string filenameTag="");
	  void writeSmoothedSimsUsed(my_string filenameTag="");
	  double getArea(ColumnVector posterior, ColumnVector theta);
	  float getCommonSurfaceRetainedPosterior(int p);
	  void calculateSmoothedRetainedMatrix();
	  void getPosteriorHDI(int param, float q, float* lower, float* upper);
	  void calculateMarginalDensitiesOfRetained();
};


#endif
