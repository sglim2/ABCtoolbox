//---------------------------------------------------------------------------
#ifndef TPriorH
#define TPriorH
#include "global.h"
#include "my_cstring.h"
#include <strstream>
//------------------------------------------------------------------------------
class TPrior{
   public:
	  my_string name;
	  bool isInt;
	  double curValue, oldValue;

	  TPrior(my_string Name, bool isInt);
	  void makeCurValueInt();
	  void writeValue(const double& value, ofstream& file);
	  void writeCurValue(ofstream& file);
	  void writeHyperPriorGamma(const double& arg, ofstream& file);
	  void writeHyperPriorBeta(ofstream& file);
	  void writeHyperPriorNormal(const double& stdev, ofstream& file);
	  void writeHyperPriorNormalPositive(const double& stdev, ofstream& file);
	  void writeHyperPriorLognormal(const double& stdev, ofstream& file);
	  void writeHyperPriorLognormalParametersInLog10Scale(const double& stdev, ofstream& file); // base 10!!!
	  void saveOldValue();
	  void resetOldValue();
	  void setCurValue(double newValue);
};
class TSimplePrior:public TPrior{
   public:
	  double mcmcStep;
	  double upperLimit, lowerLimit;

	  TSimplePrior(my_string Name, double min, double max, bool IsInt);
	  virtual void changeCurrentValue();
	  virtual double getPriorDensity();
	  virtual void changeCurrentValueWithLimits(const double& Min, const double& Max);
	  virtual double getPriorDensityFromValue(const double& value);
};
class TCombinedPrior: public TPrior{
   public:
	  my_string equation;
	  TCombinedPrior(my_string Name, my_string Equation, bool IsInt, ofstream* gotLogFile);
};
//------------------------------------------------------------------------------
class TUniformPrior:public TSimplePrior{
   public:
	  TUniformPrior(my_string Name, double min, double max, bool IsInt);
	   void changeCurrentValueWithLimits(const double& Min, const double& Max);
	   void changeCurrentValue();
	   double getPriorDensity();
	   double getPriorDensity(const double& value);
	   double getPriorDensityFromValue(const double& value);
};
class TLogUniformPrior:public TSimplePrior{
   public:
      double inverseTruncatedArea;
	  TLogUniformPrior(my_string Name, double min, double max, bool IsInt);
	  void changeCurrentValue();
	  double getPriorDensity();
	  void changeCurrentValueWithLimits(const double& Min, const double& Max);
	  double getPriorDensity(const double& value);
	  double getPriorDensityFromValue(const double& value);
};
class TNormalPrior:public TSimplePrior{
   public:
   double mean, sigma, inverseTruncatedArea;
	  TNormalPrior(my_string Name, double min, double max, double Mean, double Sigma, bool IsInt);
	  void changeCurrentValue();
	  double cumulativeDistributionFunction(double x);
	  double complementaryErrorFunction(double x);
	  double complementaryErrorCheb(double x);
	  double getPriorDensity();
	  void changeCurrentValueWithLimits(const double& Min, const double& Max);
	  double getPriorDensity(const double& value);
	  double getPriorDensityFromValue(const double& value);
};
class TLogNormalPrior:public TNormalPrior{
   public:
	  TLogNormalPrior(my_string Name, double min, double max, double Mean, double Sigma, bool IsInt);
	  void changeCurrentValue();
	  double getPriorDensity();
	  void changeCurrentValueWithLimits(const double& Min, const double& Max);
	  double complementaryErrorFunction(double x);
	  double complementaryErrorCheb(double x);
	  double getPriorDensity(const double& value);
	  double getPriorDensityFromValue(const double& value);
};
class TGammaPrior:public TSimplePrior{
   public:
	  TGammaPrior(my_string Name, double FirstParameter, double SecondParameter, bool IsInt);
	  void changeCurrentValue();
	  double getPriorDensity();
	  void changeCurrentValueWithLimits(const double& Min, const double& Max);
	  double getPriorDensity(const double& value);
	  double getPriorDensityFromValue(const double& value);
};
//------------------------------------------------------------------------------
#endif
