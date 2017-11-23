//---------------------------------------------------------------------------


#include "TPrior.h"
#include <strstream>
//------------------------------------------------------------------------------
 //TPrior
TPrior::TPrior(my_string Name, bool IsInt){
   name=Name;
   isInt=IsInt;
   curValue=0;
   oldValue=0;
}
void TPrior::makeCurValueInt(){
   if(curValue-(int)curValue<0.5) curValue=(int)curValue;
   else                           curValue=1+(int)curValue;
   if(curValue<1) curValue=1;
}
void TPrior::writeValue(const double& value, ofstream& file){
   //attention about sicientific notations (integers cant read them)
   if(isInt) file << (long) value;
   else file << value;
}
void TPrior::writeCurValue(ofstream& file){
   writeValue(curValue, file);
}
void TPrior::writeHyperPriorGamma(const double& arg, ofstream& file){
   writeValue(gamma_dev(arg, arg/curValue), file);
}
void TPrior::writeHyperPriorBeta(ofstream& file){
	if(curValue<0.001) writeValue(0.0, file);
	else {
	   if(curValue>1) throw TException("Hyperprior '"+name+"': mean of the beta distribution > 1!", _FATAL_ERROR);
	   double a=0.5+199*curValue;
	   writeValue(BetaRandom(a, a*(1-curValue)/curValue), file);
	}
}
void TPrior::writeHyperPriorNormal(const double& stdev, ofstream& file){
   writeValue(NormalRandom(curValue, stdev), file);
}
void TPrior::writeHyperPriorNormalPositive(const double& stdev, ofstream& file){
	float val=NormalRandom(curValue, stdev);
	while(val<=0) val=NormalRandom(curValue, stdev);
    writeValue(val, file);
}
void TPrior::writeHyperPriorLognormal(const double& stdev, ofstream& file){
	   if(stdev<=0.) throw TException ("Hyperprior '"+name+"': stdev of the log normal distribution <=0!", _FATAL_ERROR);
	   double mean=log(curValue) - 0.5 * log(1 + ( (stdev*stdev)/(curValue*curValue) ));
	   double sigma=sqrt( log( ((stdev*stdev)/(curValue*curValue)) +1 ) );
	   writeValue(exp(NormalRandom(mean, sigma)), file);
}
void TPrior::writeHyperPriorLognormalParametersInLog10Scale(const double& stdev, ofstream& file){ //base 10!!!!
	if(stdev<=0.) throw TException ("Hyperprior '"+name+"': stdev of the log normal distribution with parameters in log scale <=0!", _FATAL_ERROR);
	writeValue(pow(10,NormalRandom(curValue, stdev)), file);
}

void TPrior::saveOldValue(){
   oldValue=curValue;
}
void TPrior::resetOldValue(){
   curValue=oldValue;
}
void TPrior::setCurValue(double newValue){
   curValue=newValue;
   if (isInt) makeCurValueInt();
}

//------------------------------------------------------------------------------
//TCombinedPrior
TCombinedPrior::TCombinedPrior(my_string Name, my_string Equation, bool IsInt, ofstream* gotLogFile):TPrior(Name, IsInt){
  equation=Equation;
}
//------------------------------------------------------------------------------
//TSimplePiror
TSimplePrior::TSimplePrior(my_string Name, double min, double max, bool IsInt):TPrior(Name, IsInt){
   lowerLimit=min;
   upperLimit=max;
   if(lowerLimit>upperLimit) throw TException ("Prior '" + name + "' initialized with min > max!", _FATAL_ERROR);
}
void TSimplePrior::changeCurrentValue(){
   changeCurrentValueWithLimits(lowerLimit, upperLimit);
}
double TSimplePrior::getPriorDensity(){
   return getPriorDensityFromValue(curValue);
}
void TSimplePrior::changeCurrentValueWithLimits(const double& Min, const double& Max){}
double TSimplePrior::getPriorDensityFromValue(const double& value){
   return 0.0;
}
//------------------------------------------------------------------------------
//uniform prior
TUniformPrior::TUniformPrior(my_string Name, double min, double max, bool IsInt):TSimplePrior(Name, min, max, IsInt){}
//standard function with no arguments
void TUniformPrior::changeCurrentValue(){
   changeCurrentValueWithLimits(lowerLimit, upperLimit);
}
double TUniformPrior::getPriorDensity(){
   return getPriorDensityFromValue(curValue);
}
void TUniformPrior::changeCurrentValueWithLimits(const double& Min, const double& Max){
   curValue=UniformRandom(Min, Max);
   if (isInt) makeCurValueInt();
}
double TUniformPrior::getPriorDensityFromValue(const double& value){
   if(value>upperLimit || value<lowerLimit) return 0.0;
   return (1/(upperLimit-lowerLimit));
}
//------------------------------------------------------------------------------
//loguniform prior
TLogUniformPrior::TLogUniformPrior(my_string Name, double min, double max, bool IsInt):TSimplePrior(Name, min, max, IsInt){
   if( lowerLimit < 1e-15) lowerLimit= 1e-15;
}
void TLogUniformPrior::changeCurrentValue(){
   changeCurrentValueWithLimits(lowerLimit, upperLimit);
}
double TLogUniformPrior::getPriorDensity(){
   return getPriorDensityFromValue(curValue);
}
void TLogUniformPrior::changeCurrentValueWithLimits(const double& Min, const double& Max){
   curValue=exp(UniformRandom(log(Min), log(Max)));
   if (isInt) makeCurValueInt();
}
double TLogUniformPrior::getPriorDensityFromValue(const double& value){
   if(value>upperLimit || value<lowerLimit) return 0.0;
   return 1.0/(value*(log(upperLimit) - log(lowerLimit)));
}
//------------------------------------------------------------------------------
//normal prior
TNormalPrior::TNormalPrior(my_string Name, double min, double max, double Mean, double Sigma, bool IsInt):TSimplePrior(Name, min, max, IsInt){
   mean=Mean;
   sigma=Sigma;
   if(sigma<=0.) throw TException ("Normal prior '" + name + "' initialized with stdev <=0!", _FATAL_ERROR);
   //now take care of the truncation
   inverseTruncatedArea=0;
   //lower truncation
   if(lowerLimit>0) inverseTruncatedArea+=cumulativeDistributionFunction(lowerLimit);
   //upper truncation
   inverseTruncatedArea+=1-cumulativeDistributionFunction(lowerLimit);
   inverseTruncatedArea=1/inverseTruncatedArea;
}
double TNormalPrior::cumulativeDistributionFunction(double x){
  if(x==0.) return 0.;
  return 0.5*complementaryErrorFunction(-0.707106781186547524*(x-mean)/sigma);
}
double TNormalPrior::complementaryErrorFunction(double x){
   //see book "numerical recipes"
   if(x>=0.) return complementaryErrorCheb(x);
   else return 2.0- complementaryErrorCheb(-x);
}
double TNormalPrior::complementaryErrorCheb(double x){
   //see book "numerical recipes"
   double coef[28]={-1.3026537197817094, 6.4196979235649026e-1,
   1.9476473204185836e-2,-9.561514786808631e-3,-9.46595344482036e-4,
   3.66839497852761e-4,4.2523324806907e-5,-2.0278578112534e-5,
   -1.624290004647e-6,1.303655835580e-6,1.5626441722e-8,-8.5238095915e-8,
   6.529054439e-9,5.059343495e-9,-9.91364156e-10,-2.27365122e-10,
   9.6467911e-11, 2.394038e-12,-6.886027e-12,8.94487e-13, 3.13092e-13,
   -1.12708e-13,3.81e-16,7.106e-15,-1.523e-15,-9.4e-17,1.21e-16,-2.8e-17};
   int j;
   double t, ty, tmp, d=0., dd=0.;
   t=2./(2.+x);
   ty=4.*t-2;
   for(j=27;j>0;--j){
	  tmp=d;
	  d=ty*d-dd+coef[j];
	  dd=tmp;
   }
   return t*exp(-x*x+0.5*(coef[0]+ty*d)-dd);
}
void TNormalPrior::changeCurrentValue(){
   changeCurrentValueWithLimits(lowerLimit, upperLimit);
}
double TNormalPrior::getPriorDensity(){
   return getPriorDensityFromValue(curValue);
}
void TNormalPrior::changeCurrentValueWithLimits(const double& Min, const double& Max){
   do{
	  curValue=NormalRandom(mean, sigma);
   }while(curValue<Min || curValue>Max);
   if (isInt) makeCurValueInt();
}
double TNormalPrior::getPriorDensityFromValue(const double& value){
   if(value>upperLimit || value<lowerLimit) return 0.;
   //standardize value
   double y = ( value - mean ) / sigma;
   //no calculate density according to standart formulae, where 0.3989...=1/sqrt(2*PI)
   return 0.398942280401433 * exp ( -0.5 * y * y );
}
//------------------------------------------------------------------------------
//lognormal prior
TLogNormalPrior::TLogNormalPrior(my_string Name, double min, double max, double Mean, double Sigma, bool IsInt):TNormalPrior(Name, min, max, Mean, Sigma, IsInt){
   if(lowerLimit<0.) throw TException ("Log normal prior '" + name + "' initialized with min < 0!", _FATAL_ERROR);
   //mean and sigma provided is in x-space, but we need it in the log(x) space!
   if(sigma<=0.) throw TException ("Log normal prior '" + name + "' initialized with stdev <=0!", _FATAL_ERROR);
   mean=log(Mean) - 0.5 * log(1 + ( (Sigma*Sigma)/(Mean*Mean) ));
   sigma=sqrt( log( ((Sigma*Sigma)/(Mean*Mean)) +1 ) );
   //now take care of the truncation
   inverseTruncatedArea=0;
   //lower truncation
   if(lowerLimit>0) inverseTruncatedArea+=cumulativeDistributionFunction(log(lowerLimit));
   //upper truncation
   inverseTruncatedArea+=1-cumulativeDistributionFunction(log(lowerLimit));
   inverseTruncatedArea=1/inverseTruncatedArea;
}
void TLogNormalPrior::changeCurrentValue(){
   changeCurrentValueWithLimits(lowerLimit, upperLimit);
}
double TLogNormalPrior::getPriorDensity(){
   return getPriorDensityFromValue(curValue);
}
void TLogNormalPrior::changeCurrentValueWithLimits(const double& Min, const double& Max){
   do{
	  curValue=exp(NormalRandom(mean, sigma));
   }while(curValue<Min || curValue>Max);
   if (isInt) makeCurValueInt();
}
double TLogNormalPrior::getPriorDensityFromValue(const double& value){
	if(value>upperLimit || value<lowerLimit) return 0.0;

	double y = ( log(value) - sigma ) / mean;
	return inverseTruncatedArea*((1/value)*exp ( -0.5 * y * y ));
}
//------------------------------------------------------------------------------
//gamma prior
/*
TGammaPrior::TGammaPrior(my_string Name, double min, double max, bool IsInt):TSimplePrior(Name, min, max, IsInt){
   throw TException("Gamma Prior not tested!!!!", _FATAL_ERROR);
}
void TGammaPrior::changeCurrentValue(){
   changeCurrentValueWithLimits(lowerLimit, upperLimit);
}
double TGammaPrior::getPriorDensity(){
   return getPriorDensityFromValue(curValue);
}
void TGammaPrior::changeCurrentValueWithLimits(const double& Min, const double& Max){
   do{
	  curValue=gamma_dev(firstParameter, secondParameter);
   }while(curValue<Min || curValue>Max);
   if (isInt) makeCurValueInt();
}
double TGammaPrior::getPriorDensityFromValue(const double& value){
   if(value>upperLimit || value<lowerLimit) return 0.0;
   return(DFGamma(value, secondParameter, firstParameter));
}
*/

//------------------------------------------------------------------------------




