//---------------------------------------------------------------------------


#include "TEstimation.h"
//---------------------------------------------------------------------------

TEstimation::TEstimation(TParameters* parameters){
   //read some parameters
   posteriorDensityPoints=parameters->getdoubleParameter("posteriorDensityPoints", false);
   if(posteriorDensityPoints==0) posteriorDensityPoints=100;
   writeSmoothedsimsUsed=parameters->getdoubleParameter("writeRetained", false);
   outputPrefix=parameters->getParameter("outputPrefix", false);
   if(outputPrefix=="") outputPrefix="ABC_GLM_";
   standardizeStats=parameters->getdoubleParameter("stadardizeStats", false);
   //read Data
   myObsData=new TObsData(parameters->getParameter("obsName"));
   mySimData=new TSimData(parameters->getParameter("simName"), parameters->getParameter("params"), parameters->getdoubleParameter("maxReadSims"), myObsData);

   if(parameters->parameterExists("numRetained")){
  	   numToRetain=parameters->getdoubleParameter("numRetained");
  	   retainedDefined=true;
  	   if(numToRetain<10) throw TException("The number of simulations to retain is set to a value <10!", _FATAL_ERROR);
   } else retainedDefined=false;
   if(!parameters->parameterExists("tolerance")) thresholdDefined=false;
   else {
  	   thresholdDefined=true;
  	   threshold=parameters->getdoubleParameter("tolerance", false);
   }

   diracPeakWidth=parameters->getdoubleParameter("diracPeakWidth", false);
   if(diracPeakWidth==0){
	   if(retainedDefined) diracPeakWidth=1/numToRetain;
	   else diracPeakWidth=0.001;
   }
   writeRetainedSimulations=parameters->getdoubleParameter("writeRetained", false);

   calcObsPValue=parameters->getdoubleParameter("obsPValue", false);

    //calculate distance from different stats?
    if(parameters->getParameter("distSimFile", false)!=""){
   	  if(parameters->getParameter("distObsFile", false)!=""){
   		 cout << endl << "Files to calculate distances:" << endl << "-----------------------------" << endl;
   		 distFromDifferentStats=true;
   		 myObsDataForDist=new TObsData(parameters->getParameter("distObsFile"));
   		 mySimDataForDistForDist=new TSimData(parameters->getParameter("distSimFile"), 0, mySimData->numReadSims, myObsDataForDist);
   		 //check
   		 if(myObsDataForDist->numObsDataSets!=myObsData->numObsDataSets)
   			throw TException("The number of observed data sets differs between the file used for estimation and the one used to calculate distances!", _FATAL_ERROR);
   	  } else {
   		 cout << "The parameter 'distObsFile' is missing in the input file" << endl
   			  << "   --> calculating distances from the stats defined in the simulation file!" << endl;
   			  distFromDifferentStats=false;
   	  }
      } else {
    	  if(parameters->getParameter("distObsFile", false)!=""){
    		  cout << "The parameter 'distSimFile' is missing in the input file" << endl
   			  << "   --> calculating distances from the stats defined in the simulation file!" << endl;
    	  }
    	  distFromDifferentStats=false;
      }

	  //open file to write marginal density
      my_string filename=outputPrefix+"marginalDensity.txt";
      cout << filename.c_str() << endl;
      //mdFile.open(filename.c_str());
      filename=outputPrefix+"L1DistancePriorPosterior.txt";
      propCommonSurfacePriorPosterior.open(filename.c_str());
      //write header
      for(int p=0; p<mySimData->numParams;++p){
    	  propCommonSurfacePriorPosterior << "\t" << mySimData->paramNames[p] << "\tdensity";
      }
}
//---------------------------------------------------------------------------
void TEstimation::standardize(){
	//standardize params
	cout << "Standardizing parameters ...";
    mySimData->standardizeParameters();
    cout << "done!" << endl;
    //standaridze stats
    if(standardizeStats){
    	cout << "Standardizing statistics ...";
    	mySimData->standardizeStatistics();
    	myObsData->standardizeObservedValues(mySimData->statMeans, mySimData->statSDs);
    	cout << " done!" << endl;
    }
}
//---------------------------------------------------------------------------
void TEstimation::preparePosteriorDensityPoints(){
	 //prepare vector for points at which the density (also prior!!!) is estimated
	   theta=ColumnVector(posteriorDensityPoints);
	   //float step=1.2/(posteriorDensityPoints-1);
	   //for(int j=1; j<=posteriorDensityPoints; ++j) theta(j)=-0.1+(j-1)*step;
	   //since the parameters are standardized between 0 and 1, we can easy calculate the posterior between 0 and 1
	   step=1.0/(posteriorDensityPoints-1.0);
	   for(int j=0; j<posteriorDensityPoints; ++j) theta(j+1)=(j)*step;
}
//---------------------------------------------------------------------------
void TEstimation::performEstimations(){
   cout << endl << "No estimation procedure defined for the base class!" << endl;
}
//---------------------------------------------------------------------------
void TEstimation::preparePrior(){
	cout << "   - using smoothed empirical distribution of parameters as truncated prior ...";
	priorMatrix=Matrix(mySimData->numUsedSims, mySimData->numParams+1);
	priorMatrix=mySimData->paramMatrix;
	cout << "done!" << endl;
}
//---------------------------------------------------------------------------
void TEstimation::calculatePosterior(){
	//prepare posterior Matrix
	posteriorMatrix=Matrix(posteriorDensityPoints, mySimData->numParams);
	posteriorMatrix=0.0;
	Matrix T_j=calculate_T_j();
	ColumnVector theta_j, v_j;
	//now loop over all parameters
	cout << "   - calculating posterior densities for " << mySimData->numParams << " parameters ...";
	//first calculate all Qj in order to shift them -> otherwise maybe too large numbers...
	//the discrepancy is constant and therefore integrated out by normArea
	//loop over all dirac peaks
	double* exponents=new double[mySimData->numUsedSims];
	double meanExponent=0;
	for(int d=1; d<=mySimData->numUsedSims;++d){
		//theta_j is a vector with the parameters of one simulation / from the uniform matrix
		theta_j=priorMatrix.submatrix(d,d,2,mySimData->numParams+1).as_column();
		v_j=calculate_v_j(theta_j);
		//removed (obsValues-c_zero).t()*SigmaSInv*(obsValues-c_zero) from the exponent since it is constant -> corrected by normArea
		Matrix exponent=theta_j.t()*SigmaThetaInv*theta_j-v_j.t()*T_j*v_j;
		exponents[d-1]=-0.5*exponent(1,1);
		meanExponent+=exponents[d-1];
	}

	//shift the Qj
	meanExponent=meanExponent/mySimData->numUsedSims;
	for(int d=0; d<mySimData->numUsedSims;++d) exponents[d]=exponents[d]-meanExponent;
	//now calculate the posterior densities
	for(int p=1; p<=mySimData->numParams;++p){
		float tau=T_j(p,p);
		//loop over all dirac peaks
		for(int d=1; d<=mySimData->numUsedSims;++d){
			//theta_j is a vector with the parameters of one simulation / from the uniform matrix
			theta_j=priorMatrix.submatrix(d,d,2,mySimData->numParams+1).as_column();
			v_j=calculate_v_j(theta_j);
			double c_j=exp(exponents[d-1]);
			ColumnVector t_j=T_j*v_j;

			for(int k=1;k<=posteriorDensityPoints;++k){
				posteriorMatrix(k,p)=posteriorMatrix(k,p)+c_j*exp(-pow(theta(k)-t_j(p),2.)/(2.*tau));
			}
		}
	}
	cout << " done!" << endl;
}
//---------------------------------------------------------------------------
Matrix TEstimation::calculate_T_j(){
	   return (CHat.t()*SigmaSInv*CHat+SigmaThetaInv).i();
}
//---------------------------------------------------------------------------
ColumnVector TEstimation::calculate_v_j(ColumnVector theta_j){
	return CHat.t()*SigmaSInv*(obsValues-c_zero)+SigmaThetaInv*theta_j;
}
//---------------------------------------------------------------------------
void TEstimation::performLinearRegression(){
   //perform regression
   cout << "   - performing local linear regression ...";
   C=mySimData->paramMatrix.t()*mySimData->paramMatrix;
   C=C.i();
   C=mySimData->paramMatrix*C;
   C=mySimData->statMatrix*C;
   c_zero=C.column(1);
   CHat = C.columns(2,C.ncols());
   RHat = mySimData->statMatrix.t()-mySimData->paramMatrix*C.t();
   SigmaS=1.0/(mySimData->numUsedSims-mySimData->numParams)*RHat.t()*RHat;
   SigmaSInv=SigmaS.i();
   cout << " done!" << endl;
}
//---------------------------------------------------------------------------
void TEstimation::prepareDiracPeaks(){
	//approximation of priors with Dirac peaks
	cout << "   - approximation of the prior with Dirac peaks ...";
	DiagonalMatrix diracPeakWidthMatrix(mySimData->numParams);
	diracPeakWidthMatrix=diracPeakWidth;
	//SigmaTheta=pow(1./numToRetain, 2./mySimData->numParams)*diracPeakWidthMatrix;
	SigmaTheta=diracPeakWidthMatrix;
	SigmaThetaInv=SigmaTheta.i();
	cout << " done!" << endl;
}
//---------------------------------------------------------------------------
float TEstimation::theatParameterScale(int param, int k){
	return mySimData->paramMinima[param-1]+theta(k)*(mySimData->paramMaxima[param-1]-mySimData->paramMinima[param-1]);
}
//---------------------------------------------------------------------------
void TEstimation::writePosteriorFile(my_string filenameTag){
   //write file with posterior
   cout << "   - write file with posterior estimates ...";
   ofstream output;
   my_string filename=outputPrefix+"PosteriorEstimates"+filenameTag+".txt";
   output.open(filename.c_str());
   //write Header
   output << "number";
   for(int p=0; p<mySimData->numParams;++p){
	  output << "\t" << mySimData->paramNames[p] << "\tdensity";
   }
   output << endl;
   //write values for each parameter: norm area!
   ColumnVector area(mySimData->numParams);
   for(int p=1; p<=mySimData->numParams;++p){
	  area(p)=getArea(posteriorMatrix.column(p), theta)*(mySimData->paramMaxima[p-1]-mySimData->paramMinima[p-1]);
   }
   for(int k=1;k<=posteriorDensityPoints;++k){
	  output << k;
	  for(int p=1; p<=mySimData->numParams;++p){
		 output << "\t" << theatParameterScale(p, k) << "\t" << posteriorMatrix(k,p)/area(p);
	  }
	  output << endl;
   }
   output.close();
   cout << " done!" << endl;
}
//---------------------------------------------------------------------------
float TEstimation::getPositionofMode(int param){
	float max=0;
	float pos;
	for(int k=1;k<=posteriorDensityPoints;++k){
		if(posteriorMatrix(k,param)>max){
			max=posteriorMatrix(k,param);
			pos=k;
		}
	}
	return theatParameterScale(param, pos);
}
//---------------------------------------------------------------------------
float TEstimation::getPosteriorMean(int param){
	//norm area!!
	float areaPosterior=getArea(posteriorMatrix.column(param), theta)*(mySimData->paramMaxima[param-1]-mySimData->paramMinima[param-1]);
	float mean=0;
	for(int k=1;k<=posteriorDensityPoints;++k){
			mean+=posteriorMatrix(k,param)*theatParameterScale(param, k);
	}
	return mean*step/areaPosterior;
}
//---------------------------------------------------------------------------
float TEstimation::getPosteriorquantile(int param, float q){
	//norm area!!
	float areaPosterior=getArea(posteriorMatrix.column(param), theta);
	float area=0;
	float a;
	for(int k=1;k<=posteriorDensityPoints;++k){
		a=step*posteriorMatrix(k,param)/areaPosterior;
		if((area+a)>q){
			float add=step*(q-area)/a - step/2;
			add=add*(mySimData->paramMaxima[param-1]-mySimData->paramMinima[param-1]);
			return theatParameterScale(param, k)+add;
		}
		else area+=a;
	}
	return 0;
}
//---------------------------------------------------------------------------
void TEstimation::getPosteriorHDI(int param, float q, float* lower, float* upper){
	//get shortest interval of fraction q
	//I assume an unimodal distribution, therefore, just adding bins from the largest to the smallest until q is reached
	//since two values have to be returned, the results are written to a passed variable
	ColumnVector d=posteriorMatrix.column(param).AsColumn();
	SortDescending(d);
	float areaPosterior=getArea(posteriorMatrix.column(param), theta);
	float area=0;
	float criticalDensity=0;
	for(int k=1;k<=posteriorDensityPoints;++k){
		criticalDensity=d(k);
		area+=step*d(k)/areaPosterior;
		if(area>q) break;
	}

	//now find coordinates
	float minK=0;
	for(int k=1;k<=posteriorDensityPoints;++k){
		if(posteriorMatrix(k,param)>criticalDensity){
			minK=k;
			break;
		}
	}
	float maxK=0;
	for(int k=posteriorDensityPoints;k>0;--k){
		if(posteriorMatrix(k,param)>criticalDensity){
			maxK=k;
			break;
		}
	}
	//now get area within interval
	if(minK==maxK){
		cout << "\nWARNING: Problems calculating the HDI for parameter '" + mySimData->paramNames[param-1] + "'!";
		*lower=0;
		*upper=0;
	} else {
		area=0;
		for(int k=(minK+1);k<maxK;++k) area+=step*posteriorMatrix(k,param);
		area=area/areaPosterior;
		float dif=area-q;
		//how much to add?
		float missing=dif/((posteriorMatrix(minK,param)+posteriorMatrix(maxK,param))/areaPosterior);
		*lower=theatParameterScale(param, minK)-missing*step;
		*upper=theatParameterScale(param, maxK)+missing*step;
	}
}
//---------------------------------------------------------------------------
void TEstimation::writePosteriorCharacteristics(my_string filenameTag){
	   //write file with posterior characteristics, such as mean, mode etc.
	   cout << "   - write file with posterior characteristics ...";

	   ofstream output;
	   my_string filename=outputPrefix+"PosteriorCharacteristics"+filenameTag+".txt";
	   output.open(filename.c_str());
	   //write Header
	   output << "what";
	   for(int p=1; p<=mySimData->numParams;++p){
		  output << "\t" << mySimData->paramNames[p-1];
	   }
	   output << endl;

	   //get mode of the posterior
	   output << "mode";
	   for(int p=1; p<=mySimData->numParams;++p){
		   output << "\t" << getPositionofMode(p);
	   }
	   output << endl;
	   //get mean of the posterior
	   output << "mean";
	   for(int p=1; p<=mySimData->numParams;++p){
   		   output << "\t" << getPosteriorMean(p);
   	   }
   	   output << endl;
   	   //get median of the posterior
	   output << "median";
	   for(int p=1; p<=mySimData->numParams;++p){
   		   output << "\t" << getPosteriorquantile(p, 0.5);
   	   }
   	   output << endl;


   	   //get Quantiles  of the posterior
   	   output << "quantile_50_lower_bound";
   	   for(int p=1; p<=mySimData->numParams;++p){
   		   output << "\t" << getPosteriorquantile(p, 0.25);
   	   }
   	   output << endl;
   	   output << "quantile_50_upper_bound";
   	   for(int p=1; p<=mySimData->numParams;++p){
   		   output << "\t" << getPosteriorquantile(p, 0.75);
   	   }
   	   output << endl;
   	   output << "quantile_90_lower_bound";
   	   for(int p=1; p<=mySimData->numParams;++p){
   		   output << "\t" << getPosteriorquantile(p, 0.05);
   	   }
   	   output << endl;
   	   output << "quantile_90_upper_bound";
   	   for(int p=1; p<=mySimData->numParams;++p){
   		   output << "\t" << getPosteriorquantile(p, 0.95);
   	   }
   	   output << endl;
   	   output << "quantile_95_lower_bound";
   	   for(int p=1; p<=mySimData->numParams;++p){
   		   output << "\t" << getPosteriorquantile(p, 0.025);
   	   }
   	   output << endl;
   	   output << "quantile_95_upper_bound";
   	   for(int p=1; p<=mySimData->numParams;++p){
   		   output << "\t" << getPosteriorquantile(p, 0.975);
   	   }
   	   output << endl;
   	   output << "quantile_99_lower_bound";
   	   for(int p=1; p<=mySimData->numParams;++p){
   		   output << "\t" << getPosteriorquantile(p, 0.005);
   	   }
   	   output << endl;
   	   output << "quantile_99_upper_bound";
   	   for(int p=1; p<=mySimData->numParams;++p){
   		   output << "\t" << getPosteriorquantile(p, 0.995);
   	   }
   	   output << endl;

	   //get HDI
	   float lower, upper;
	   float* upperStore=new float[mySimData->numParams];
	   output << "HPD_50_lower_bound";
	   for(int p=1; p<=mySimData->numParams;++p){
		   getPosteriorHDI(p, 0.5, &lower, &upper);
		   upperStore[p-1]=upper;
		   output << "\t" << lower;
	   }
	   output << endl;
	   output << "HPD_50_upper_bound";
	   for(int p=1; p<=mySimData->numParams;++p){
		   output << "\t" << upperStore[p-1];
	   }
	   output << endl;
	   output << "HPD_90_lower_bound";
	   for(int p=1; p<=mySimData->numParams;++p){
		   getPosteriorHDI(p,0.9, &lower, &upper);
		   upperStore[p-1]=upper;
		   output << "\t" << lower;
	   }
	   output << endl;
	   output << "HPD_90_upper_bound";
	   for(int p=1; p<=mySimData->numParams;++p){
		   output << "\t" << upperStore[p-1];
	   }
	   output << endl;
	   output << "HPD_95_lower_bound";
	   for(int p=1; p<=mySimData->numParams;++p){
		   getPosteriorHDI(p,0.95, &lower, &upper);
		   upperStore[p-1]=upper;
		   output << "\t" << lower;
	   }
	   output << endl;
	   output << "HPD_95_upper_bound";
	   for(int p=1; p<=mySimData->numParams;++p){
		   output << "\t" << upperStore[p-1];
	   }
	   output << endl;
	   output << "HPD_99_lower_bound";
	   for(int p=1; p<=mySimData->numParams;++p){
		   getPosteriorHDI(p,0.99, &lower, &upper);
		   upperStore[p-1]=upper;
		   output << "\t" << lower;
	   }
	   output << endl;
	   output << "HPD_99_upper_bound";
	   for(int p=1; p<=mySimData->numParams;++p){
		   output << "\t" << upperStore[p-1];
	   }
	   output << endl;

	   //done!

	   output.close();
	   cout << " done!" << endl;
}
//---------------------------------------------------------------------------
void TEstimation::writeSmoothedSimsUsed(my_string filenameTag){
   //write file with prior
   cout << "   - write file with prior densities ...";
   ofstream output;
   my_string filename=outputPrefix + "TruncatedPrior";
   filename+= filenameTag;
   filename+=".txt";
   output.open(filename.c_str());
   //write Header
   output << "number";
   for(int p=0; p<mySimData->numParams;++p){
	  output << "\t" << mySimData->paramNames[p] << "\tdensity";
   }
   output << endl;
   //write values for each parameter: norm area!
   ColumnVector area(mySimData->numParams);
   for(int p=1; p<=mySimData->numParams;++p){
	   area(p)=getArea(retainedMatrix.column(p), theta)*(mySimData->paramMaxima[p-1]-mySimData->paramMinima[p-1]);
   }
   for(int k=1;k<=posteriorDensityPoints;++k){
	  output << k;
	  for(int p=1; p<=mySimData->numParams;++p){
		 output << "\t" << theatParameterScale(p, k) << "\t" << retainedMatrix(k,p)/area(p);
	  }
	  output << endl;
   }
   output.close();
   cout << " done!" << endl;
}
//---------------------------------------------------------------------------
float TEstimation::getCommonSurfaceRetainedPosterior(int p){
	//calculate common surface: norm area
	float areaPrior;
	float areaPosterior;
	float commonSurface=0;
	float posteriorSurface=0;
	areaPrior=getArea(retainedMatrix.column(p), theta)*(mySimData->paramMaxima[p-1]-mySimData->paramMinima[p-1]);
	areaPosterior=getArea(posteriorMatrix.column(p), theta)*(mySimData->paramMaxima[p-1]-mySimData->paramMinima[p-1]);
	for(int k=1;k<=posteriorDensityPoints;++k){
		commonSurface+=min(retainedMatrix(k,p)/areaPrior, posteriorMatrix(k,p)/areaPosterior);
		posteriorSurface+=posteriorMatrix(k,p)/areaPosterior;
	}
	return commonSurface/posteriorSurface;
}
//---------------------------------------------------------------------------
double TEstimation::getArea(ColumnVector densities, ColumnVector theta){
   double area=0;
   for(int i=1; i<densities.size();++i){
	  area+=(theta(i+1)-theta(i))*(densities(i+1)+densities(i))/2;
   }
   return area;
}
//---------------------------------------------------------------------------
void TEstimation::calculateSmoothedRetainedMatrix(){
	//prepare retained Matrix -> same points as for the posterior!
	retainedMatrix=Matrix(posteriorDensityPoints, mySimData->numParams);
	retainedMatrix=0.0;
	//calculate the retained densities
	ColumnVector theta_j;
	for(int p=1; p<=mySimData->numParams;++p){
		//loop over all dirac peaks
		for(int d=1; d<=mySimData->numUsedSims;++d){
			//theta_j is a vector with the parameters of one simulation / from the uniform matrix
			theta_j=priorMatrix.submatrix(d,d,2,mySimData->numParams+1).as_column();
			for(int k=1;k<=posteriorDensityPoints;++k){
				retainedMatrix(k,p)=retainedMatrix(k,p)+exp(-pow(theta(k)-theta_j(p),2.)/(2.*diracPeakWidth));
			}
		}
	}
}
//---------------------------------------------------------------------------
void TEstimation::calculateMarginalDensitiesOfRetained(){
   //ofstream fmfile;
   //fmfile.open("fmfile.txt");
	if(calcObsPValue>mySimData->numUsedSims) calcObsPValue=mySimData->numUsedSims;
	cout << "   - calculating the marginal density for " << calcObsPValue << " used sims to get the p-value of the observed data." << endl;
   fmFromRetained=new float[mySimData->numUsedSims];
   Matrix D = SigmaS+CHat*SigmaTheta*CHat.t();
   Matrix D_inv=D.i();
   ColumnVector theta_j, m_j, simStats;
   for(int s=1; s<=calcObsPValue;++s){
	   simStats=mySimData->statMatrix.column(s);
	   double f_M=0.0;
	   //loop over all dirac peaks
	   for(int d=1; d<=mySimData->numUsedSims;++d){
		   theta_j=mySimData->paramMatrix.submatrix(d,d,2,mySimData->numParams+1).as_column();
		   m_j=c_zero+CHat*theta_j;
		   Matrix temp=-0.5*(simStats-m_j).t()*D_inv*(simStats-m_j);
		   f_M+=exp(temp(1,1));
	   }
	   f_M=f_M/(mySimData->numUsedSims*sqrt((2*3.14159265358979323846*D).determinant()));
	   fmFromRetained[s-1]=f_M*((double)mySimData->numUsedSims/(double)mySimData->numReadSims);
	  // fmfile << fmFromRetained[s-1] << endl;
   }
   //fmfile.close();
}
//---------------------------------------------------------------------------

