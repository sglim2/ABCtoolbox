/*
 * glm.cpp
 *
 *  Created on: Mar 2, 2009
 *      Author: wegmannd
 */

//---------------------------------------------------------------------------
#include "TException.h"
#include <time.h>
#include "my_cstring.h"
#include <fstream>
#include <iostream>
#include <vector>
#include "newmat.h"
#include "newmatap.h"
#include "newmatio.h"
//---------------------------------------------------------------------------------------
long get_randomSeedFromCurrentTime(){
   time_t seconds;
   seconds = time (NULL);
   return clock()+seconds - 37*365*24*3600; //substract 37 years...
}
//---------------------------------------------------------------------------------------
#define MBIG 1000000000L
#define MSEED 161803398L
#define MZ 0
#define FAC (1.0/MBIG)

double ran3(long *idum)
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;

	if ((*idum < 0) || (iff == 0) ) {
		iff=1;
		mj=MSEED-(*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i<=54;++i) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1;k<=4;++k)
			for (i=1;i<=55;++i) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		inext=0;
		inextp=31;
		*idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;
	return mj*FAC;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
//---------------------------------------------------------------------------------------
void initialize_ran3(long& _Idum){
	_Idum=get_randomSeedFromCurrentTime();
	_Idum=_Idum * (101 - _Idum % 101) / 50;
	if(_Idum>0) _Idum=-_Idum;
	if(_Idum<-161803398) _Idum=_Idum/100;
	if(_Idum<-161803398) _Idum=_Idum/100;
	if(_Idum<-161803398) _Idum=892374;
	//cout << "Random generator initialized with seed = " << -_Idum << endl;
	ran3(&_Idum);
  	return;
}
//---------------------------------------------------------------------------
/* Returns a Normal random variate based on a unit variate,
   using a random generator as a source of uniform deviates.
   Adapted from the algorithm described in the book Numerical Recipes by
   Press et al.
*/
double NormalRandom (const double& dMean, const double& dStdDev, long& lidum){
  double w, x1, x2, dTemp3;

   do {
      x1 = 2. * ran3(&lidum) - 1.;
      x2 = 2. * ran3(&lidum) - 1.;
      w = x1*x1 + x2*x2;
   } while (w >= 1. || w < 1E-30);

   w = sqrt((-2.*log(w))/w);
   x1 *= w;
   return (x1 * dStdDev + dMean);
} /* NormalRandom */
//---------------------------------------------------------------------------

void showExplanations(){
	cout << "\n***************************************************************************************"
	     << "\n* This program calculates n summary statistics S under a general linear model (GLM)   *"
	     << "\n* of the form S=C*P+c0+e where P is a vector of m parameter values and e is the error *"
	     << "\n* term. It is designed to test simCoalABC or related program since it can be used as  *"
	     << "\n* a simulation program. The following arguments are needed:                           *"
		 << "\n*           1 : Name of the file with the definitions of the matrices, the c0 vector  *"
		 << "\n*               and the variance used for the normal error distribution e. This file  *"
		 << "\n*               simply starts with the tag \"C\" followed by the matrix C beginning on  *"
		 << "\n*               a new line and containing only numbers organized in m columns on n    *"
		 << "\n*               lines. After a line with the tag \"c0\" a single line with the n values *"
		 << "\n*               of the c0 vector. Finally, after a line with the tag \"Sigma\", n lines *"
		 << "\n*               with n components of the Variance-Covariance matrix Sigma should be   *"
		 << "\n                given. Note that Sigma must be a symmetric matrix."
		 << "\n*           2 : The name of the output file (e.g. \"summary_statistics.temp\")          *"
		 << "\n*           3 : The value of the first parameter                                      *"
		 << "\n*   (3 to m+2): The values of the remaining parameters                                *"
		 << "\n***************************************************************************************\n\n";
}

int main(int argc, char* argv[]){
	try{
		if (argc<4) throw TException("You have specified " + (my_string) (argc-1) + " parameters instead of 3 or more!", _FATAL_ERROR);

		//read file with matrix definitions
		ifstream is;
		my_string filename=argv[1];
		is.open(filename.c_str());
		if(!is) throw TException("The file with the definitions of the matrices '"+filename+"' could not be opend!", _FATAL_ERROR);
		my_string temp, buf;
		buf.read_line(is);
		buf.trim_blanks();
		if(buf!="C") throw TException("The file with the definitions of the matrices '"+filename+"' does not start with matrix C (tag missing?)!", _FATAL_ERROR);

		//read matrix C as a vector of arrays
		vector<double> firstLine;
		vector<double*> Ctemp;
		int numParams;
		//read first line
		buf.read_line(is);
		buf.trim_blanks();
		while(!buf.empty()){
			temp=buf.extract_sub_str_before_ws();
			temp.trim_blanks();
			firstLine.push_back(temp.toDouble());
		}
		numParams=firstLine.size();
		Ctemp.push_back(new double[numParams]);
		for(int i=0;i<numParams;++i) Ctemp[0][i]=firstLine[i];
		//read all the other lines
		buf.read_line(is);
		while(!buf.contains("c0")){
			buf.trim_blanks();
			if(!buf.empty()){
				Ctemp.push_back(new double[numParams]);
				int i=0;
				while(!buf.empty()){
					if(i==numParams) throw TException("The file with the definitions of the matrices '"+filename+"' contains unequal number of values on the different lines specifying the matrix C!", _FATAL_ERROR);
					temp=buf.extract_sub_str_before_ws();
					temp.trim_blanks();
					Ctemp[Ctemp.size()-1][i]=temp.toDouble();
					++i;
				}
				if(i<numParams) throw TException("The file with the definitions of the matrices '"+filename+"' contains unequal number of values on the different lines specifying the matrix C!", _FATAL_ERROR);
			}
			buf.read_line(is);
		}
		int numStats=Ctemp.size();

		Matrix C(numStats, numParams);
		for(int i=0; i<numStats; ++i){
			for(int j=0; j<numParams; ++j){
				C.element(i,j)=Ctemp[i][j];
			}
		}

		//read the c0 vector
		ColumnVector c0(numStats);
		buf.read_line(is);
		buf.trim_blanks();
		int i=0;
		while(!buf.empty()){
			if(i==numStats) throw TException("The file with the definitions of the matrices '"+filename+"' contains too many values for c0!", _FATAL_ERROR);
			temp=buf.extract_sub_str_before_ws();
			temp.trim_blanks();
			c0.element(i)=temp.toDouble();
			++i;
		}
		if(i<numStats) throw TException("The file with the definitions of the matrices '"+filename+"' contains too few values for c0!", _FATAL_ERROR);


		//read variances
		buf.read_line(is);
		buf.trim_blanks();
		if(buf!="Sigma") throw TException("The file with the definitions of the matrices '"+filename+"' does not contain the matrix Sigma (tag missing?)!", _FATAL_ERROR);
		SymmetricMatrix Sigma(numStats);


		for(int line=0; line<numStats;++line){
			buf.read_line(is);
			buf.trim_blanks();
			if(buf.empty()) throw TException("The file with the definitions of the matrices '"+filename+"' contains too few rows for Sigma!", _FATAL_ERROR);
			i=0;
			while(!buf.empty()){
				if(i==numStats) throw TException("The file with the definitions of the matrices '"+filename+"' contains too many values for Sigma on line "+line+"!", _FATAL_ERROR);
				temp=buf.extract_sub_str_before_ws();
				temp.trim_blanks();
				Sigma.element(line, i)=temp.toDouble();
				++i;
			}
			if(i<numStats) throw TException("The file with the definitions of the matrices '"+filename+"' contains too few values for Sigma on line "+line+"!", _FATAL_ERROR);
		}
		buf.read_line(is);
		buf.trim_blanks();
		if(!buf.empty()) throw TException("The file with the definitions of the matrices '"+filename+"' contains too few rows for Sigma!!", _FATAL_ERROR);
		is.close();

		//read params from input
		if(argc<numParams+3) throw TException("Too few parameter values passed!", _FATAL_ERROR);
		if(argc>numParams+3) throw TException("Too many parameter values passed!", _FATAL_ERROR);
		ColumnVector P(numParams);
		for(i=0;i<numParams;++i){
			temp=argv[i+3];
			P.element(i)=temp.toDouble();
		}

		//write header to output file
		ofstream out;
		filename=argv[2];
		out.open(filename.c_str());
		if(!out) throw TException("The output file '"+filename+"' could not be opened!", _FATAL_ERROR);
		out << "Stat_1";
		for(i=1;i<numStats;++i){
			out << "\t" << "Stat_" << i+1;
		}
		out << endl;

		//initialize random generator
		long _Idum=1.0;
		initialize_ran3(_Idum);

		//prepare epsilon
		Matrix A;
		try{
			A=Cholesky(Sigma);
		} catch (...){
			throw TException("problem solving the Cholesky decomposition of Sigma!", _FATAL_ERROR);
		}
		ColumnVector e(numStats);
		for(int i=0;i<numStats;++i) e.element(i)=NormalRandom(0,1, _Idum);
		e=A*e;

		//compute stats
		ColumnVector s(numStats);
		s=C*P+c0+e;

		//write stats
		for(i=0;i<numStats;++i){
			out << s.element(i);
		}
		out << endl;
		out.close();

	 }
	 catch (TException & error){
		 cout << "ERROR: " << error.getMessage() << endl;
		 showExplanations();
	 }
	 catch (...){
		 cout <<"ERROR: unhandeled error!" << endl;
		 showExplanations();
	 }
	 return 1;
}
