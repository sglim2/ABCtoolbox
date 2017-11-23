//---------------------------------------------------------------------------

#ifndef TLinearCombBoxCoxH
#define TLinearCombBoxCoxH
//---------------------------------------------------------------------------
#include <strstream>
#include <vector>
#include "my_cstring.h"
#include "TLinearComb.h"
#include <math.h>

class TLinearCombBoxCox:public TLinearComb{
	public:
	   double* lambda;
	   double* gm;
	   double* min;
	   double* max;

	   TLinearCombBoxCox();
	   TLinearCombBoxCox(my_string gotLinearCombFileName, vector<my_string> NamesInInputFiles);
	   virtual void readLinearCombFile();
	   virtual void calcObsLineraComb(vector<double> obsData);
	   virtual void calcSimDataPCA(double** simData);
	   virtual void calcSimDataPCA(double* simData);
	   double getBoxCox(double simData, int stat);
};
#endif
