//---------------------------------------------------------------------------

#ifndef TParametersH
#define TParametersH
#include <map>
#include "TException.h"
//---------------------------------------------------------------------------

class TParameters{
	public:
		map<my_string, my_string> mapParameter;
		map<my_string, bool> parameterUsed;
		map<my_string, my_string>::iterator curParameter, endParameter;


		TParameters(){}
		TParameters(my_string fileName, int numCommandLineParams=0, my_string* commandLineParams=NULL);
		~TParameters(){
			mapParameter.clear();
			parameterUsed.clear();
		}
		bool parameterExists(my_string my_name);
		my_string getParameter(my_string my_name, bool mandatory=1);
		double getdoubleParameter(my_string my_name, bool mandatory=1);
		my_string getListOfUnusedParameters();

};
#endif
