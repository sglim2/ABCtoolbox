//---------------------------------------------------------------------------
#include "TParameters.h"
#include <strstream>
//---------------------------------------------------------------------------
TParameters::TParameters(my_string fileName, int numCommandLineParams, my_string* commandLineParams){
	ifstream is (fileName.c_str());
	if(!is) throw TException("Input file '" + fileName + "' could not be opened!", _FATAL_ERROR);
	my_string buf, my_name;
	my_string my_value;
	while(is.good() && !is.eof()){
		strstream isLine;                   // create the new stream
		buf.read_line(is);                  // read the line
		buf=buf.extract_before_doubleSlash();   // remove the commentated part
		if(!buf.empty()){
		// allocate the line to a new stream for easy reading
           isLine << buf;
			my_name.read_to_delim(isLine);
			isLine >> ws;
			my_value.read_to_delim(isLine);
			if(!my_name.empty() && !my_value.empty()){
				mapParameter[my_name]= my_value;
			}
		}
	}
	//now check passed parameters and overwrite earlier definitions (if available) or add new parameters
	for(int i=0;i<numCommandLineParams;++i){
		my_name=commandLineParams[i].extract_before('=');
		//if(parameterExists(my_name)) mapParameter[my_name]= commandLineParams[i].extract_after('=');
		mapParameter[my_name]= commandLineParams[i].extract_after('=');
	}
	//prepare map to store if a parameter was used
	curParameter=mapParameter.begin();
	endParameter=mapParameter.end();
	for(;curParameter!=endParameter;++curParameter){
		parameterUsed[curParameter->first]=false;
	}
	endParameter=mapParameter.end();
}
//---------------------------------------------------------------------------
bool TParameters::parameterExists(my_string my_name){
	curParameter=mapParameter.begin();
	for(;curParameter!=endParameter;++curParameter){
		if(curParameter->first==my_name){
			parameterUsed[curParameter->first]=true;
			return true;
		}
	}
    return false;
}
//---------------------------------------------------------------------------
my_string TParameters::getParameter(my_string my_name, bool mandatory){
	curParameter=mapParameter.begin();
	for(;curParameter!=endParameter;++curParameter){
		if(curParameter->first==my_name){
			parameterUsed[curParameter->first]=true;
			return curParameter->second;
		}
	}
	if(mandatory)
		throw TException("The parameter '" + my_name + "' is not defined in the inputfile!", _FATAL_ERROR);
   else return "";
}
//---------------------------------------------------------------------------
double TParameters::getdoubleParameter(my_string my_name, bool mandatory){
	curParameter=mapParameter.begin();
	for(;curParameter!=endParameter;++curParameter){
		if(curParameter->first==my_name){
			parameterUsed[curParameter->first]=true;
			return curParameter->second.toDouble();
		}
	}
	if(mandatory){
		throw TException("The parameter '" + my_name + "' is not defined in the inputfile!", _FATAL_ERROR);
   } else {
      return false;
   }
}
//---------------------------------------------------------------------------
my_string TParameters::getListOfUnusedParameters(){
	my_string parameterList="";
	map<my_string, bool>::iterator cur;
	cur=parameterUsed.begin();
	curParameter=mapParameter.begin();
	for(;cur!=parameterUsed.end() && curParameter!=endParameter;++cur, ++curParameter){
		if(!cur->second){
			if(parameterList!="") parameterList+=", ";
			parameterList+=curParameter->first;
		}
	}
	return parameterList;
}
