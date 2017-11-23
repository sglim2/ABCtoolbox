//---------------------------------------------------------------------------
#include "TException.h"
#include "TStandardEstimation.h"
//---------------------------------------------------------------------------
int main(int argc, char* argv[]){
    cout << endl << "  ABCestimator  " << endl << "***********" << endl << endl;
   try{
		if (argc <2){
		   //Show explanation if there are not enough parameters
		   cout << "You have specified " << argc << " arguments instead of 2 or more!"
				<< "\n* The following arguments are needed:"
				<< "\n*  1 :  this executable (default, normally already included)"
				<< "\n*  2 :  InputFile"
		        << "\n* (3):  optional arguments of the form NAME=VALUE\n";
           throw TException("Wrong number of arguments!", _FATAL_ERROR);
		}
		//read inputfile
		cout << "Filenames:" << endl << "----------" << endl;
		cout << "- Reading inputfile '" << argv[1] << "' ...";
		//read parameters from the command line
		TParameters* myParameters;
		if(argc>2){
			my_string* params=new my_string[argc-2];
			for(int i=2;i<argc;++i) params[i-2]=argv[i];
			myParameters = new TParameters(argv[1], argc-2, params);
		} else myParameters = new TParameters(argv[1]);
		cout << " done!" << endl;
		//create estimation object -> switch between types
		TEstimation* myEstimation;
		my_string estimationType=myParameters->getParameter("estimationType",0);
		if(estimationType=="standard" || estimationType==""){
			cout << "- Performing a standard estimation." << endl;
			myEstimation=new TStandardEstimation(myParameters);
		}
		else throw(TException("The estimation type '"+estimationType+"' does not exist!", _FATAL_ERROR));
		myEstimation->performEstimations();
		my_string unusedParams=myParameters->getListOfUnusedParameters();
		if(unusedParams!="") cout << "\nThe following parameters were not used: " << unusedParams << "!" << endl;
		//delete myParameters;
		//delete myEstimation;
    }
	catch (TException & error){
	  cout <<"\nERROR: " << error.getMessage() << endl;
	}
	// catch exceptions thrown by newmat11
	catch(BaseException) {
	   //cout << BaseException::what() << endl;
	}
	catch(std::exception & error){
		cout << error.what();
	}
	catch (...){
	  cout <<"\nERROR: unhandeled error!" << "\nProgram is finished!"  << endl;
	}
	cout << "\nEnd of Program!!\n";
    return 0;
}
//---------------------------------------------------------------------------
