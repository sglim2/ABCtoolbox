//---------------------------------------------------------------------------
#pragma hdrstop
#include <unistd.h>
#include <strstream>
#include <vector>
#include "my_cstring.h"
#include "TLinearComb.h"
#include "TLinearCombBoxCox.h"
#include "TException.h"
//---------------------------------------------------------------------------

#pragma argsused

//Show explanation if there are not enought parameters
void showExplenations(){
			 cout << "\n* The following parameters are needed:"
				  << "\n*  1 :  name of the file in which the linear combinations are defined"
				  << "\n*  2 :  input filename"
				  << "\n*  3 :  output filename"
				  << "\n* (4):  if a Box-Cox transformation is requested add the tag \"boxcox\" here"
				  << "\n\nThe Linear-Combination-File is a file containing a table with rows for each statistic and teh follwoing columns:"
				  << "\n\tStatistics Name"
				  << "\n\tmean of the statistic"
				  << "\n\tsd  of the statistic"
				  << "\n\tcontribution to Linear-Combination_1"
				  << "\n\t..."
				  << "\n\tcontribution to Linear-Combination_X"
				  << "\nAll lines not containing stats have to be commented by //."
				  << "\n\nIf a Box-cox transformation is required the Linear-Combination-File contains these columns:"
				  << "\n\tStatistics Name"
				  << "\n\tmax of the statistic"
				  << "\n\tmin of the statistic"
				  << "\n\tlambda  for this statistic"
				  << "\n\tgeometric mean for this statistic"
  				  << "\n\tmean of the statistic after Box-Cox transformation"
  				  << "\n\tsd  of the statistic after Box-Cox transformation"
  				  << "\n\tcontribution to Linear-Combination_1"
  				  << "\n\t..."
  				  << "\n\tcontribution to Linear-Combination_X"
				  << "\n\nThe Output will consist of all columns prior to the first column containing a summary statistics mentioned in the Linear-Combination-File and the X specified Linear-Combinations\n"
				  << "\nCheck the manual of the ABC-package for more information\n\n";
};
//---------------------------------------------------------------------------
int main(int argc, char* argv[]){
	my_string buf;
	my_string LinearComb_fileName, inputFileName, outputFileName;
	int numParams, numValuesPerLine;
	int i;
	bool boxcox=false;

	try{
		if (argc!=4 && argc!=5){
			throw TException("You have specified "+ (my_string) argc + " parameters instead of 4!\n", _FATAL_ERROR);
		}
		//read parameters
		LinearComb_fileName=argv[1];
		LinearComb_fileName.trim_blanks();
		inputFileName=argv[2];
		inputFileName.trim_blanks();
		outputFileName=argv[3];
		outputFileName.trim_blanks();
		if(argc==5){
			my_string temp=argv[4];
			temp.trim_blanks();
			if(temp=="boxcox") boxcox=true;
			else throw TException("Unknown tag '" + temp + "'!\n", _FATAL_ERROR);
		}

		//open file to perfom the linear combinations on
		ifstream inputfile(inputFileName.c_str());
		if(!inputfile){
		   throw TException("The input-File '" + inputFileName + "' can not be read!\n", _FATAL_ERROR);
		}
		//read header line with the names
		vector<my_string> InputFileNamesVector;
		strstream isLine;
		my_string name;
		buf.read_line(inputfile);
		if(!buf.empty()){
		   isLine << buf;
		   name.read_to_delim(isLine);
		   while(!name.empty()){
			  InputFileNamesVector.push_back(name);
			  name.read_to_delim(isLine);
		   }
		}
		numValuesPerLine=InputFileNamesVector.size();

		//read Linear-Combination File
		TLinearComb* myLC;
		if(boxcox) myLC = new TLinearCombBoxCox(LinearComb_fileName, InputFileNamesVector);
		else myLC = new TLinearComb(LinearComb_fileName, InputFileNamesVector);

		//check stats
		my_string missingStat=myLC->getFirstMissingStat();
		if(!missingStat.empty()) throw TException("The in the Linear-Combination-File '" + LinearComb_fileName + "' speciefied stat '" + missingStat +"' is missing in the Inputfile '" + inputFileName + "'!\n", _FATAL_ERROR);

		//create Output-File
		ofstream outputfile(outputFileName.c_str());
		if(!outputfile){
		   throw TException("The output-File '" + outputFileName + "' can not be created!\n", _FATAL_ERROR);
		}

		//get the number of coulmns of parameters
		numParams=myLC->getNumberOfFirstUsedStat();

		//calculate LinearCombinations
		//first write the header line: parameters and then the stats
		for(int i=0; i<numParams; ++i){
			if(i>0) outputfile << "\t";
			outputfile << InputFileNamesVector[i];
		}
		myLC->writeHeader(outputfile);
		outputfile << endl;

		//now go through the input file an copy the lines
		vector<double> tempInput;
		vector<double>::iterator curTempInput;
		double* simData=new double[numValuesPerLine];
		long lineIterator=1;
		long copiedLines=0;
		while(inputfile.good() && !inputfile.eof()){
			strstream testLine;
			tempInput.clear();
			buf.read_line(inputfile);
			++lineIterator;
			if(!buf.empty()){
				testLine << buf;
				for(int i=0;i<numValuesPerLine;++i){
					name.read_to_delim(testLine);
					name.trim_blanks();
					if(name.empty()) throw TException("Unequal number of columns in header line and Line " + (my_string) (int) lineIterator + "!\n", _FATAL_ERROR);
					simData[i]=name.toDouble();
				}
				name.read_to_delim(testLine);
				name.trim_blanks();
				if(!name.empty()) throw TException("Unequal number of columns in header line and Line " + (my_string) (int) lineIterator + "!\n", _FATAL_ERROR);

				//write parameters to outputfile
				for(int i=0; i<numParams; ++i){
					outputfile << simData[i];
					if(i<numParams-1) outputfile << "\t";
				}

				//calculate and write Linear Combinations
				myLC->calcSimDataPCA(simData);
				myLC->writeSimPCA(outputfile);
				outputfile << endl;
				++copiedLines;
			}
		}
	 cout << "Successful Program termination, " <<  copiedLines << " copied.\n";
	} catch (TException & error){
		 cout <<"\nERROR: " << error.getMessage() << endl;
		 showExplenations();
		 return 0;
	}

}
//---------------------------------------------------------------------------
