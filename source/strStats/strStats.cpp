/*
 * strStats.cpp
 *
 *  Created on: Feb 5, 2009
 *      Author: wegmannd
 */

//---------------------------------------------------------------------------
#include "TException.h"
#include "TPopData.h"
//---------------------------------------------------------------------------

int main(int argc, char* argv[]){
	try{
		if (argc!=6){
			cout << "You have specified " << (argc-1) << " parameters instead of 3!"
				 << "\n* The following parameters are needed:"
				 << "\n*  1 :  name of the arlequin file (.arp)"
				 << "\n*  2 :  name of the output file"
				 << "\n*  3 :  write header (0/1)"
				 << "\n*  4 :  overwrite or add (0/1)"
				 << "\n*  5 :  for each locus individually (1) or mean over loci (0)\n";
			throw TException("Wrong number of parameters", _FATAL_ERROR);
		}
		//read arlequin file
		TPopDataVector* myData=new TPopDataVector(argv[1]);
		//prepare output file
		ofstream output;
		my_string outname=argv[2];
		outname.trim_quotes_and_blanks();
		if(outname.empty()) throw TException("No valid output filename!", _FATAL_ERROR);
		my_string temp=argv[4];
		if(temp.toInt()==1) output.open(outname.c_str(), std::ios::app);
		else output.open(outname.c_str());
		temp=argv[3];
		if(temp.toInt()==1) myData->writeHeader(&output);
		temp=argv[5];
		if(temp.toInt()==1) myData->writeStatsForEachLocusIndividually(&output);
		else myData->writeStats(&output);
		output.close();
	 }
	 catch (TException & error){
		 cout <<"\nERROR: " << error.getMessage() << endl;
	 }
	 catch (...){
		 cout <<"\nERROR: unhandeled error!" << "\n"  << endl;
	 }
	 return 0;
}

