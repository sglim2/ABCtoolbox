/*
 * cumuldens.cpp
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
//---------------------------------------------------------------------------
void showExplanations(){
	cout << "\n**********************************************************************************"
	     << "\n* This program calculates the cumulative density at a given location in from     *"
	     << "\n* samples of the distributions. The following arguments are needed:              *"
		 << "\n*           1 : name of the file with samples from the distribution.             *"
		 << "\n*               The file must contain a header line.                             *"
		 << "\n*           2 : name of the file with the locations of interest.                 *"
		 << "\n*               The file must contain a header line.                             *"
		 << "\n*          (3): print header (just type \"header\")                                *"
		 << "\n*                                                                                *"
		 << "\n* The program writes the cumulative densities of the locations of interest to    *"
		 << "\n* the standard output in the same order as in the file with the locations of     *"
		 << "\n* interest. The file with the samples may contain additional columns and any     *"
		 << "\n* order of columns.                                                              *"
		 << "\n**********************************************************************************\n\n";
}

int main(int argc, char* argv[]){
	try{
		if (argc!=3 && argc!=4) throw TException("You have specified " + (my_string) (argc-1) + " parameters instead of 3 or 4!", _FATAL_ERROR);

		my_string temp, buf;
		bool header=false;
		if(argc==4){
			temp=argv[3];
			temp.trim_blanks();
			if(temp=="header") header=true;
			else throw TException("Unknown argument '"+temp+"'!", _FATAL_ERROR);
		}

		//read file with locations of interest
		temp=argv[2];
		ifstream lf;
		lf.open(temp.c_str());
		//first line is header....
		vector<my_string> lfHeader;
		buf.read_line(lf);
		buf.trim_blanks();
		while(!buf.empty()){
			temp=buf.extract_sub_str_before_ws();
			temp.trim_blanks();
			lfHeader.push_back(temp);
		}
		double* lfValues=new double[lfHeader.size()];
		buf.read_line(lf);
		buf.trim_blanks();
		for(int i=0;i<lfHeader.size();++i){
			if(buf.empty()) throw TException("Less values than header columns in file with locations of interest!", _FATAL_ERROR);
			temp=buf.extract_sub_str_before_ws();
			temp.trim_blanks();
			lfValues[i]=temp.toDouble();
		}
		buf.trim_blanks();
		if(!buf.empty()) throw TException("More values than header columns in file with locations of interest!", _FATAL_ERROR);
		lf.close();

		//read file with samples...
		temp=argv[1];
		ifstream sf;
		sf.open(temp.c_str());
		//first line is header....
		vector<my_string> sfHeader;
		buf.read_line(sf);
		buf.trim_blanks();
		while(!buf.empty()){
			temp=buf.extract_sub_str_before_ws();
			temp.trim_blanks();
			sfHeader.push_back(temp);
		}
		//find corresponding headers
		int* locationToSamples=new int[lfHeader.size()];
		for(int i=0;i<lfHeader.size();++i){
			for(int j=0; j<sfHeader.size()+1;++j){
				if(j==sfHeader.size()) throw TException("Column '"+lfHeader[i]+"' not found in the file with samples!", _FATAL_ERROR);
				if(lfHeader[i]==sfHeader[j]){
					locationToSamples[i]=j;
					break;
				}
			}
		}

		//now go through the file with samples and store the first smaller and first above counts....
		double* justBelow=new double[lfHeader.size()];
		double* justAbove=new double[lfHeader.size()];
		int* numAbove=new int[lfHeader.size()];
		int* numBelow=new int[lfHeader.size()];
		for(int i=0;i<lfHeader.size();++i){
			justBelow[i]=-999999;
			justAbove[i]=999999;
			numAbove[i]=0;
			numBelow[i]=0;
		}

		vector<double> values;
		int line=0;
		while(sf.good() && !sf.eof()){
			buf.read_line(sf); // read the line
			buf.trim_blanks();
			if(!buf.empty()){
				++line;
				values.clear();
				while(!buf.empty()){
					temp=buf.extract_sub_str_before_ws();
					temp.trim_blanks();
					values.push_back(temp.toDouble());
				}
				if(values.size()!=sfHeader.size()) throw TException("In the file with samples, wrong number of values on line "+(my_string) line+"!", _FATAL_ERROR);
				//test values...
				for(int i=0;i<lfHeader.size();++i){
					if(values[locationToSamples[i]]>lfValues[i]){
						++numAbove[i];
						if(values[locationToSamples[i]]<justAbove[i]) justAbove[i]=values[locationToSamples[i]];
					}
					if(values[locationToSamples[i]]<lfValues[i]){
						++numBelow[i];
						if(values[locationToSamples[i]]>justBelow[i]) justBelow[i]=values[locationToSamples[i]];
					}
				}
			}
		}

		//calculate cumulative density by interpolating and output
		if(header){
			cout << lfHeader[0];
			for(int i=1;i<lfHeader.size();++i){
				cout << "\t" << lfHeader[i];
			}
			cout << endl;
		}
		for(int i=0;i<lfHeader.size();++i){
			double dens=((lfValues[i]-(double)justBelow[i])/((double)justAbove[i]-(double)justBelow[i]));
			dens=dens*(double)(line-numBelow[i]-numAbove[i]+1)/(double)line;
			dens+=(double) numBelow[i]/(double)line;
			if(i>0) cout << "\t";
			cout << dens;
		}
		cout << endl;
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
