/*
 * TPopData.cpp
 *
 *  Created on: Feb 5, 2009
 *      Author: wegmannd
 */

#include "TPopData.h";

TLocus::TLocus(int gotSize){
	size=gotSize;
	alleles=new int[size];
	readAlleles=0;
	meanCalculated=false;
	varCalculated=false;
}

void TLocus::add(int allele){
	if(readAlleles==size) throw TException("More samples in population than expected!", _FATAL_ERROR);
	alleles[readAlleles]=allele;
	++readAlleles;
}

double TLocus::getMean(){
	if(!meanCalculated){
		mean=0;
		numNotMissing=0;
		for(int i=0; i<size; ++i){
			if(alleles[i]>0){
				mean+=alleles[i];
				++numNotMissing;
			}
		}
		mean=mean/numNotMissing;
		meanCalculated=true;
	}
	return mean;
}

double TLocus::getVariance(){
	if(!varCalculated){
		if(!meanCalculated) getMean();
		var=0;
		for(int i=0; i<size; ++i){
			if(alleles[i]>0) var+=(alleles[i]-mean)*(alleles[i]-mean);
		}
		var=var/(numNotMissing-1);
		varCalculated=true;
	}
	return var;
}

//------------------------------------------------------------------------
TPopData::TPopData(my_string gotName, int gotSize, bool gotIsGenotypic){
	if(gotName.empty()) throw TException("Population initialized with empty name!", _FATAL_ERROR);
	if(gotSize<1) throw TException("Population initialized with sampel size <1!", _FATAL_ERROR);
	name=gotName;
	sampleSize=gotSize;
	isGenotypic=gotIsGenotypic;
	nextPhase=0;
	readLines=0;
	numLoci=0;
}

void TPopData::readArpfileLine(my_string line){
	strstream isLine;  // create the new stream
	isLine << line;     //allocate the line to a new stream for easy reading
	my_string temp;
	int readLoci=0;
	if(nextPhase==0){
		temp.read_to_delim(isLine);
		temp.read_to_delim(isLine);
	}
	while(isLine){
		temp.read_to_delim(isLine);
		temp.trim_blanks();
		if(!temp.empty()){
			if(readLines==0){
				loci.push_back(new TLocus(sampleSize));
				++numLoci;
			}
			if(readLoci==numLoci) throw TException("Too many loci on one line!", _FATAL_ERROR);
			loci[readLoci]->add(temp.toInt());
			++readLoci;
		}
	}
	++readLines;
	if(isGenotypic) nextPhase=1-nextPhase;
}

double TPopData::getMean(){
	double mean=0;
	for(int i=0;i<numLoci;++i) mean+=loci[i]->getMean();
	return mean/numLoci;
}

double TPopData::getVariance(){
	double var=0;
	for(int i=0;i<numLoci;++i) var+=loci[i]->getVariance();
	return var/numLoci;
}

int TPopData::getNumLoci(){
	return numLoci;
}

double TPopData::getMeanOneLocus(int locus){
	if(locus>=numLoci) throw TException("This locus does not exist!", _FATAL_ERROR);
	return loci[locus]->getMean();
}

double TPopData::getVarianceOneLocus(int locus){
	if(locus>=numLoci) throw TException("This locus does not exist!", _FATAL_ERROR);
	return loci[locus]->getVariance();
}
//------------------------------------------------------------------------
TPopDataVector::TPopDataVector(my_string filename){
	ifstream is (filename.c_str()); // opening the file for reading
	if(!is) throw TException("Data file '" + filename + "' could not be opened!", _FATAL_ERROR);
	my_string buf;
	strstream isLine;  // create the new stream
	buf.read_line(is); // read the line
	isLine << buf;     //allocate the line to a new stream for easy reading

	//read file
	bool withinSample=false;
	my_string name;
	int size;
	while(is.good() && !is.eof()){
	  isLine.clear();
	  buf.read_line(is); // read the line
	  buf=buf.extract_before('#');
	  buf=buf.extract_before_doubleSlash();
	  if(!buf.empty()){
		  if(!withinSample){
			  if(buf.contains("MissingData")){
				  missingDataString=buf.extract_after('=');
				  missingDataString.trim_quotes_and_blanks();
			  }
			  if(buf.contains("GenotypicData")){
				  buf=buf.extract_after('=');
				  buf.trim_quotes_and_blanks();
				  isGenotypic=buf.toInt();
			  }
			  if(buf.contains("SampleSize")){
				  buf=buf.extract_after('=');
				  buf.trim_quotes_and_blanks();
				  size=buf.toInt();
			  }
			  if(buf.contains("SampleName")){
				  name=buf.extract_after('=');
				  name.trim_quotes_and_blanks();
			  }
			  if(buf.contains("SampleData")){
				  myPopData.push_back(new TPopData(name, size, isGenotypic));
				  withinSample=true;
			  }
		  } else{
			  if(buf.contains('}')){
				  withinSample=false;
				  name="";
				  size=0;
			  } else  myPopData.back()->readArpfileLine(buf);
		  }
	  }
   }
   is.close();
}

void TPopDataVector::writeStats(ofstream* output){
	*output << myPopData[0]->getVariance();
	for(int i=1; i<myPopData.size();++i) *output << "\t" << myPopData[i]->getVariance();
		for(int i=1; i<myPopData.size();++i){
			for(int j=0; j<(myPopData.size()-1);++j){
				if(i!=j) *output << "\t" << fabs(myPopData[i]->getMean()-myPopData[j]->getMean());
			}
		}
		*output << endl;
}

void TPopDataVector::writeStatsForEachLocusIndividually(ofstream* output){
	for(int l=0; l<myPopData[0]->getNumLoci(); ++l){
		*output << myPopData[0]->getVarianceOneLocus(l);
		for(int i=1; i<myPopData.size();++i) *output << "\t" << myPopData[i]->getVarianceOneLocus(l);
		for(int i=1; i<myPopData.size();++i){
			for(int j=0; j<(myPopData.size()-1);++j){
				if(i!=j) *output << "\t" << fabs(myPopData[i]->getMeanOneLocus(l)-myPopData[j]->getMeanOneLocus(l));
			}
		}
		*output << endl;
	}
}

void TPopDataVector::writeHeader(ofstream* output){
	*output << "var_1";
	for(int i=1; i<myPopData.size();++i) *output << "\t" << "var_" << i+1;
	for(int i=1; i<myPopData.size();++i){
		for(int j=0; j<(myPopData.size()-1);++j){
			if(i!=j) *output << "\t" << "meanDif_" << i+1 << "_" << j+1;
		}
	}
	*output << endl;
}


