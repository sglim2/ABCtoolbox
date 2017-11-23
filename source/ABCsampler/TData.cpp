//------------------------------------------------------------------------------

#pragma hdrstop

#include "TData.h"
//------------------------------------------------------------------------------
#pragma package(smart_init)
//------------------------------------------------------------------------------
TData::TData(my_string obsname, my_string simDataName, ofstream* gotlogfile){
	simDataFileName=simDataName;
	logFile=gotlogfile;
	obsFileName=obsname;
	readObsfile();
	checkParametersPassedToSumStatProgram=false;
	checkParametersPassedToScriptAfterSSCalc=false;
	checkParametersPassedToScriptBeforeSSCalc=false;
}
TData::TData(my_string obsname, ofstream* gotlogfile){
	simDataFileName="";
	logFile=gotlogfile;
	obsFileName=obsname;
	readObsfile();
	checkParametersPassedToSumStatProgram=false;
	checkParametersPassedToScriptAfterSSCalc=false;
	checkParametersPassedToScriptBeforeSSCalc=false;
}

void TData::readObsfile(){
	*logFile << "Reading observed data file '" << obsFileName << "' ... ";
      //Vector
      vector<my_string> vecLine;

      // open file
	  ifstream is (obsFileName.c_str()); // opening the file for reading
	  if(!is) throw TException("The .obs file '" + obsFileName + "' could not be opened!", _FATAL_ERROR);

      // Reading the file
      my_string buf;
      my_string data_temp;
      double my_data;
      strstream headerLine, dataLine;     // create the new stream for the header Line

      //Read observed Header
      buf.read_line(is);        // read the line
      headerLine << buf;        // allocate the line to a new stream for easy reading
      while(headerLine){
         data_temp.read_to_delim(headerLine);
         data_temp.trim_blanks();
         if(!data_temp.empty()) vecLine.push_back(data_temp);
      }
		numObsData=(int)vecLine.size();

		*logFile << "done, read " << numObsData << " Data\n";
        obsDataName = new my_string[numObsData];
		obsData = new double[numObsData];
		oldData = new double[numObsData];
		obsDataVariance=new double[numObsData];
		simDataPointers = new double*[numObsData];

      for(int i=0; i<numObsData; ++i){
			obsDataName[i]= vecLine[i];
      }

      //read observed Data
      buf.read_line(is);   // read the line
	  dataLine << buf;        // allocate the line to a new stream for easy reading
      for(int i=0; i<numObsData; ++i){
         data_temp.read_to_delim(dataLine);
         obsData[i]=data_temp.toDouble();
      }
	  is.close();
};
//------------------------------------------------------------------------------
void TData::initializeScriptBeforeSSCalc(my_string gotScript, my_string gotParameter, my_string sumStatsTempFileName){
	vector<my_string> param;
	while(!gotParameter.empty()){
		param.push_back(gotParameter.extract_sub_str(*"#"));
	    gotParameter.remove(0,1);
	}
	valuesToPassToScriptBeforeSSCalc=new my_string[param.size()];
	numValuesToPassToScriptBeforeSSCalc=param.size();
	for(int i=0;i<param.size();++i){
		//check tags
		if(param[i]=="SIMDATANAME"){
			valuesToPassToScriptBeforeSSCalc[i]=simDataFileName;
		} else {
			if(param[i]=="SSFILENAME"){
				valuesToPassToScriptBeforeSSCalc[i]=sumStatsTempFileName;
			} else {
				if(param[i]=="SIMNUM") checkParametersPassedToScriptBeforeSSCalc=true;
				valuesToPassToScriptBeforeSSCalc[i]=param[i];
			}
		}
	}
	if(checkParametersPassedToScriptBeforeSSCalc) scriptBeforeSSCalc=new TExecuteProgram(gotScript.c_str());
	else scriptBeforeSSCalc=new TExecuteProgram(gotScript.c_str(), valuesToPassToScriptBeforeSSCalc, numValuesToPassToScriptBeforeSSCalc);
}
//------------------------------------------------------------------------------
void TData::initializeScriptAfterSSCalc(my_string gotScript, my_string gotParameter, my_string sumStatsTempFileName){
	vector<my_string> param;
	while(!gotParameter.empty()){
		param.push_back(gotParameter.extract_sub_str(*"#"));
	    gotParameter.remove(0,1);
	}
	valuesToPassToScriptAfterSSCalc=new my_string[param.size()];
	numValuesToPassToScriptAfterSSCalc=param.size();
	for(int i=0;i<param.size();++i){
		//check tags
		if(param[i]=="SIMDATANAME"){
			valuesToPassToScriptAfterSSCalc[i]=simDataFileName;
		} else {
			if(param[i]=="SSFILENAME"){
				valuesToPassToScriptAfterSSCalc[i]=sumStatsTempFileName;
			} else {
				if(param[i]=="SIMNUM") checkParametersPassedToScriptAfterSSCalc=true;
				valuesToPassToScriptAfterSSCalc[i]=param[i];
			}
		}
	}
	if(checkParametersPassedToScriptAfterSSCalc) scriptAfterSSCalc=new TExecuteProgram(gotScript.c_str());
	else scriptAfterSSCalc=new TExecuteProgram(gotScript.c_str(), valuesToPassToScriptAfterSSCalc, numValuesToPassToScriptAfterSSCalc);
}
//------------------------------------------------------------------------------
void TData::initializeSumstatProgramm(my_string gotSumstatprogram, my_string gotParameter, my_string sumStatsTempFileName){
	vector<my_string> param;
	while(!gotParameter.empty()){
		param.push_back(gotParameter.extract_sub_str(*"#"));
	    gotParameter.remove(0,1);
	}
	valuesToPassToSumStatProgram=new my_string[param.size()];
	numValuesToPassToSumStatProgram=param.size();
	for(int i=0;i<param.size();++i){
		//check tags
		if(param[i]=="SIMDATANAME"){
			valuesToPassToSumStatProgram[i]=simDataFileName;
		} else {
			if(param[i]=="SSFILENAME"){
				valuesToPassToSumStatProgram[i]=sumStatsTempFileName;
			} else {
				if(param[i]=="SIMNUM") checkParametersPassedToSumStatProgram=true;
				valuesToPassToSumStatProgram[i]=param[i];
			}
		}
	}
	if(checkParametersPassedToSumStatProgram) sumStatProgram=new TExecuteProgram(gotSumstatprogram.c_str());
	else sumStatProgram=new TExecuteProgram(gotSumstatprogram.c_str(), valuesToPassToSumStatProgram, numValuesToPassToSumStatProgram);
}
//------------------------------------------------------------------------------
int TData::launchScriptBeforeSSCalc(int simnum){
	//adjust parameters if priors are among them
	if(checkParametersPassedToScriptBeforeSSCalc){
	   my_string* values=new my_string[numValuesToPassToScriptBeforeSSCalc];
	   for(int i=0;i<numValuesToPassToScriptBeforeSSCalc;++i){
		   if(valuesToPassToScriptBeforeSSCalc[i]=="SIMNUM") values[i]=simnum;
		   else values[i]=valuesToPassToScriptBeforeSSCalc[i];
	   }
	   //run simulation
	   	return scriptBeforeSSCalc->execute(values, numValuesToPassToScriptBeforeSSCalc);
	} else return scriptBeforeSSCalc->execute();

}
//------------------------------------------------------------------------------
int TData::launchScriptAfterSSCalc(int simnum){
	//adjust parameters if priors are among them
	if(checkParametersPassedToScriptAfterSSCalc){
	   my_string* values=new my_string[numValuesToPassToScriptAfterSSCalc];
	   for(int i=0;i<numValuesToPassToScriptAfterSSCalc;++i){
		   if(valuesToPassToScriptAfterSSCalc[i]=="SIMNUM") values[i]=simnum;
		   else values[i]=valuesToPassToScriptAfterSSCalc[i];
	   }
	   //run simulation
	   	return scriptAfterSSCalc->execute(values, numValuesToPassToScriptAfterSSCalc);
	} else return scriptAfterSSCalc->execute();

}
//------------------------------------------------------------------------------
void TData::calculateSumStats(int simnum){
	if(checkParametersPassedToSumStatProgram){
		my_string* values=new my_string[numValuesToPassToSumStatProgram];
		for(int i=0;i<numValuesToPassToSumStatProgram;++i){
		   if(valuesToPassToSumStatProgram[i]=="SIMNUM"){
			   values[i]=simnum;
		   }
		   else values[i]=valuesToPassToSumStatProgram[i];
		}
		sumStatProgram->execute(values, numValuesToPassToSumStatProgram);
	} else sumStatProgram->execute();
}
//------------------------------------------------------------------------------
void TData::readInitialSimData(my_string simFilename){
	cout << "Reading simulated datafile '" << simFilename << "' ... ";
      //Vector
		vector<my_string> vecLine;

      // open file
      ifstream is (simFilename.c_str()); // opening the file for reading
	  if(!is) throw TException("The file containing the computed summary statistics '" + simFilename + "' could not be opened!", _FATAL_ERROR);

	  // Reading the file
	  my_string buf;
	  my_string data_temp;
	  double my_data;
	  strstream headerLine, dataLine;     // create the new stream for the header Line

	  //Read Header
	  buf.read_line(is);        // read the line
	  headerLine << buf;        // allocate the line to a new stream for easy reading
	  while(headerLine){
		 data_temp.read_to_delim(headerLine);
		 data_temp.trim_blanks();
		 if(!data_temp.empty()) vecLine.push_back(data_temp);
		}
	  numSimData=vecLine.size();
	  cout << "done " << numSimData << " Data\n";
	  simDataName = new my_string[numSimData];
	  simDataInput = new double[numSimData];
	  for(int i=0; i<numSimData; ++i){
		simDataName[i]= vecLine[i];
	  }

	  //read simulated Data
      buf.read_line(is);   // read the line
      dataLine << buf;        // allocate the line to a new stream for easy reading
      readSimDataLine(dataLine);
      is.close();

		//fill simDataPointers with pointers to simDataInput where the value is requiered for distance calculations (is present in obsData)
		*logFile << "  -number of observed Statistics: " << numObsData << endl;
		*logFile << "  -number of simulated Statistics: " << numSimData << endl;
		int j;
		for(int i=0; i<numObsData;++i){
			for(j=0; j<(numSimData+1); ++j){
				if(j<numSimData && obsDataName[i]==simDataName[j]){
					simDataPointers[i]=&simDataInput[j];
					break;
				}
				if(j==numSimData)
				   throw TException("Column name '" + obsDataName[i] + "' in the file containing the computed summary statistics '" + simFilename + "' missing!", _FATAL_ERROR);
			}
		}
}
//------------------------------------------------------------------------------
void TData::readSimData(my_string simFilename){
      // open file
	  ifstream is (simFilename.c_str()); // opening the file for reading
	  if(!is) throw TException("The file containing the computed summary statistics '" + simFilename + "' could not be opened!", _FATAL_ERROR);

	  // Reading the file
      my_string buf;
      my_string data_temp;
      double my_data;
	  strstream dataLine;     // create the new stream for the header Line

	  //read simulated Data
	  buf.read_line(is); //read headerLine
      buf.read_line(is);   // read the line
      dataLine << buf;        // allocate the line to a new stream for easy reading
      readSimDataLine(dataLine);
	  is.close();
	//*logFile << "done!" << endl;
}
//------------------------------------------------------------------------------
void
TData::readSimDataLine(strstream & line){
   for(int i=0; i<numSimData; ++i){
			line >> simDataInput[i];
		}
}
//------------------------------------------------------------------------------
void
TData::writeData(ofstream& ofs){
   for(int i=0; i<numObsData;++i){
	  ofs << "\t" << *simDataPointers[i];
   }
}
//------------------------------------------------------------------------------
void
TData::writeObsData(ofstream& ofs){
   for(int i=0; i<numObsData;++i) ofs << "\t" << obsData[i];
}

//------------------------------------------------------------------------------
void
TData::writeHeader(ofstream& ofs){
	for(int i=0; i<numObsData;++i) ofs << "\t" << obsDataName[i];
}
//------------------------------------------------------------------------------
void
TData::saveOldValues(){
	for(int i=0; i<numObsData;++i) oldData[i]=*simDataPointers[i];
}
//------------------------------------------------------------------------------
void
TData::resetOldValues(){
	for(int i=0; i<numObsData;++i) *simDataPointers[i]=oldData[i];
}

//------------------------------------------------------------------------------
int TData::getObsDataIndexFromName(my_string name){
	for(int i=0; i<numObsData;++i) if(name==simDataName[i]) return i;
	return false;
}

