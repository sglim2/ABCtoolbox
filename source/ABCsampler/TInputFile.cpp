//---------------------------------------------------------------------------

#pragma hdrstop

#include "TInputFile.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
//---------------------------------------------------------------------------
TInputFile::TInputFile(my_string siminputname, TPriorVector* gotPriors){
	//siminputname may contain several files....
	my_string buf, extension;
	while(!siminputname.empty()){
		buf=siminputname.extract_sub_str('#');
		siminputname.remove(0,1);
		buf.trim_blanks();
		if(!buf.empty()){
			simulationProgramInputFilenames.push_back(buf);
			extension=buf.get_file_extension();
			buf.remove_file_extension();
			buf+="-temp";
			buf+=extension;
			newSimulationProgramInputFilenames.push_back(buf);
		}
	}
   priors=gotPriors;
   checkParametersPassedToSimulationProgram=false;
   checkParametersPassedToScriptAfterSimulations=false;
}
//---------------------------------------------------------------------------
void TInputFile::writeNamesToLogfile(ofstream* logfile){
	for(int i=0;i<simulationProgramInputFilenames.size();++i){
		*logfile << " -'" <<  simulationProgramInputFilenames[i] << "'\n";
	}
}
//---------------------------------------------------------------------------
void TInputFile::writeNewNamesToLogfile(ofstream* logfile){
	for(int i=0;i<simulationProgramInputFilenames.size();++i){
		*logfile << " -'" <<  newSimulationProgramInputFilenames[i] << "'\n";
	}
}
//---------------------------------------------------------------------------
/** creats a new input file *.par */
void TInputFile::createNewInputFile(){
	//if several files are in the list parse all of them
	char delim;
	my_string buf;
	double value;

	for(int i=0;i<simulationProgramInputFilenames.size();++i){
		//open the simulation program input file SPIF
		ifstream spif(simulationProgramInputFilenames[i].c_str());
		if (!spif) throw TException("The simulation program input file '" + simulationProgramInputFilenames[i] + "' could not be read!", _FATAL_ERROR);
		//create a new file
		ofstream newSpif(newSimulationProgramInputFilenames[i].c_str());
		while(spif.good() && !spif.eof()){
			delim=buf.read_to_delim_plus(spif);
			if(priors->writeCurValueToFileFromName(buf,newSpif)) newSpif << delim;
			else if(buf!="") newSpif << buf << delim;         // if anything
		}
		// close the files
		spif.close();
		newSpif.close();
	}
} // createNewInputFile
//------------------------------------------------------------------------------
void TInputFile::initializeSimulationProgramm(my_string gotSimulationprogram, my_string gotParameter){
	vector<my_string> simParam;
	while(!gotParameter.empty()){
	    simParam.push_back(gotParameter.extract_sub_str(*"#"));
	    gotParameter.remove(0,1);
	}
	valuesToPassToSimulationProgramm=new my_string[simParam.size()];
	numValuesToPassToSimulationProgramm=simParam.size();
	for(int i=0;i<(int)simParam.size();++i){
		//check tags
		if(simParam[i].contains("SIMINPUTNAME")){
			my_string a=simParam[i].extract_after(11);
			if(a.empty()) valuesToPassToSimulationProgramm[i]=newSimulationProgramInputFilenames[0];
			else valuesToPassToSimulationProgramm[i]=newSimulationProgramInputFilenames[a.toInt()];
		} else {
			if(simParam[i]=="SIMNUM"){
				checkParametersPassedToSimulationProgram=true;
				break;
			}
			else {
				//check priors
				if(priors->isPrior(simParam[i])){
					checkParametersPassedToSimulationProgram=true;
					break;
				}
			}
			valuesToPassToSimulationProgramm[i]=simParam[i];
		}
	}
	if(gotSimulationprogram=="INTERNALGLM"){
		if(!checkParametersPassedToSimulationProgram) throw TException("The INTERNALGLM program requires model parameter tag as arguments!", _FATAL_ERROR);
		simulationProgram=new TGLM(valuesToPassToSimulationProgramm, numValuesToPassToSimulationProgramm);
	} else {
		if(checkParametersPassedToSimulationProgram) simulationProgram=new TExecuteProgram(gotSimulationprogram.c_str());
		else simulationProgram=new TExecuteProgram(gotSimulationprogram.c_str(), valuesToPassToSimulationProgramm, numValuesToPassToSimulationProgramm);
	}
}
//------------------------------------------------------------------------------
void TInputFile::initializeScriptBeforeSimulation(my_string gotScript, my_string gotParameter){
	vector<my_string> scriptParam;
	while(!gotParameter.empty()){
		scriptParam.push_back(gotParameter.extract_sub_str(*"#"));
		gotParameter.remove(0,1);
	}
	valuesToPassToScriptBeforeSimulations=new my_string[scriptParam.size()];
	numValuesToPassToScriptBeforeSimulations=scriptParam.size();
	for(int i=0;i<scriptParam.size();++i){
		//check tags
		if(scriptParam[i].contains("SIMINPUTNAME")){
			my_string a=scriptParam[i].extract_after(11);
			if(a.empty()) valuesToPassToScriptBeforeSimulations[i]=newSimulationProgramInputFilenames[0];
			else valuesToPassToScriptBeforeSimulations[i]=newSimulationProgramInputFilenames[a.toInt()];
		} else {
			valuesToPassToScriptBeforeSimulations[i]=scriptParam[i];
			if(scriptParam[i]=="SIMNUM") checkParametersPassedToScriptBeforeSimulations=true;
			else {
				//check priors
				if(priors->isPrior(scriptParam[i])) checkParametersPassedToScriptBeforeSimulations=true;
			}
			valuesToPassToScriptBeforeSimulations[i]=scriptParam[i];
		}
	}
	if(checkParametersPassedToScriptBeforeSimulations) scriptBeforeSimulations=new TExecuteProgram(gotScript.c_str());
	else scriptBeforeSimulations=new TExecuteProgram(gotScript.c_str(), valuesToPassToScriptBeforeSimulations, numValuesToPassToScriptBeforeSimulations);

}
//------------------------------------------------------------------------------
void TInputFile::initializeScriptAfterSimulation(my_string gotScript, my_string gotParameter){
	vector<my_string> scriptParam;
	while(!gotParameter.empty()){
		scriptParam.push_back(gotParameter.extract_sub_str(*"#"));
		gotParameter.remove(0,1);
	}
	valuesToPassToScriptAfterSimulations=new my_string[scriptParam.size()];
	numValuesToPassToScriptAfterSimulations=scriptParam.size();
	for(int i=0;i<scriptParam.size();++i){
		if(scriptParam[i].contains("SIMINPUTNAME")){
			my_string a=scriptParam[i].extract_after(11);
			if(a.empty()) valuesToPassToScriptAfterSimulations[i]=newSimulationProgramInputFilenames[0];
			else valuesToPassToScriptAfterSimulations[i]=newSimulationProgramInputFilenames[a.toInt()];
		} else {
			if(scriptParam[i]=="SIMNUM") checkParametersPassedToScriptAfterSimulations=true;
			else {
				//check priors
				if(priors->isPrior(scriptParam[i])) checkParametersPassedToScriptAfterSimulations=true;
			}
			valuesToPassToScriptAfterSimulations[i]=scriptParam[i];
		}
	}
	if(checkParametersPassedToScriptAfterSimulations) scriptAfterSimulations=new TExecuteProgram(gotScript.c_str());
	else scriptAfterSimulations=new TExecuteProgram(gotScript.c_str(), valuesToPassToScriptAfterSimulations, numValuesToPassToScriptAfterSimulations);

}
//------------------------------------------------------------------------------
int TInputFile::performSimulation(int simnum){
	//adjust parameters if priors are among them
	if(checkParametersPassedToSimulationProgram){
	   my_string* values=new my_string[numValuesToPassToSimulationProgramm];
	   for(int i=0;i<numValuesToPassToSimulationProgramm;++i){
		   if(valuesToPassToSimulationProgramm[i]=="SIMNUM") values[i]=simnum;
		   else {
			   float val=priors->getValue(valuesToPassToSimulationProgramm[i]);
			   if(val!=_nan) values[i]=val;
			   else values[i]=valuesToPassToSimulationProgramm[i];
		   }
	   }
	   //run simulation
	   	return simulationProgram->execute(values, numValuesToPassToSimulationProgramm);
	} else return simulationProgram->execute();

}
//------------------------------------------------------------------------------
int TInputFile::launchScriptBeforeSimulations(int simnum){
	//adjust parameters if priors are among them
	if(checkParametersPassedToScriptBeforeSimulations){
	   my_string* values=new my_string[numValuesToPassToScriptBeforeSimulations];
	   for(int i=0;i<numValuesToPassToScriptBeforeSimulations;++i){
		   if(valuesToPassToScriptBeforeSimulations[i]=="SIMNUM") values[i]=simnum;
		   else {
			   float val=priors->getValue(valuesToPassToScriptBeforeSimulations[i]);
			   if(val!=_nan) values[i]=val;
			   else values[i]=valuesToPassToScriptBeforeSimulations[i];
		   }
	   }
	   //run simulation
	   	return scriptBeforeSimulations->execute(values, numValuesToPassToScriptBeforeSimulations);
	} else return scriptBeforeSimulations->execute();

}
//------------------------------------------------------------------------------
int TInputFile::launchScriptAfterSimulations(int simnum){
	//adjust parameters if priors are among them
	if(checkParametersPassedToScriptAfterSimulations){
	   my_string* values=new my_string[numValuesToPassToScriptAfterSimulations];
	   for(int i=0;i<numValuesToPassToScriptAfterSimulations;++i){
		   if(valuesToPassToScriptAfterSimulations[i]=="SIMNUM") values[i]=simnum;
		   else {
			   float val=priors->getValue(valuesToPassToScriptAfterSimulations[i]);
			   if(val!=_nan) values[i]=val;
			   else values[i]=valuesToPassToScriptAfterSimulations[i];
		   }
	   }
	   //run simulation
	   	return scriptAfterSimulations->execute(values, numValuesToPassToScriptAfterSimulations);
	} else return scriptAfterSimulations->execute();

}
//------------------------------------------------------------------------------
