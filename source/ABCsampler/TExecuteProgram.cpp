//---------------------------------------------------------------------------

#pragma hdrstop

#include "TExecuteProgram.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
//---------------------------------------------------------------------------
TExecuteProgram::TExecuteProgram(my_string exe, my_string* parameter, int numParameter, int passedParamIdentifier){
   myExe=exe;
   myParameter=new char*[numParameter+2];
   myParameter[0]=myExe.c_str();
   for(int i=0; i<numParameter; ++i){
	   myParameter[i+1]=parameter[i].c_str();
   }
   myParameter[numParameter+1]=NULL;
   myPassedParamIdentifier=passedParamIdentifier+1;
}
TExecuteProgram::TExecuteProgram(my_string exe){
   myExe=exe;
   myParameter=new char*[2];
   myParameter[0]=myExe.c_str();
   myParameter[1]=NULL;
   myPassedParamIdentifier=-1;
}

//---------------------------------------------------------------------------
int TExecuteProgram::execute(){
	double returnValue=0.0;
		try{
			#ifdef _GCC_
			// here should be a test if the file exists.....
			int pid, status;
			pid = fork();
			if (pid == -1)
			   throw TException("Calling program '" + myExe + "': not able to create fork!", _FATAL_ERROR);
			if (pid == 0)	execv(myParameter[0], myParameter);

			// wait that the process terminates
			if (wait (&status) != pid) {
				throw TException("Error when executing '" + myExe +"'!", _FATAL_ERROR);
			}
			if (WIFEXITED (status)) {
				  returnValue=WEXITSTATUS(status);
			} else {
	           returnValue=0;
			}
			#else
			returnValue=spawnv(P_WAIT, myParameter[0], myParameter);
			if(returnValue==-1)
				throw TException("Error when executing '" + myExe +"'!" + myExe + "'!", _FATAL_ERROR);
			#endif
		  }
		  catch(TException error){
		     throw;
		  }
		  catch (...){
		     throw TException("Error when executing '" + myExe +"'!", _FATAL_ERROR);
		  }
		return returnValue;

} // executeProgram
//---------------------------------------------------------------------------
int TExecuteProgram::execute(my_string passedParam){
	// parameters to pass have to be char*!!!!!!
	if(myPassedParamIdentifier>0 && passedParam!="") myParameter[myPassedParamIdentifier]=passedParam.c_str();
	return execute();
}
//---------------------------------------------------------------------------
int TExecuteProgram::execute(my_string* passedParams, int numPassedParams){
	myParameter=new char*[numPassedParams+2];
	myParameter[0]=myExe.c_str();
	for(int i=0; i<=numPassedParams; ++i) myParameter[i+1]=passedParams[i].c_str();
	myParameter[numPassedParams+1]=NULL;
	return execute();
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
TGLM::TGLM(my_string* passedParams, int numPassedParams){
	if(numPassedParams<3) throw TException("The INTERNALGLM program requires at least three parameters!", _FATAL_ERROR);

	//the first argument is the name of the matrix file
	my_string matricesFilename=passedParams[0];
	//the second argument is the output filename
	outputFilename=passedParams[1];

	//read file with matrix definitions
	ifstream is;
	is.open(matricesFilename.c_str());
	if(!is) throw TException("The file with the definitions of the matrices '"+matricesFilename+"' could not be opend!", _FATAL_ERROR);
	my_string temp, buf;
	buf.read_line(is);
	buf.trim_blanks();
	if(buf!="C") throw TException("The file with the definitions of the matrices '"+matricesFilename+"' does not start with matrix C (tag missing?)!", _FATAL_ERROR);

	//read matrix C as a vector of arrays
	vector<double> firstLine;
	vector<double*> Ctemp;
	//read first line
	buf.read_line(is);
	buf.trim_blanks();
	while(!buf.empty()){
		temp=buf.extract_sub_str_before_ws();
		temp.trim_blanks();
		firstLine.push_back(temp.toDouble());
	}
	numParams=firstLine.size();
	Ctemp.push_back(new double[numParams]);
	for(int i=0;i<numParams;++i) Ctemp[0][i]=firstLine[i];
	//read all the other lines
	buf.read_line(is);
	while(!buf.contains("c0")){
		buf.trim_blanks();
		if(!buf.empty()){
			Ctemp.push_back(new double[numParams]);
			int i=0;
			while(!buf.empty()){
				if(i==numParams) throw TException("The file with the definitions of the matrices '"+matricesFilename+"' contains unequal number of values on the different lines specifying the matrix C!", _FATAL_ERROR);
				temp=buf.extract_sub_str_before_ws();
				temp.trim_blanks();
				Ctemp[Ctemp.size()-1][i]=temp.toDouble();
				++i;
			}
			if(i<numParams) throw TException("The file with the definitions of the matrices '"+matricesFilename+"' contains unequal number of values on the different lines specifying the matrix C!", _FATAL_ERROR);
		}
		buf.read_line(is);
	}
	numStats=Ctemp.size();

	C.ReSize(numStats, numParams);
	for(int i=0; i<numStats; ++i){
		for(int j=0; j<numParams; ++j){
			C.element(i,j)=Ctemp[i][j];
		}
	}

	//read the c0 vector
	c0.ReSize(numStats);
	buf.read_line(is);
	buf.trim_blanks();
	int i=0;
	while(!buf.empty()){
		if(i==numStats) throw TException("The file with the definitions of the matrices '"+matricesFilename+"' contains too many values for c0!", _FATAL_ERROR);
		temp=buf.extract_sub_str_before_ws();
		temp.trim_blanks();
		c0.element(i)=temp.toDouble();
		++i;
	}
	if(i<numStats) throw TException("The file with the definitions of the matrices '"+matricesFilename+"' contains too few values for c0!", _FATAL_ERROR);


	//read variances
	buf.read_line(is);
	buf.trim_blanks();
	if(buf!="Sigma") throw TException("The file with the definitions of the matrices '"+matricesFilename+"' does not contain the matrix Sigma (tag missing?)!", _FATAL_ERROR);
	Sigma.ReSize(numStats);

	for(int line=0; line<numStats;++line){
		buf.read_line(is);
		buf.trim_blanks();
		if(buf.empty()) throw TException("The file with the definitions of the matrices '"+matricesFilename+"' contains too few rows for Sigma!", _FATAL_ERROR);
		i=0;
		while(!buf.empty()){
			if(i==numStats) throw TException("The file with the definitions of the matrices '"+matricesFilename+"' contains too many values for Sigma on line "+line+"!", _FATAL_ERROR);
			temp=buf.extract_sub_str_before_ws();
			temp.trim_blanks();
			Sigma.element(line, i)=temp.toDouble();
			++i;
		}
		if(i<numStats) throw TException("The file with the definitions of the matrices '"+matricesFilename+"' contains too few values for Sigma on line "+line+"!", _FATAL_ERROR);
	}
	buf.read_line(is);
	buf.trim_blanks();
	if(!buf.empty()) throw TException("The file with the definitions of the matrices '"+matricesFilename+"' contains too few rows for Sigma!!", _FATAL_ERROR);
	is.close();
	//DONE reading matrix file....

	//check the number of passed parameters
	if(numPassedParams<(2+numParams)) throw TException("Too few parameters passed to the INTERNALGLM program!", _FATAL_ERROR);
	if(numPassedParams>(2+numParams)) throw TException("Too many parameters passed to the INTERNALGLM program!", _FATAL_ERROR);

	//prepare some matrices
	try{
		A=Cholesky(Sigma);
	} catch (...){
		throw TException("INTERNALGLM program: problem solving the Cholesky decomposition of Sigma!", _FATAL_ERROR);
	}
	e.ReSize(numStats);
	P.ReSize(numParams);
	s.ReSize(numStats);
}
//---------------------------------------------------------------------------
int TGLM::execute(my_string* passedParams, int numPassedParams){
	//read params from input --> first two are filenames
	for(int i=0;i<numParams;++i){
		P.element(i)=passedParams[i+2].toDouble();
	}

	//write header to output file
	ofstream out;
	out.open(outputFilename.c_str());
	if(!out) throw TException("INTERNALGLM program: the output file '"+outputFilename+"' could not be opened!", _FATAL_ERROR);
	out << "Stat_1";
	for(int i=1;i<numStats;++i){
		out << "\t" << "Stat_" << i+1;
	}
	out << endl;

	for(int i=0;i<numStats;++i) e.element(i)=NormalRandom(0,1);
	e=A*e;

	//compute stats
	s=C*P+c0+e;

	//write stats
	for(int i=0;i<numStats;++i){
		out << s.element(i);
	}
	out << endl;
	out.close();
	return 1;
}
