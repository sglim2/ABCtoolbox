//---------------------------------------------------------------------------

#ifndef TExecuteProgramH
#define TExecuteProgramH

#include "global.h"
#include "my_cstring.h"
#include <fstream>
#include <iostream.h>
#include <vector>
#include "newmat.h"
#include "newmatap.h"
#include "newmatio.h"
//-----------------------
#ifdef _GCC_
   #include <unistd.h>
   #include <sys/types.h>
   #include <sys/wait.h>
   #include <errno.h>
#else
	#include <process.h>
    //#include <sys/wait.h>
#endif
//---------------------------------------------------------------------------


class TExecuteProgram{
private:
	  my_string myExe;
	  char** myParameter;
	  int myPassedParamIdentifier;

public:
	TExecuteProgram(){};
	TExecuteProgram(my_string exe, my_string* parameter, int numParameter, int passedParamIdentifier=-1);
	TExecuteProgram(my_string exe);
	virtual int execute();
	virtual int execute(my_string passedParam);
	virtual int execute(my_string* passedParams, int numPassedParams);
};

class TGLM:public TExecuteProgram{
private:
	Matrix C;
	ColumnVector c0;
	SymmetricMatrix Sigma;
	Matrix A;
	ColumnVector e;
	ColumnVector P;
	ColumnVector s;
	my_string outputFilename;
	int numStats;
	int numParams;

public:
	TGLM(my_string* passedParams, int numPassedParams);
	virtual int execute(my_string* passedParams, int numPassedParams);
};

#endif
