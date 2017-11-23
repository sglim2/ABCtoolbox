//---------------------------------------------------------------------------


#include "TRule.h"
//---------------------------------------------------------------------------
TRule::TRule(TSimplePrior* First, char Sign, TSimplePrior* Second){
		 first=First;
		 sign=Sign;
		 second=Second;
		 isEquation=false;
	  }
TEquationRule::TEquationRule(TSimplePrior* First, char Sign, my_string Equation){
		 first=First;
		 sign=Sign;
		 equation=Equation;
		 isEquation=true;
	  }

