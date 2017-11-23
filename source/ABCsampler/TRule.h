//---------------------------------------------------------------------------

#ifndef TRuleH
#define TRuleH
#include <vector>
#include "TPrior.h"

//---------------------------------------------------------------------------


class TRule{
   public:
	  TSimplePrior *first, *second;
	  char sign;
	  bool isEquation;
	  my_string equation;
	  TRule(TSimplePrior* First, char Sign, TSimplePrior* Second);
	  TRule(){};
	  ~TRule(){}
};
class TEquationRule:public TRule{
   public:
	   TEquationRule(TSimplePrior* First, char Sign, my_string Equation);
	  ~TEquationRule(){}
};
//---------------------------------------------------------------------------
class TRuleVector{
   public:
	  vector<TRule> vecRule;
	  vector<TRule>::iterator curRule, endRule;

	  TRuleVector(){}
      ~TRuleVector(){ vecRule.clear();  }
};
//---------------------------------------------------------------------------
#endif
