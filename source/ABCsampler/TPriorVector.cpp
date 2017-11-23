//---------------------------------------------------------------------------



#include "TPriorVector.h"
//---------------------------------------------------------------------------

istream&
operator >> (istream& is, TPriorVector& S) {

   return is;
} // operator >>
//------------------------------------------------------------------------------
TPriorVector::TPriorVector(my_string fileName, ofstream* gotLogFile){
   logFile=gotLogFile;
   numPrior=0;
   numSimplePrior=0;
   numCombinedPrior=0;
   mapSimplePrior.empty();
   mapCombinedPrior.empty();
   readPriorsAndRules(fileName);
   _Idum=1L;
}
//---------------------------------------------------------------------------
/** read the *est file and fill the array with the priors */
void TPriorVector::readPriorsAndRules(my_string fileName){
   *logFile << "Reading priors and rules '" << fileName << "' ... ";
	// open file
   ifstream is (fileName.c_str()); // opening the file for reading
   if(!is) throw TException("The .est file '" + fileName + "' could not be opened!", _FATAL_ERROR);

   vector<my_string> simplePriorsNames, combinedPriorsNames;
   vector<my_string>::iterator curName, endName;
   // Reading the file
	  my_string buf;
	  my_string my_name, my_min, my_max, my_firstParameter, my_secondParameter, my_sign, my_equation;
	  my_string myType;
	  my_string my_first, my_second;
	  int my_isInt;

	  int curSection=0; // 0: none, 1: Priors, 2:Rules

	  // read line by line
	  while(is.good() && !is.eof()){
		 // go to the prior section
		 strstream isLine;                   // create the new stream
		 buf.read_line(is);                  // read the line
		 buf=buf.extract_before_doubleSlash();   // remove the commentated part
		 if(!buf.empty()){
			if(buf.contains("[PARAMETERS]")) { curSection=1; continue;}
			if(buf.contains("[RULES]")){
				if(curSection==0) throw TException("Section '[PARAMETERS]]' is missing in the .est file '"+fileName+"'!", _FATAL_ERROR);
				curSection=2;
				continue;
			}
			if(buf.contains("[COMPLEX PARAMETERS]")){
				if(curSection==0) throw TException("Section '[PARAMETERS]]' is missing in the .est file '"+fileName+"'!", _FATAL_ERROR);
				curSection=3;
				continue;
			}

			// allocate the line to a new stream for easy reading
			isLine << buf;
			switch (curSection){
			   case 1: //read the priors
					 isLine >> my_isInt;              // is it a integer parameter?
					 isLine >> ws;
					 my_name.read_to_delim(isLine);   // name of the parameter
                     myType.read_to_delim(isLine);   // prior distribution
					 my_min.read_to_delim(isLine);                // read lower limit
					 my_max.read_to_delim(isLine);                // read upper limit
					 my_firstParameter.read_to_delim(isLine);                // read lower limit
					 my_secondParameter.read_to_delim(isLine);                // read upper limit

					 // if the values are ok put save the set as prior
					 if(!my_name.empty() && (my_isInt==0 || my_isInt==1)){
						if(myType=="unif")          mapSimplePrior[my_name]= new TUniformPrior(my_name, my_min.toDouble(), my_max.toDouble(), (bool)my_isInt);
						else if(myType=="logunif")  mapSimplePrior[my_name]= new TLogUniformPrior(my_name, my_min.toDouble(), my_max.toDouble(), (bool)my_isInt);
						else if(myType=="norm")     mapSimplePrior[my_name]= new TNormalPrior(my_name, my_min.toDouble(), my_max.toDouble(), my_firstParameter.toDouble(), my_secondParameter.toDouble(), (bool)my_isInt);
						else if(myType=="lognorm")  mapSimplePrior[my_name]= new TLogNormalPrior(my_name, my_min.toDouble(), my_max.toDouble(), my_firstParameter.toDouble(), my_secondParameter.toDouble(), (bool)my_isInt);
						//else if(myType=="gamma")    mapSimplePrior[my_name]= new TGammaPrior(my_name, my_min.toDouble(), my_max.toDouble(), (bool)my_isInt);
						else throw TException("Unknown prior type (Line: '" + buf + "')!", _FATAL_ERROR);
						simplePriorsNames.push_back(my_name);
					 } else throw TException("Problems reading prior (Line:' " + buf + "')!", _FATAL_ERROR);
					 break;
			   case 2: // read the rules
					 my_first.read_to_delim(isLine);     // name of the first parameter
					 isLine >> ws;
					 my_sign.read_to_delim(isLine);      // sign
					 isLine >> ws;
					 my_second.read_to_delim(isLine);    // name of the second parameter

					 // if the values are ok put save the set as a rule
					 if(!my_first.empty() && !my_second.empty() && (my_sign=="<" || my_sign==">")){
						TSimplePrior* myFirstPrior=getSimplePriorFromName(my_first);
						if(!myFirstPrior) throw TException("Problems reading rule: parameter '" + my_first + "' does not exist!", _FATAL_ERROR);
						TSimplePrior* mySecondPrior=getSimplePriorFromName(my_second);
						if(!mySecondPrior){
							rules.vecRule.push_back(TEquationRule(myFirstPrior, my_sign[0], my_second));
						}
						else{
							rules.vecRule.push_back(TRule(myFirstPrior, my_sign[0], mySecondPrior));
						}
						//throw TException("Problems reading rule: prior '" + my_second + "' does not exist!\n", _FATAL_ERROR);


					 } else throw TException("Problems reading rule (Line: '" + buf + "')!", _FATAL_ERROR);
					 break;
			   case 3: // read the complex parameters
							isLine >> my_isInt;
							isLine >> ws;
							my_name.read_to_delim(isLine);      // name of the parameter
					 isLine >> ws;
					 my_sign.read_to_delim(isLine);      // sign
					 isLine >> ws;
					 my_equation.read_line(isLine);      // the entire equation

					 // if the values are ok put save the set as combined paramter
							if(!my_name.empty() && !my_equation.empty() && my_sign=="="){
								mapCombinedPrior[my_name]= new TCombinedPrior(my_name, my_equation, (bool)my_isInt, logFile);
								combinedPriorsNames.push_back(my_name);
							}
					 else throw TException("Problems reading complex parameter (Line: '" + buf + "')!", _FATAL_ERROR);
					 break;
			}
		 }
	  }
	//create iterators
	endSimplePrior=mapSimplePrior.end();
	curSimplePrior=mapSimplePrior.begin();
	endCombinedPrior=mapCombinedPrior.end();
	curCombinedPrior=mapCombinedPrior.begin();
	rules.endRule=rules.vecRule.end();
	rules.curRule=rules.vecRule.begin();
	numSimplePrior=simplePriorsNames.size();
	numCombinedPrior=combinedPriorsNames.size();

	//create Arrays of simple and combined Prior pointers for fast reference
	simplePriors = new TSimplePrior*[numSimplePrior];
	combinedPriors = new TCombinedPrior*[numCombinedPrior];
    //keep order of est File: combined priors have to be in the same order as in the est file!!!
	endName=simplePriorsNames.end();
	int i=0;
	for(curName=simplePriorsNames.begin(); curName!=endName; ++curName, ++i){
	   simplePriors[i]=mapSimplePrior[*curName];
	}
	i=0;
	endName=combinedPriorsNames.end();
	for(curName=combinedPriorsNames.begin(); curName!=endName; ++curName, ++i)
	   combinedPriors[i]=mapCombinedPrior[*curName];
	is.close();
	*logFile << "done:" << endl << " -" << mapSimplePrior.size() << " Priors" << endl << " -" << mapCombinedPrior.size() << " Combined Priors" << endl << " -" << rules.vecRule.size() << " Rules\n";
} // readPriorsAndRules
//------------------------------------------------------------------------------
void TPriorVector::writeHeader(ofstream& ofs){
   writeHeaderSimplePriors(ofs);
   writeHeaderCombinedPriors(ofs);
}
void TPriorVector::writeHeaderSimplePriors(ofstream& ofs){
	for(curSimplePrior=mapSimplePrior.begin(); curSimplePrior!=endSimplePrior; ++curSimplePrior)
	   ofs << "\t" << curSimplePrior->second->name;
}
void TPriorVector::writeHeaderCombinedPriors(ofstream& ofs){
	for(curCombinedPrior=mapCombinedPrior.begin(); curCombinedPrior!=endCombinedPrior; ++curCombinedPrior)
	   ofs << "\t" << curCombinedPrior->second->name;
}
//------------------------------------------------------------------------------
void TPriorVector::writeParameters(ofstream& ofs){
   writeParametersSimplePriors(ofs);
   writeParametersCombinedPriors(ofs);
}
void TPriorVector::writeParametersSimplePriors(ofstream& ofs){
	for(curSimplePrior=mapSimplePrior.begin(); curSimplePrior!=endSimplePrior; ++curSimplePrior){
	   ofs << "\t";
	   curSimplePrior->second->writeCurValue(ofs);
    }
}
void TPriorVector::writeParametersCombinedPriors(ofstream& ofs){
	for(curCombinedPrior=mapCombinedPrior.begin(); curCombinedPrior!=endCombinedPrior; ++curCombinedPrior){
	   ofs << "\t";
	   curCombinedPrior->second->writeCurValue(ofs);
    }
}
//------------------------------------------------------------------------------
TPrior* TPriorVector::getPriorFromName(const my_string& name){
   map<my_string,TSimplePrior*>::iterator myCurSimplePrior=mapSimplePrior.find(name);
   if(myCurSimplePrior!=mapSimplePrior.end()) return myCurSimplePrior->second;
   map<my_string,TCombinedPrior*>::iterator myCurCombinedPrior=mapCombinedPrior.find(name);
   if(myCurCombinedPrior!=mapCombinedPrior.end()) return myCurCombinedPrior->second;
   return 0;
}
TSimplePrior* TPriorVector::getSimplePriorFromName(const my_string& name){
   map<my_string,TSimplePrior*>::iterator myCurSimplePrior=mapSimplePrior.find(name);
   if(myCurSimplePrior!=mapSimplePrior.end()) return myCurSimplePrior->second;
   return 0;
}
TCombinedPrior* TPriorVector::getCombinedPriorFromName(const my_string& name){
   map<my_string,TCombinedPrior*>::iterator myCurCombinedPrior=mapCombinedPrior.find(name);
   if(myCurCombinedPrior!=mapCombinedPrior.end()) return myCurCombinedPrior->second;
   return 0;
}
int TPriorVector::getNumberOfSimplePriorFromName(const my_string& name){
   for(int i=0; i<numSimplePrior;++i) if(simplePriors[i]->name==name) return i;
   return -1;
}
//------------------------------------------------------------------------------
bool TPriorVector::writeCurValueWithHyperprior(const my_string& name, ofstream& file){
   //check for hyperprior
	char key=name[0];
    my_string param=name;
	switch (key){
			case '%':{ // gamma_deviation
				int pos2 = param.find_first_of('%',1);
				my_string arg = param.extract_between(1, pos2);
				double value=getValue(arg);   // check if it is a prior name
				if(value==_nan) value=arg.toDouble();
				curPrior=getPriorFromName(param.extract_after(pos2));
				if(!curPrior) return false;
				curPrior->writeHyperPriorGamma(value, file);
				return true;
                }
			case '&': // beta_deviation
				curPrior=getPriorFromName(param.extract_after(0));
				if(!curPrior) return false;
				curPrior->writeHyperPriorBeta(file);
				return true;
			case '$': {// normal
				int pos2 = param.find_first_of('$',1);
				my_string arg = param.extract_between(1, pos2);
				double value=getValue(arg);   // check if it is a prior name
				if(value==_nan) value=arg.toDouble();
				curPrior=getPriorFromName(param.extract_after(pos2));
				if(!curPrior) return false;
				curPrior->writeHyperPriorNormal(value, file);
				return true;}
			case '!': {// normal truncated at 0 --> always positive!
				int pos2 = param.find_first_of('!',1);
				my_string arg = param.extract_between(1, pos2);
				double value=getValue(arg);   // check if it is a prior name
				if(value==_nan) value=arg.toDouble();
				curPrior=getPriorFromName(param.extract_after(pos2));
				if(!curPrior) return false;
				curPrior->writeHyperPriorNormalPositive(value, file);
				return true;}
			case '#': {//lognormal
				int pos2 = param.find_first_of('#',1);
				my_string arg = param.extract_between(1, pos2);
				double value=getValue(arg);   // check if it is a prior name
				if(value==_nan) value=arg.toDouble();
				curPrior=getPriorFromName(param.extract_after(pos2));
				if(!curPrior) return false;
				curPrior->writeHyperPriorLognormal(value, file);
				return true;}
			case '@': { //lognorm base 10, but values given in logscale!
				int pos2 = param.find_first_of('@',1);
				my_string arg = param.extract_between(1, pos2);
				double value=getValue(arg);   // check if it is a prior name
				if(value==_nan) value=arg.toDouble();
				curPrior=getPriorFromName(param.extract_after(pos2));
				if(!curPrior) return false;
				curPrior->writeHyperPriorLognormalParametersInLog10Scale(value, file);
				return true;}
	}
	return false;
}
//------------------------------------------------------------------------------
bool TPriorVector::writeCurValueToFileFromName(const my_string& name, ofstream& file){
   if(writeCurValueWithHyperprior(name, file)) return true;
   curPrior=getPriorFromName(name);
   if(curPrior){
	  curPrior->writeCurValue(file);
	  return true;
   }
   return false;
}
//------------------------------------------------------------------------------
// changes the parameter set controlling for the rules
/*
void TPriorVector::getNewValues(){
	// change all simple prior values
	for(curSimplePrior=mapSimplePrior.begin(); curSimplePrior!=endSimplePrior; ++curSimplePrior)
	   curSimplePrior->second->changeCurrentValue();
	// check the rules if the rules do not fit change the FIRST parameter
	TPrior* firstPrior;
	TPrior* secondPrior;
	double secondValue;
	vector<TRule>::iterator curRule = rules.vecRule.begin();
	for(; curRule!=rules.vecRule.end(); ++curRule){
		if(curRule->isEquation) secondValue=calcEquation(curRule->equation);
		else secondValue=curRule->second->curValue;
	  switch ((*curRule).sign){
		 case '<': if(curRule->first->curValue >= secondValue)
					  curRule->first->changeCurrentValueWithLimits(curRule->first->lowerLimit, secondValue);
				   break;
		 case '>': if(curRule->first->curValue <= secondValue)
					  curRule->first->changeCurrentValueWithLimits(secondValue, curRule->first->upperLimit);
				   break;
	  }
   }
   // calc all combined parameters
   updateCombinedParameters();
} // getNewValues
*/
void TPriorVector::getNewValues(){
	bool checkRules=false;
	while(!checkRules){
		// change all simple prior values
		for(curSimplePrior=mapSimplePrior.begin(); curSimplePrior!=endSimplePrior; ++curSimplePrior)
			   curSimplePrior->second->changeCurrentValue();
		//check the rules if the rules do not fit change ALL parameters until the rules fit
		double secondValue;
		checkRules=true;
		vector<TRule>::iterator curRule = rules.vecRule.begin();
		for(; curRule!=rules.vecRule.end(); ++curRule){
			if(curRule->isEquation) secondValue=calcEquation(curRule->equation);
			else secondValue=curRule->second->curValue;
			switch ((*curRule).sign){
				case '<': if(curRule->first->curValue >= secondValue) checkRules=false; break;
				case '>': if(curRule->first->curValue <= secondValue) checkRules=false; break;
			}
		}
	}

   // calc all combined parameters
   updateCombinedParameters();
}
//------------------------------------------------------------------------------
//functions to update parameters via MCMC
void TPriorVector::getNewValuesMcmc(){
	//first get new values
	double* oldValues=new double[numSimplePrior];
	for(int i=0;i<numSimplePrior;++i) oldValues[i]=simplePriors[i]->curValue;
	bool rulesok=false;
	while(!rulesok){
		for(int i=0;i<numSimplePrior;++i){
			simplePriors[i]->curValue=oldValues[i]+UniformRandom(0, simplePriors[i]->mcmcStep) - simplePriors[i]->mcmcStep/2;
			//reflect, if necessary
			if(simplePriors[i]->curValue < simplePriors[i]->lowerLimit) simplePriors[i]->curValue=simplePriors[i]->lowerLimit+(simplePriors[i]->lowerLimit-simplePriors[i]->curValue);
			if(simplePriors[i]->curValue > simplePriors[i]->upperLimit) simplePriors[i]->curValue=simplePriors[i]->upperLimit+(simplePriors[i]->upperLimit-simplePriors[i]->curValue);
		}
		//check the rules if the rules do not fit change ALL parameters until the rules fit
    	double secondValue;
    	rulesok=true;
    	vector<TRule>::iterator curRule = rules.vecRule.begin();
    	for(; curRule!=rules.vecRule.end(); ++curRule){
    		if(curRule->isEquation) secondValue=calcEquation(curRule->equation);
    		else secondValue=curRule->second->curValue;
    		switch ((*curRule).sign){
    			case '<': if(curRule->first->curValue >= secondValue) rulesok=false; break;
    			case '>': if(curRule->first->curValue <= secondValue) rulesok=false; break;
    		}
    	}
	}
	// calc all combined parameters
	updateCombinedParameters();
}
void TPriorVector::getNewValuesMcmcUpdateOnePriorOnly(){
   //select prior
   int r=numSimplePrior*ran3(&_Idum)-0.5;
   getNewValuesMcmc(simplePriors[r]);
   // calc all combined parameters
   updateCombinedParameters();
}
void TPriorVector::getNewValuesMcmc(const int& priorNumber){
	  getNewValuesMcmc(simplePriors[priorNumber]);
	   // calc all combined parameters
	   updateCombinedParameters();
}
void TPriorVector::getNewValuesMcmc(TSimplePrior* thisSimplePrior){
   double newValue;
   double MCMC_max=thisSimplePrior->upperLimit;
   double MCMC_min=thisSimplePrior->lowerLimit;
   double secondValue;
   rules.curRule=rules.vecRule.begin();
   for(; rules.curRule!=rules.endRule; ++rules.curRule){
	  if(rules.curRule->first->name==thisSimplePrior->name){
		  if(rules.curRule->isEquation) secondValue=calcEquation(rules.curRule->equation);
		  else secondValue=rules.curRule->second->curValue;
		  switch(rules.curRule->sign){
			case '<': if(secondValue < MCMC_max) MCMC_max = secondValue; break;
			case '>': if(secondValue > MCMC_min) MCMC_min = secondValue; break;
		}
	 }
	 if(rules.curRule->second->name==thisSimplePrior->name){
		switch(rules.curRule->sign){
		   case '<': if(rules.curRule->first->curValue > MCMC_min)
						MCMC_min = rules.curRule->first->curValue;
					 break;
		   case '>': if(rules.curRule->first->curValue < MCMC_max)
						MCMC_max = rules.curRule->first->curValue;
					 break;
		}
	 }
  }
  //Now we have the range the new prior lay in. We use a uniform prior
  //within the range current value +/- mcmcStep, which is assinged in the calibrating step
  newValue=thisSimplePrior->curValue+UniformRandom(0, thisSimplePrior->mcmcStep) - thisSimplePrior->mcmcStep/2;
  //test with normal distribution:
  //curValue=curValue+NormalRandom(0, mcmcStep);
  //If we move out of the possible range
  //given by MCMC_min and MCMC_max, we mirror the move on this border. As this
  //is a symmetrical move, we do not have to adjust the transition Kernel exept for the
  //case of the border itself, where it is 2 times more probable to move to any point
  //than to move from this point to the border.
  if(newValue < MCMC_min) newValue=MCMC_min+(MCMC_min-newValue);
  if(newValue > MCMC_max) newValue=MCMC_max+(MCMC_max-newValue);
  thisSimplePrior->setCurValue(newValue);
}
//------------------------------------------------------------------------------
//functions to update parameters via PMC
bool TPriorVector::getNewValuesPMC(double* newParams){
	//save new params into prior objects...
	for(int i=0; i<numSimplePrior; ++i) simplePriors[i]->setCurValue(newParams[i]);

	//check the rules if the rules do not fit change ALL parameters until the rules fit
	double secondValue;
	vector<TRule>::iterator curRule = rules.vecRule.begin();
	for(; curRule!=rules.vecRule.end(); ++curRule){
		if(curRule->isEquation) secondValue=calcEquation(curRule->equation);
		else secondValue=curRule->second->curValue;
		switch ((*curRule).sign){
			case '<': if(curRule->first->curValue >= secondValue) return false;
			case '>': if(curRule->first->curValue <= secondValue) return false;
		}
	}

   // calc all combined parameters
   updateCombinedParameters();
   return true;
}

//------------------------------------------------------------------------------
void TPriorVector::resetOldValues(){
	for(curSimplePrior=mapSimplePrior.begin(); curSimplePrior!=endSimplePrior; ++curSimplePrior)
	   curSimplePrior->second->resetOldValue();
	for(curCombinedPrior=mapCombinedPrior.begin(); curCombinedPrior!=endCombinedPrior; ++curCombinedPrior)
	   curCombinedPrior->second->resetOldValue();
};
//------------------------------------------------------------------------------
void TPriorVector::saveOldValues(){
   for(curSimplePrior=mapSimplePrior.begin(); curSimplePrior!=endSimplePrior; ++curSimplePrior)
	   curSimplePrior->second->saveOldValue();
   for(curCombinedPrior=mapCombinedPrior.begin(); curCombinedPrior!=endCombinedPrior; ++curCombinedPrior)
	   curCombinedPrior->second->saveOldValue();
};
//------------------------------------------------------------------------------
void TPriorVector::updateCombinedParameters(){
	//same order as in the est file!!!
	for(int i=0; i<numCombinedPrior; ++i){
	   combinedPriors[i]->setCurValue(calcEquation(combinedPriors[i]->equation));
	}
};
//------------------------------------------------------------------------------
double TPriorVector::getPriorDensity(){
   double density=1.0;
   for(curSimplePrior=mapSimplePrior.begin(); curSimplePrior!=endSimplePrior; ++curSimplePrior)
	  density*=curSimplePrior->second->getPriorDensity();
   return density;
};
double TPriorVector::getOldPriorDensity(){
   double density=1.0;
   for(curSimplePrior=mapSimplePrior.begin(); curSimplePrior!=endSimplePrior; ++curSimplePrior)
	  density*=curSimplePrior->second->getPriorDensity();
   return density;
};
double TPriorVector::getPriorDensity(double* values){
   double density=1.0;
   for(int i=0; i<numSimplePrior; ++i){
	  density*=simplePriors[i]->getPriorDensityFromValue(values[i]);
   }
   return density;
};
//------------------------------------------------------------------------------
double TPriorVector::getValue(const my_string& name){
   curPrior=getPriorFromName(name);
   if(!curPrior) return _nan;
   return curPrior->curValue;
}
//------------------------------------------------------------------------------
bool TPriorVector::isPrior(const my_string& name){
   curPrior=getPriorFromName(name);
   if(!curPrior) return false;
   return true;
}
//------------------------------------------------------------------------------
/** Change the current value of the combined parameter, e.g. resolve the equation.
  * The normal priors have to be already taken.
  */
double TPriorVector::calcEquation(my_string equat){
   my_string curString;
   double first=_nan,
		 second=_nan,
		 temp_value;
   char  sign='0';
   bool  stop;
   char  c;
   strstream equ;      // create the new stream
   equ << equat;

   // for the whole string
   while(equ.good() && !equ.eof()){
	  // reset paramters
	  curString="";
	  temp_value=_nan;
	  // remove heading spaces
	  do{
		 equ.get(c);
      } while (c && !equ.eof() && isspace(c));
	  if(equ.eof()) return first;   // if the stream was empty

	  // get token (check for brackets!)
	  do{
		 stop=false;
		 switch (c){
            default:
			   curString += c;
               break;

            // if it is a sign put it back
			case '+':
            case '-':
			case '*':
            case '/':
			   if(!curString.empty()) equ.putback(c);
               else curString += c;
			   stop=true;   // to stop the loop
               break;

            // if it is a bracket
			case '(':
               // find the closing bracket and calculate the contents recursively
				//take care of brackets within brackets!!!
			   my_string bracket;
			   int level=0;
               while((c=equ.get()) && !equ.eof() && ( c!=')' || level>0)){
            	   if(c=='(') ++level;
            	   if(c==')') --level;
            	   bracket += c;
			   }
			   if(equ.eof())    // if there is no closing bracket
				  throw TException("Problems solving equation '" + equat + "' (Closing bracket missing)!", _FATAL_ERROR);
			   // calc recursively the bracket
			   temp_value=calcEquation(bracket);

               // check if there is a function associated to the bracket
			   if(!curString.empty()){
				  if(curString=="log")          temp_value = log(temp_value);
				  else if(curString=="log10")   temp_value = log10(temp_value);
				  else if(curString=="exp")     temp_value = exp(temp_value);
				  else if(curString=="pow10")     temp_value = pow(10, temp_value);
				  else if(curString=="abs")     temp_value = fabs(temp_value);
				  else if(curString=="ispositive"){
					  if(temp_value>0) temp_value=1; else temp_value=0;
				  }
				  else if(curString=="isnegative"){
					  if(temp_value<0) temp_value=1; else temp_value=0;
				  }
				  else throw TException("Problems solving equation '" + equat + "': function '" + curString + "' not known)!", _FATAL_ERROR);
			   }
               stop=true;   // to stop the loop

         }
	  }while(!stop && equ.get(c) && !equ.eof()  && !isspace(c));

	  // search first value
	  if(first==_nan){
		 if(temp_value != _nan) first=temp_value;
			else{
			   first=getValue(curString);   // check if it is a prior name
			   if(first==_nan) first=curString.toDouble();
		 }
	  }
	  else if(sign=='0') sign=curString[0];
      else{
		 // search second value
		 if(temp_value != _nan) second=temp_value;
         else{
			second=getValue(curString);   // check if it is a prior name
				if(second==_nan) second=curString.toDouble();
         }
      }

      // if we have all values perform the calculation
		if (second !=_nan && sign!=0 && first!=_nan){
			switch(sign){
				case '+': first+=second;  break;
				case '-': first-=second;  break;
				case '*': first*=second;  break;
				case '/': if(!second) throw TException("Problems solving equation '" + equat +  "': devision by zero!", _FATAL_ERROR);
						  first/=second;
						  break;
				default:  throw TException("Problems solving equation '" + equat + "': no known sign used!", _FATAL_ERROR);
			}

			// reset the values
			sign='0';
			second=_nan;
		}
	}
	return first;
}
//------------------------------------------------------------------------------
void TPriorVector::setSimplePriorValue(const int& priorNumber, const double& value){
   simplePriors[priorNumber]->curValue=value;
   simplePriors[priorNumber]->oldValue=value;
}
