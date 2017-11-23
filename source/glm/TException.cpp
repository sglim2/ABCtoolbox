//---------------------------------------------------------------------------

#pragma hdrstop

#include <iostream>
using namespace std;

#include "TException.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
/** Output */
ostream&
operator << (ostream& os, TException& E) {
   switch(E.getIntensity()){
      default:
      case _WARNING      : os << "\nWARNING: ";      break;
      case _ERROR        : os << "\nERROR: ";        break;
      case _FATAL_ERROR  : os << "\nFATAL ERROR: ";  break;
   }
   os << E.getMessage();
   return os;
}


