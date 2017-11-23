//---------------------------------------------------------------------------

#ifndef TExceptionH
#define TExceptionH

#include "my_cstring.h"



//---------------------------------------------------------------------------
enum EXCEPTION {_WARNING, _ERROR, _FATAL_ERROR};

class TException{
//   friend istream& operator >> (istream& is, TException& E);
   friend ostream& operator << (ostream& os, TException& E);
   private:
      my_string message;
      EXCEPTION intensity; //0: warning; 1: error; 2: fatal error

   public:
      TException(){
         message="No error message!";
         intensity=_WARNING;
      }

      TException(my_string error){
         message = error;
         intensity=_WARNING;
      }

      TException(my_string error, int i){
         message = error;
         intensity=(EXCEPTION) i;
      }

      TException(my_string error, EXCEPTION i){
         message = error;
         intensity= i;
      }

      ~TException(){}

      EXCEPTION getIntensity(){return intensity;}
      my_string getMessage(){return message;}
};

#endif
