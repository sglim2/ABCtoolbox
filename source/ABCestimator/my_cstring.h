#ifndef _MY_Cstring_HPP_
#define _MY_Cstring_HPP_

#include <stdlib.h>
#include <string.h>   //for functions strlen strcmp strchr
#include <iostream>  //for << operator
#include <fstream>   //for read_token
#include <stdio.h>    //for sprint
using namespace std;

const int NPOS=-1;

char *my_strcpy(char *dest, const char *src);
char *my_strcat(char *dest, const char *src);
void my_strtok_right(char *dest, const char src);   //token from right to left !!!


//Marker(Insert_End)
class my_string{
	private:
   	char* _data;
      int _del_inc;
      int _nb_char;
      int _tot_space;
      char* _reserve;
   protected:
   	int tot_space(){return _tot_space;};
   public:

   	my_string(){
      	_data=new char[10];
         _data[0]='\0';
         _del_inc=10;
         _tot_space=10;
      	_nb_char=0;
         _reserve=NULL;
      };

   	my_string(const char& c){
      	_data=new char[10];
         _data[0]=c;
         _data[1]='\0';
         _del_inc=10;
         _tot_space=10;
      	_nb_char=1;
         _reserve=NULL;
      };

      my_string(int i){              // changed by Samuel 12.08.2004 due to cast (my_string) "int"
      	_data=new char[22];
         _del_inc=10;
         _tot_space=22;
         _reserve=NULL;

         sprintf(_data,"%d",i);
         _nb_char=strlen(_data);

      }

      my_string(float f){              // changed by Samuel 12.08.2004 due to cast (my_string) "int"
      	_data=new char[22];
         _del_inc=10;
         _tot_space=22;
         _reserve=NULL;

         sprintf(_data,"%f",f);
         _nb_char=strlen(_data);

      }

      my_string(int i, int a){       // changed by Samuel 12.08.2004 (see above) second param is fictive
      	_data=new char[i+2];
         _data[0]='\0';
         _del_inc=10;
         _tot_space=i+2;
      	_nb_char=0;
         _reserve=NULL;
      };

/*      my_string(const my_string& x) {
         _nb_char=x._nb_char;
         _del_inc=x._del_inc;
         _tot_space=x._tot_space;
         _reserve=NULL;

      	_data=new char[_tot_space];
         for(int i=0; i<_tot_space; ++i) _data[i]=x._data[i];
      };
*/
      my_string(const my_string& x) {
      	_data=new char[10];
         _data[0]='\0';
         _del_inc=10;
         _tot_space=10;
      	_nb_char=0;
      	*this=x;
         _reserve=NULL;
      };


      my_string(char* x) {
         _nb_char=strlen(x);
         _del_inc=10;
         if(_nb_char<10000){
            _data=new char[_nb_char+2];
            my_strcpy(_data,x);
            _tot_space=_nb_char+2;
         }
         else{
      		_data=new char[10];
         	_data[0]='\0';
         	_del_inc=10;
         	_tot_space=10;
      		_nb_char=0;
         }
      	_reserve=NULL;
      };

      ~my_string(){if(_data) delete[] _data; if(_reserve) delete[] _reserve;};

      char* c_str() const {return _data;};
      int length() const {return _nb_char;};
      void Lower();
      void Upper();
      void read_token(ifstream& is);
      //void read_token(strstream& is);
      my_string& operator=(const my_string& x);
      my_string& operator=(char* x);
      my_string& operator=(const char& x);
      my_string& operator=(int x);
      my_string& operator=(float x);
      my_string& operator+=(const my_string& x);
      my_string& operator+=(char* x);
      my_string& operator+=(char x); // new
      my_string& operator+=(const int& x); // new Samuel 12.08.2004
      char operator[](int pos) const {return _data[pos];};  //new
      char&  operator[](int pos) {return _data[pos];};  //new
      friend my_string operator+(const my_string& x, const my_string& y);
      friend my_string operator+(const my_string& x ,char* y);
      friend my_string operator+(char* x ,const my_string& y);
      friend my_string operator+(const my_string& x ,const int& y);
      friend my_string operator+(const int& x ,const my_string& y);
      friend my_string operator+(const my_string& x ,const float& y);
      friend ostream& operator<<(ostream& os, const my_string& x);     // new Samuel 12.08.2004
      friend istream& operator>>(istream& is, my_string& s);
      int operator==(const my_string& x) const{return !(strcmp(_data, x.c_str()));};
      int operator==(char* x) const {return !(strcmp(_data, x));};
      int operator!=(const my_string& x) const {return (strcmp(_data, x.c_str()));};
      int operator!=(char* x)  const {return (strcmp(_data, x));};
      int operator<(const my_string& x)const {return ((strcmp(_data, x.c_str()))< 0) ;};  //new
      int operator>(const my_string& x)const {return ((strcmp(_data, x.c_str()))> 0) ;};  //new
      //added since 21_3_97:
      void to_lower();      //converts my_string to lower case
      void to_upper();      //converts my_string to upper case
      void read_to_delim(istream& is, char delim); //read to delimiter
      void read_to_delim(istream& is); //read to delimiter
      char read_to_delim_plus(istream& is); //read to delimiter and return the delim
      my_string extract_sub_str(char delim); //extracts until the next occurence of delim or end of string //new 11_11_97
      my_string& remove( int pos, int n );
      void remove_file_extension();
      int find_first_of(const char& s, int pos ) const;
      int find_first_of(const char& s) const;
      int find_last_of(const char& s) const;
      int find_last_of(const char& s, const int& pos) const;
      int read_line(istream& is);
      bool empty();
      bool contains(const char* pat) const;
      bool contains(char pat) const;
      bool contains(const my_string& s) const;
      void assign( const my_string& s ) {*this=s;};
      int find( const my_string& s ); //return the position of patter s if not found return -1
      void rm_path(); //stef_28_3_98 removes path info from a file name
      void extract_path();   //stef_4_12_98  //extract path before a complete file name
      my_string extract_filename(); // Samuel 26-07-04 // extract the filename of a path
      my_string extract_after(char delim); //returns the part of the string after the char delim
                                            //if delim does not exist, it returns the whole string
      my_string extract_after(int index);
      my_string extract_before(char delim); //returns the part of the string before the char delim
                                            //if delim does not exist, it returns the whole string
      my_string extract_before_doubleSlash(); // returns the part of the string before the double slash "//"
      my_string extract_before(int index);
      my_string extract_between(int first, int last);
      void remove_blanks(); //rm all ' ' and '\t'    //stef_20_11_98
      int  GetNextQuota(ifstream&  ifs); //in a stream if you have "q1  l"  "q2  k" ... reads
      								// [q1  l], used in the reading of Population labels in "a_mantel.cpp"
      //Loro_14_12_98
      void convert_to_slash();
      void convert_to_html_slash(); //stef_20_4_99  converts ':' into / suitable for netscape
      void convert_to_underline();
      my_string convert_blank_to_html_blank();

      int toInt();         // Samuel 30_03_2005
      double toDouble();   // Samuel 30_03_2005

};

#endif








