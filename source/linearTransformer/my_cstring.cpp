#include "my_cstring.h"


//------------------------------------------------------------------------------
//diverse functions
char *my_strcpy(char *dest, const char *src){
	char *p_src=(char*) src;
   char *p_dest=dest;
   for(;;){
   	*p_dest=*p_src;
      if(*p_src=='\0') break;
      p_dest++;
      p_src++;
   }
   return dest;
};

char *my_strcat(char *dest, const char *src){

	char *p_src=(char*) src;
   char *p_dest=dest;
   //goto end of dest
   for(;*p_dest!='\0';){
   	p_dest++;
   }
   //then append:

   for(;;){
   	*p_dest=*p_src;
      if(*p_src=='\0') break;
      p_dest++;
      p_src++;
   }
   return dest;
};

void my_strtok_right(char *str, const char del){   //stef_3_8_99
	int len=strlen(str);
   int i;
   for(i=len-1; i>=0 && str[i]!=del ; --i);
   str[i]='\0';  //stop the string
};

//my_string member functions
void
my_string::Lower(){
	int len=strlen(_data);
  	for(int i=0; i<len; i++){
   	if(_data[i]=='A') _data[i]='a';
      if(_data[i]=='B') _data[i]='b';
      if(_data[i]=='C') _data[i]='c';
      if(_data[i]=='D') _data[i]='d';
      if(_data[i]=='E') _data[i]='e';
      if(_data[i]=='F') _data[i]='f';
      if(_data[i]=='G') _data[i]='g';
      if(_data[i]=='H') _data[i]='h';
      if(_data[i]=='I') _data[i]='i';
      if(_data[i]=='J') _data[i]='j';
      if(_data[i]=='K') _data[i]='k';
      if(_data[i]=='L') _data[i]='l';
      if(_data[i]=='M') _data[i]='m';
      if(_data[i]=='N') _data[i]='n';
      if(_data[i]=='O') _data[i]='o';
      if(_data[i]=='P') _data[i]='p';
      if(_data[i]=='Q') _data[i]='q';
      if(_data[i]=='R') _data[i]='r';
      if(_data[i]=='S') _data[i]='s';
      if(_data[i]=='T') _data[i]='t';
      if(_data[i]=='U') _data[i]='u';
      if(_data[i]=='V') _data[i]='v';
      if(_data[i]=='W') _data[i]='w';
      if(_data[i]=='X') _data[i]='x';
      if(_data[i]=='Y') _data[i]='y';
      if(_data[i]=='Z') _data[i]='z';
   }
};

void
my_string::Upper(){ //convert to upper case
	int len=strlen(_data);
   for(int i=0; i<len; i++){
   	if(_data[i]=='a') _data[i]='A';
      if(_data[i]=='b') _data[i]='B';
      if(_data[i]=='c') _data[i]='C';
      if(_data[i]=='d') _data[i]='D';
      if(_data[i]=='e') _data[i]='E';
      if(_data[i]=='f') _data[i]='F';
      if(_data[i]=='g') _data[i]='G';
      if(_data[i]=='h') _data[i]='H';
      if(_data[i]=='i') _data[i]='I';
      if(_data[i]=='j') _data[i]='J';
      if(_data[i]=='k') _data[i]='K';
      if(_data[i]=='l') _data[i]='L';
      if(_data[i]=='m') _data[i]='M';
      if(_data[i]=='n') _data[i]='N';
      if(_data[i]=='o') _data[i]='O';
      if(_data[i]=='p') _data[i]='P';
      if(_data[i]=='q') _data[i]='Q';
      if(_data[i]=='r') _data[i]='R';
      if(_data[i]=='s') _data[i]='S';
      if(_data[i]=='t') _data[i]='T';
      if(_data[i]=='u') _data[i]='U';
      if(_data[i]=='v') _data[i]='V';
      if(_data[i]=='w') _data[i]='W';
      if(_data[i]=='x') _data[i]='X';
      if(_data[i]=='y') _data[i]='Y';
      if(_data[i]=='z') _data[i]='Z';
   }
};

void
my_string::to_lower(){
  Lower();
};

void
my_string::to_upper(){ //convert to upper case
	Upper();
};

void   // There are bugs in this routine: do not use or debug it first
my_string::read_token(ifstream& is){ //not tested
//1) clean:
	if(_data) delete[] _data;
   _data=new char[10];
   _data[0]='\0';
   _del_inc=10;
   _tot_space=10;
   _nb_char=0;
   _reserve=NULL;
//2) skipp blanks
	char current=' ';
   for(;current==' ' || current=='\t' || current=='\n' && !is.bad();){
   	 is.get(current);  //stop at first not blank char
   }
//3)take as long no blanks
   for(;current!=' ' && current!='\t' && current!='\n';){
   	if(_nb_char>=_tot_space-1){
      	int len=_tot_space+_del_inc;
   		if(_reserve) delete[] _reserve;
      	_reserve=new char[(strlen(_data)+2)];     //make spcace for buffer
      	if(!_reserve) return ;               //memory out
      	my_strcpy(_reserve, _data);                  //copy data to buffer
   		delete[] _data;
      	_data=new char[len+2];                    //make enought space for catenate
      	_tot_space=len+2;
      	my_strcpy(_data, _reserve);                  //copy bac
      	delete[] _reserve; _reserve=NULL;         //clean buffer;
         _tot_space=len;    								//store new len
      }
      _data[_nb_char]=current;
      _nb_char++;
      _data[_nb_char]='\0';
      is.get(current);
   }
   _data[_nb_char]='\0';
   if(current=='\n') is.putback(current);     //read token is
                                              //supposed not to 'eat' newline

};

/*void
my_string::read_token(strstream& is){ //not tested
//1) clean:
	if(_data) delete[] _data;
   _data=new char[10];
   _data[0]='\0';
   _del_inc=10;
   _tot_space=10;
   _nb_char=0;
   _reserve=NULL;
//2) skipp blanks
	char current=' ';
   for(;current==' ' || current=='\t' || current=='\n';){
   	 is.get(current);  //stop at first not blank char
   }
//3)take as long no blanks
   for(;current!=' ' && current!='\t' && current!='\n';){
   	if(_nb_char>=_tot_space-1){
      	int len=_tot_space+_del_inc;
   		if(_reserve) delete[] _reserve;
      	_reserve=new char[(strlen(_data)+2)];     //make spcace for buffer
      	if(!_reserve) return ;               //memory out
      	my_strcpy(_reserve, _data);                  //copy data to buffer
   		delete[] _data;
      	_data=new char[len+2];                    //make enought space for catenate
      	_tot_space=len+2;
      	my_strcpy(_data, _reserve);                  //copy bac
      	delete[] _reserve; _reserve=NULL;         //clean buffer;
         _tot_space=len;    								//store new len
      }
      _data[_nb_char]=current;
      _nb_char++;
      _data[_nb_char]='\0';
      is.get(current);
   }
   _data[_nb_char]='\0';



};*/


//------------------------------------------------------------------------------
//copy operators
my_string&
my_string::operator=(const char& x){

   if(_tot_space < 4){
   	delete[] _data;
      _data=new char[12];
      _tot_space=12;
   }
   _nb_char=1;
   _data[0]=x;
   _data[1]='\0';
   return *this;
}

my_string&
my_string::operator=(const my_string& x){
	int len=x._nb_char;
   if(len > 10000){
		return *this;
   }
   if(_tot_space < len+2){
   	delete[] _data;
      _data=new char[len+2];
      _tot_space=len+2;
   }
   _nb_char=len;
   my_strcpy(_data, x._data);
   return *this;
}

my_string&
my_string::operator=(char* x){
	int len=strlen(x);
   if(len > 10000){
		return *this;
   }
   if(_tot_space < len+2){
   	delete[] _data;
      _data=new char[len+2];
      _tot_space=len+2;
   }
   _nb_char=len;
   my_strcpy(_data, x);
   return *this;
}

my_string&
my_string::operator=(int x){
   if(_tot_space < 20){
   	delete[] _data;
      _data=new char[22];
      _tot_space=22;
   }
   sprintf(_data,"%d",x);
   _nb_char=strlen(_data);
   return *this;
}

my_string&
my_string::operator=(float x){
   if(_tot_space < 40){
   	delete[] _data;
      _data=new char[42];
      _tot_space=42;
   }
   sprintf(_data,"%f",x);
    _nb_char=strlen(_data);
   return *this;
}

//------------------------------------------------------------------------------
//catenation operators
my_string&
my_string::operator+=(const my_string& x){
	int len=x._nb_char+_nb_char;
   if(len > 10000){
		return *this;
   }
   if(_tot_space < len+2){                   // if the array is to small to contain all
      char* reserve=_data;                   //use a temp pointer to the array
      _tot_space=len+2;
      _data=new char[_tot_space];                    //make enought space for catenate
      my_strcpy(_data, reserve);                  //copy bac
      delete[] reserve;                      //delete the old array;
   }
   _nb_char=len;
   if (_data) my_strcat(_data, x._data);
   return *this;
}

my_string&
my_string::operator+=(char* x){
	int len=strlen(x)+_nb_char;
   if(len > 10000){
		return *this;
   }
   if(_tot_space < len+2){
      char* reserve=_data;                   //use a temp pointer to the array
      _tot_space=len+2;
      _data=new char[_tot_space];                    //make enought space for catenate
      my_strcpy(_data, reserve);                  //copy bac
      delete[] reserve;                      //delete the old array;
   }
   _nb_char=len;
   if (_data) my_strcat(_data, x);
   return *this;
}


my_string&
my_string::operator+=(char x){
	int len=_nb_char+1;
   if(len > 10000){
		return *this;
   }
   if(_tot_space < len+2){
      char* reserve=_data;                   //use a temp pointer to the array
      _tot_space=len+2;
      _data=new char[_tot_space];            //make enought space for catenate
      my_strcpy(_data, reserve);             //copy the array content
      delete[] reserve;                      //delete the old array
   }

   if (_data) _data[_nb_char]=x;               //catenate character
   _data[len]='\0';
   _nb_char=len;
   return *this;

} // new

my_string&                                      // created Samuel 12.08.2004
my_string::operator+=(const int& x){
   char char_number[50];
   sprintf(char_number, "%d", x);

	int len=strlen(char_number)+_nb_char;
   if(len > 10000){
		return *this;
   }
   if(_tot_space < len+2){
      char* reserve=_data;                   //use a temp pointer to the array
      _tot_space=len+2;
      _data=new char[_tot_space];            //make enought space for catenate
      my_strcpy(_data, reserve);             //copy the array content
      delete[] reserve;                      //delete the old array
   }
   _nb_char=len;
   if (_data) my_strcat(_data, char_number);
   return *this;
}



//------------------------------------------------------------------------------
//addition
my_string
operator+(const my_string& x, const my_string& y){
	int len=x._nb_char+y._nb_char;
   if(len > 10000){
      my_string z;
      return z;
   }
   my_string z(len, 0);
   if (z._data) my_strcpy(z._data, x._data);           //copy x
   if (z._data) my_strcat(z._data, y._data);           //catenate y
   z._nb_char=len;
   return z;
}

my_string
operator+(const my_string& x, const int& y) {  // created Samuel 31.10.2003
   char char_number[50];
   sprintf(char_number, "%d", y);

	int len=x._nb_char+strlen(char_number);
   if(len > 10000){
      my_string z;
      return z;
   }
   my_string z(len, 0);
   if (z._data) my_strcpy(z._data, x._data);                    //copy x
   if (z._data) my_strcat(z._data, char_number);           //catenate y
   z._nb_char=len;
   return z;
}

my_string
operator+(const int& x, const my_string& y) {  // created Samuel 12.08.2004
   char char_number[50];
   sprintf(char_number, "%d", x);

	int len=strlen(char_number)+y._nb_char;
   if(len > 10000){
      my_string z;
      return z;
   }
   my_string z(len, 0);
   if (z._data) my_strcpy(z._data, char_number);        //copy x
   if (z._data) my_strcat(z._data, y._data);           //catenate y
   z._nb_char=len;
   return z;
}

my_string
operator+(const my_string& x, char* y){
	int len=strlen(x.c_str())+strlen(y);
   if(len > 10000){
      my_string z;
      return z;
   }
   my_string z(len+1, 0);
   if (z._data) my_strcpy(z._data, x.c_str());        //copy x
   if (z._data) my_strcat(z._data, y);           //catenate y
   z._nb_char=len;
   return z;
}

my_string
operator+(char* x, const my_string& y){
	int len=strlen(x)+strlen(y.c_str());
   if(len > 10000){
      my_string z;
      return z;
   }
   my_string z(len+1, 0);
   if (z._data) my_strcpy(z._data, x);                   //copy x
   if (z._data) my_strcat(z._data, y.c_str());           //catenate y
   z._nb_char=len;
   return z;
}

my_string
operator+(const my_string& x ,const float& y){   // created Samuel 14.09.2004
   my_string f=y;
   return x+f;
}

//------------------------------------------------------------------------------
// new since 21_3_97 //does not read the delimiter
void
my_string::read_to_delim(istream & is, char delim = '\n'){
	(*this)="";
	char curr_char;
   if(!is) return;
   is.get(curr_char);

   for(;curr_char!=EOF && curr_char!=delim && curr_char!='\n' && curr_char!='\0' && is;){
   	(*this)+=curr_char;             //add curr_char to array
      is.get(curr_char);
   }
}

void
my_string::read_to_delim(istream & is){
	char delim = ' ';
	(*this)="";
	char curr_char;
   if(!is) return;
   is.get(curr_char);

   for(;curr_char!=EOF && curr_char!=delim && curr_char!='\n' && curr_char!='\0'&& curr_char!='\t' && curr_char!='\r' && is;){
   	(*this)+=curr_char;             //add curr_char to array
      is.get(curr_char);
   }
}

// same as read_to_delim but it retuns the delimiter
char
my_string::read_to_delim_plus(istream & is){
	char delim = ' ';
	(*this)="";
	char curr_char;
   if(!is) return ' ';
   is.get(curr_char);

   for(;curr_char!=EOF && curr_char!=delim && curr_char!='\n' && curr_char!='\0'&& curr_char!='\t' && is;){
   	(*this)+=curr_char;             //add curr_char to array
      is.get(curr_char);
   }
   return curr_char;
}



my_string&
my_string::remove( int pos, int n ){  //first position =0 !!!

	if(n<=0 || pos<0) return *this;

	if(n+pos >= _nb_char){  //all is removed
   	_data[pos]='\0';
   }
   else{
      for(int i=pos; i<_nb_char-n; i++){
      	_data[i]=_data[i+n];
      }
      _data[_nb_char-n]='\0';
   }
   _nb_char=strlen(_data);
   return *this;
}

void
my_string::remove_file_extension(){    // Samuel 28-10-2004
   int pos=find_last_of('.');
   _data[pos]='\0';
   _nb_char=strlen(_data);
}
my_string
my_string::get_file_extension(){   //Phäntu 29-08-2008
	int pos=find_last_of('.');
	my_string s;
	for(int i=pos;i<_nb_char;++i){
		s+=_data[i];
	}
	//delete[] line;
	return s;
}


my_string
my_string::extract_sub_str(char delim){ //extracts until the next occurence of delim or end of string //new 11_11_97

	my_string sub_str;
   int nb_removed=0;
   if(!_nb_char) return sub_str;   //return empty string.
   char* buffer=new char[_nb_char+1];   // one char was missing
   //first remove some delimiters at the beginning
   int i=0;
   for(;i<_nb_char && _data[i]==delim; ){
   	i++;
      nb_removed++;
   }
   //then copy the next series of charaters into the buffer until delim occures.
   int pos=0;
   for(;i<_nb_char && _data[i]!=delim; ){
      buffer[pos]=_data[i];
      i++;
      nb_removed++;
      pos++;
   }
   buffer[pos]='\0';
   if(!pos){
      //now remove the first chars.
   	for(i=0;i<_nb_char-nb_removed; i++){
   		_data[i]=_data[i+nb_removed];
   	}
   	_data[_nb_char-nb_removed]='\0';
   	_nb_char-= nb_removed;
		delete[] buffer;
   	return sub_str;   //return empty string.
   }
   sub_str=buffer;   //copy buffer into sub_string;
   //now remove the first chars.
   for(i=0;i<_nb_char-nb_removed; i++){
   	_data[i]=_data[i+nb_removed];
   }
   _data[_nb_char-nb_removed]='\0';
   _nb_char-= nb_removed;
	delete[] buffer;
   return sub_str;
}

my_string
my_string::extract_sub_str_before_ws(){ //extracts until the next occurence of a ws or a tab or end of string
	int a=find_first_of('\t');
	int b=find_first_of(' ');
	if(a==-1) return extract_sub_str(' ');
	if(b==-1) return extract_sub_str('\t');
	if(b<a) return extract_sub_str(' ');
	else return extract_sub_str('\t');
}


/** find the position of the first occurance of the deliminater after position pos */
int
my_string::find_first_of(const char& s, int pos) const{    // changed Samuel 31.10.2003
	for(int i=pos;i<_nb_char;i++){
   	if(_data[i]==s) return i;
   }
   return NPOS;
}

int
my_string::find_first_of(const char& s) const{
	for(int i=0; i<_nb_char; ++i){
      if(_data[i]==s) return i;
   }
   return NPOS;
}

/** find the position of the first occurance of the deliminater after position pos */
int
my_string::find_last_of(const char& s, const int& pos) const{  // created Samuel 23.04.2004
	for(int i=pos-1; i>=0; --i){
      if(_data[i]==s) return i;
   }
   return NPOS;
}

int
my_string::find_last_of(const char& s) const{  // created Samuel 31.10.2003
	for(int i=_nb_char-1; i>=0; --i){
      if(_data[i]==s) return i;
   }
   return NPOS;
}

int
my_string::read_line(istream& is){
	(*this)="";
	char curr_char;
	if(is) is.get(curr_char);
	else return 0;

/*	#ifdef _MAC_FILE_SYST_ //stef_24_6_99  \r is also a control character in MAc text files  //skipp blanks
   	for(;curr_char!=EOF && (curr_char==' ' ||  curr_char=='\t'  ||  curr_char=='\r');){
   	   if(is) is.get(curr_char);
   	   else break;
   	}
   #ielse  //the mac file format has also \r control char
      for(;curr_char!=EOF && (curr_char==' ' ||  curr_char=='\t');)  is.get(curr_char);
 	#endif   */
   //stef_15_7_99 suite to some troubles with LinuX do it any case.
   for(;curr_char!=EOF && (curr_char==' ' ||  curr_char=='\t'  ||  curr_char=='\r');){
   	if(is) is.get(curr_char);
      else break;
   }

   for(;curr_char!=EOF && curr_char!='\n';){
   	  if(is){
   		(*this)+=curr_char;             // if is is bad then stop

         if (curr_char=='\0') break;
   	  }
   	  else{
          ////return 1; (in old version)
   	    return 0;
   	  }
   	  is.get(curr_char);
   }
   return 1;
}

bool
my_string::empty(){
   return !_nb_char;
}

bool
my_string::contains(const char* pat) const{
   int pat_len=strlen(pat);
   if(pat_len==0 || pat_len> _nb_char) return 0;  //pat too long or empty
   for(int i=0; i<_nb_char-pat_len+1; i++){
   	int equal=1;
      for(int j=0;j<pat_len && equal;j++){
      	equal=(_data[j+i]==pat[j]);
      }
      if(equal) return true;
   }
	return false;
}

bool
my_string::contains(char pat) const{
   for(int i=0; i<_nb_char; i++){
      if(_data[i]==pat) return true;
   }
	return false;
}

bool
my_string::contains(const my_string& s) const{
  return contains(s._data);
}

int
my_string::find( const my_string& s ){ //return the position of pattern s if not found return -1
	int pat_len=s._nb_char;
   if(pat_len==0 || pat_len> _nb_char) return NPOS;  //pat to long or empty
   for(int i=0; i<_nb_char-pat_len+1; i++){
   	int equal=1;
      for(int j=0;j<pat_len && equal;j++){
      	equal=(_data[j+i]==(s.c_str())[j]);
      }
      if(equal) return i;         //if here equal is true we found pat and we return position
   }
	return NPOS; //not found then return -1;
}


//------------------------------------------------------------------------------
ostream&
operator<<(ostream& os, const my_string& x) {
	os	<< x._data;
   return os;
}

istream&
operator>>(istream& is, my_string& s){
	char curr_char;
   s="";
   is.get(curr_char);
   for(;curr_char!=EOF && curr_char!='\n' && (curr_char==' ' || curr_char=='\t');){ //skipp blanks before my_string
   	is.get(curr_char);
   }
                                 //skipp also \r char
   for(;curr_char!=EOF && curr_char!='\r' && curr_char!='\n' && curr_char!=' ';){  //stef_16_7_99
   	s+=curr_char;             //add curr_char to array
      is.get(curr_char);
   }
   if(curr_char=='\n') is.putback(curr_char); //it seams that this operator leafs '\n' in stream???
   return is;
}

//------------------------------------------------------------------------------
void
my_string::rm_path(){ // removes path info from a file name

	int pos=0;
	for(int i=_nb_char-1;i>=0;i--){
      #ifdef _MAC_FILE_SYST_
   		if(_data[i]==':'){
      #else
   		if(_data[i]=='/' || _data[i]=='\\'){
      #endif
      	pos=i;
         i=-1;
      }
   }
   pos++; //goto next character after '/'
   int j=0;   //where we write
   for(;pos<_nb_char;){
   	_data[j]=_data[pos];
      j++;pos++;
   }
   _data[j]='\0';  //terminate
   _nb_char=strlen(_data);
    return;
}

//------------------------------------------------------------------------------
void
my_string::extract_path(){
	int pos=0;
	for(int i=_nb_char-1;i>=0;i--){
      #ifdef _MAC_FILE_SYST_
   		if(_data[i]==':'){
      #else
   		if(_data[i]=='/' || _data[i]=='\\'){
      #endif
      	pos=i;
         i=-1;
      }
   }
   if(pos) pos++;    //goto char after / or if path exists :
   _data[pos]='\0';  //terminate
   _nb_char=strlen(_data);
    return;
}

//------------------------------------------------------------------------------
/** returns the filename of a path */
my_string
my_string::extract_filename(){
	int pos=0;
	for(int i=_nb_char-1;i>=0;i--){
      #ifdef _MAC_FILE_SYST_
   		if(_data[i]==':'){
      #else
   		if(_data[i]=='/' || _data[i]=='\\'){
      #endif
      	pos=i;
         i=-1;
      }
   }
   return extract_after(pos);
}

//------------------------------------------------------------------------------
/** returns the part of the string after the char delim, if delim does not exist,
  * it returns the whole string */
my_string
my_string::extract_after(char delim){
	my_string s;
   char* line=new char[_nb_char+1];
	int pos=0;
   for (;pos<_nb_char;){
   	if(_data[pos]==delim) break;
      pos++;
   }
   if(pos==(_nb_char)){
   	s=(*this);
   }
   else{
      pos++;
      int i=0;
      for (;pos<_nb_char;){
      	line[i]=_data[pos];
         pos++;
         i++;
      }
      line[i]='\0';
      s=line;
   }

   s._nb_char=strlen(line);
   delete[] line;
   return s;
}

//------------------------------------------------------------------------------
/** returns the part of the string after the index, if the index is bigger than
  * the string, it returns the whole string */
my_string
my_string::extract_after(int index){  // created Samuel 31.10.2003
	my_string s;
  if(!_nb_char) return s;
  char* line=new char[_nb_char+1];
	int pos=0;
   for (;pos<_nb_char;){
   	if(pos==index) break;
      pos++;
   }
   if(pos==(_nb_char)){
   	s=(*this);
   }
   else{
      pos++;
      int i=0;
      for (;pos<_nb_char;){
      	line[i]=_data[pos];
         pos++;
         i++;
      }
      line[i]='\0';
      s=line;
   }
  s._nb_char=strlen(line);
  delete[] line;
  return s;
}

//------------------------------------------------------------------------------
/** returns the part of the string before the char delim, if delim does not exist,
  * it returns the whole string */
my_string
my_string::extract_before(char delim){
 	my_string s;
   char* line=new char[_nb_char+1];
	int pos=0;
   for (;pos<_nb_char;){
   	if(_data[pos]==delim) break;
      line[pos]=_data[pos];
      pos++;
   }
   line[pos]='\0';
   s=line;
   s._nb_char=strlen(line);
   delete[] line;
   return s;
}

//------------------------------------------------------------------------------
/** returns the part of the string before the index, if the index is bigger than
  * the length of the string, it returns the whole string */
my_string
my_string::extract_before(int index){  // created Samuel 31.10.2003
 	my_string s;
   char* line=new char[_nb_char+1];
	int pos=0;
   for (;pos<_nb_char;){
   	if(pos==index) break;
      line[pos]=_data[pos];
      pos++;
   }
   line[pos]='\0';
   s=line;
   s._nb_char=strlen(line);
   delete[] line;
   return s;
}

//------------------------------------------------------------------------------
/** returns the part of the string before the double slash "//" */
my_string
my_string::extract_before_doubleSlash(){  // created Samuel 18.03.2005
 	my_string s;
   char* line=new char[_nb_char+1];
	int pos=0;
   char delim='/';
   for (;pos<_nb_char;){
   	if(_data[pos]==delim){
         if(_data[pos+1]==delim) break;
      }
      line[pos]=_data[pos];
      pos++;
   }
   line[pos]='\0';
   s=line;
   s._nb_char=strlen(line);
   delete[] line;
   return s;
}

//------------------------------------------------------------------------------
/** returns the part of the string between first and last, if there is an input
  * error (first > last), it returns an empty string (""), if the choosen range
  * is outside the string, it returns the whole string */
my_string
my_string::extract_between(int first, int last){ // created Samuel 31.10.2003

   if (first > last) return "";  // sinple range check
 	my_string s;
   char* line=new char[_nb_char+1];
	int pos=first;
   int i=0;
   for (;pos<_nb_char;){
   	if(pos==last) break;
      line[i]=_data[pos];
      pos++;
      i++;
   }
   line[i]='\0';
   s=line;
   s._nb_char=strlen(line);
   delete[] line;
   return s;
}

//------------------------------------------------------------------------------
void
my_string::remove_blanks(){ //rm all ' ' and '\t'
	char* line=new char[_nb_char+1];
   int new_len=0;
   for(int i=0;i<_nb_char;i++){
   	if(_data[i]!=' ' && _data[i]!='\t'){
      	line[new_len]=_data[i];
         new_len++;
      }
   }
   line[new_len]='\0';
   my_strcpy(_data,line);  //copy the line with no blanks
   _nb_char=new_len;
   delete[] line;


}

//------------------------------------------------------------------------------
void
my_string::trim_blanks(){ //rm all ' ' and '\t' at beginning and end
	char* line=new char[_nb_char+1];
   int new_len=0;
   int i=0;
   while(_data[i]==' ' || _data[i]=='\t') ++i;
   int j=_nb_char-1;
   while(_data[j]==' ' || _data[j]=='\t') --j;
   for(;i<=j;i++){
	   line[new_len]=_data[i];
       new_len++;
   }
   line[new_len]='\0';
   my_strcpy(_data,line);  //copy the line with no blanks
   _nb_char=new_len;
   delete[] line;
}

//------------------------------------------------------------------------------
int //stef_22_12_98  see header for description
my_string::GetNextQuota(ifstream&  ifs){
	char c;
   my_string temp;
   if(!ifs) return 0;
   ifs.get(c);
	for(;ifs && c!='"';){  //goto the first "
   	if(!ifs) return 0;
   	ifs.get(c);
      if(c=='\n') return 0;
   }
   if(!ifs) return 0;
   ifs.get(c);     //get first char after "
   if(c=='\n') return 0;
   for(;ifs && c!='"';){  //until the next "
   	temp=c;
      *this=*this+temp;     //catenate the chars
      if(!ifs) return 0;
   	ifs.get(c);
      if(c=='\n') return 0;
   }

	return 1;
}

//------------------------------------------------------------------------------
//Loro_14_12_98
//A routine to convert backslahes or colon to slashes (for file address),
//depending on file system
void
my_string::convert_to_slash() {
	for(int i=0; i<_nb_char; ++i){
   #ifdef _MAC_FILE_SYST_
   	 if (_data[i]==':') _data[i]='/';   //To_be_checked potential BUG
   #else
   	if (_data[i]=='\\') _data[i]='/';
   #endif
   }
}

//------------------------------------------------------------------------------
void
my_string::convert_to_html_slash() {
	for(int i=0; i<_nb_char; ++i){
   #ifdef _MAC_FILE_SYST_
   	if (_data[i]==':') _data[i]='/';  //To_be_checked potential BUG
   #else
   	if (_data[i]=='\\') _data[i]='/';
   #endif
   }
}

//------------------------------------------------------------------------------
void
my_string::convert_to_underline() {
	for(int i=0; i<_nb_char; ++i){
  /* 	switch (_data[i]) {
			case '\\', ':' , '/':_data[i]='_'; break;
      } */   //not ANSI C++  //stef_7_1_99
     if(_data[i] == '\\' || _data[i] == ':' || _data[i] == '/'){
     	_data[i] = '_';
     }
   }

}

//------------------------------------------------------------------------------
//Loro_19_4_99
my_string
my_string::convert_blank_to_html_blank() {
   my_string html_string("");
	for(int i=0; i<_nb_char; ++i){
     if(_data[i] == ' ') html_string+="%20";
     else html_string+=_data[i];
   }
   return html_string;
}

//------------------------------------------------------------------------------
/** returns the integer and if the string is not an integer it returns NPOS */
int
my_string::toInt(){
   int i=atoi(_data);  // atoi returns 0 if not possible
   if(!i && !this->contains('0')) return NPOS;
   return i;
}

//------------------------------------------------------------------------------
/** returns the double and if the string is not an double it returns 0 */
double my_string::toDouble(){
   char *endptr;
   double d=strtod(_data, &endptr);  // endptr is null if everything is fine
   return d;
}

//------------------------------------------------------------------------------
/** returns true if the string is a number */
bool
my_string::isNumber(){
   //double d=this->toDouble();
   char *endptr;
   double d=strtod(_data, &endptr);
   if(endptr==&_data[_nb_char]) return true;
   return false;
}




