#define READPARAMETERS_CPP

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <exception>

#include "readparameters.h"

using namespace std;

readparameters::readparameters(const char* filename)
    :filename_(filename) {
  ifstream file(filename);
  if (file.fail()) {
      readparameters_exception excep(string("Could not open ")+string(filename)+string("."));
      throw excep;
  }
  else {
   
      //reading and parsing of parameter file
      //required syntax:
      //
      // [space]name value <value ... value>[#comment]
      //
    string line;
    int linenumber=1;
    typedef string::size_type string_size;
    while (getline(file,line)) {
      string_size i=0;
      string_size size=line.size();
      //whitespace at the begining of the line
      while (i!=size && isspace(line[i])) ++i;
      string_size j=i;
      if (j!=size && line[j]!='#') {//not end of line or a comment

	//keyname
	string key;
	while (j!=size && !isspace(line[j])) ++j;
	if (i!=j) key=line.substr(i,j-i);
	i=j;
	//whitespace after keyname
	while (i!=size && isspace(line[i])) ++i;
	if (i!=size) {
	    vector<string> value;
	    if (line[i]=='"') { //string value
		++i;
		j=i;
		bool inquotes=true;
		string val="";
		for (;i!=size;++i) {
		    if (inquotes) {
			if (line[i]=='"') { //end of a string value
			    inquotes=false;
			    value.push_back(val);
			    val="";
			}
			else if (line[i]=='\\') { //escape character
			    if (i+1==size) {
				readparameters_exception excep((string(filename)+string(" line ")+int2string(linenumber)+string(": Unexpected end of line.")).c_str());
				throw excep;
				
			    }
			    else {
				if (line[i+1]=='\\') {
				    val+='\\';
				    ++i;
				}
				else if (line[i+1]=='"') {
				    val+='"';
				    ++i;
				}
				else val+='\\';
			    }
			}
			else val+=line[i]; //add character or space to string
		    }
		    else {// not in quotes
			if (line[i]=='#') break; //begin of a comment
			else if (line[i]=='"') { //start of new string value
			    inquotes=true;
			}
			else if (!isspace(line[i])) { //syntax error
			    readparameters_exception excep(string(filename)+string(" line ")+int2string(linenumber)+string(": Syntax error."));
		    throw excep;
			}			
		    }

		}//for
		if (inquotes) {
		    readparameters_exception excep(string(filename)+string(" line ")+int2string(linenumber)+string(": Unexpected end of line."));
		    throw excep;
		}
	    }
	    else { //non-string value
		while (i!=size && line[i]!='#') {
		    j=i;
		    while (i!=size && line[i]!='#' && !isspace(line[i]) ) ++i; 
		    if (i!=j) {
			value.push_back(line.substr(j,i-j));
		    }
		    while (i!=size && isspace(line[i]) ) ++i;
		}

	    }
	    tabular[key]=value;	    
	}
	else {
	    //no value specified for key
	    readparameters_exception excep(string(filename)+string(" line ")+int2string(linenumber)+string(": No value specified for '")+key+string("'."));
	    throw excep;
	}
      }//if not a comment
      ++linenumber;
    }//while getline

  }//file.fail
}

readparameters::~readparameters() {
}

template <>
string readparameters::get<string>(string name) const {
    if (exists(name)) return tabular[name][0];
    else {
	readparameters_exception excep(string("Parameter '")+name+string("' does not exist in ")+string(filename_)+string("."));
	throw excep;
    }
}

template <>
vector<string> readparameters::get_vector<string>(string name) const {
    if (exists(name)) return tabular[name];
    else {
	readparameters_exception excep(string("Parameter '")+name+string("' does not exist in ")+string(filename_)+string("."));
	throw excep;
    }
}

void readparameters::print_tabular() {
    cout<<"\ntabular:"<<endl;
    for (map<string,vector<string> >::iterator it=tabular.begin();it!=tabular.end();++it) {
       cout<<"key: "<<(*it).first<<" value: ";
       for (int i=0,N=(*it).second.size();i!=N;++i) cout<<(*it).second[i]<<" ";
       cout<<endl;
    }
    return;
}

const char* readparameters_exception::what() const throw() {
    static string ret;
    ret = string("readparameters exception: ")+what_;
    return ret.c_str();
}

//useful conversion from int to a string
string int2string(int i) {
    stringstream ss;
    string ret;
    ss << i;
    ss >> ret;
    return ret;
}

