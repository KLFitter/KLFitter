/*!
 * \class readparameters
 * \brief A class to read a parameter file.  
 * \author Jörg Meyer
 * \version 0.1
 * \date 12.2008
 * \detail This class opens a parameter file and parses variable names
 * and values. The required syntax of a line of parameter file is:\n
 * variablename value [value]... [#comment]\n
 * examples:\n
 * number 42\n
 * textarray "Hello World" "bye"\n
 * The variable values can be read by the functions get() and get_vector()\n
 * examples:\n
 * int i=rp.get<int>("number");\n
 * std::vector<std::string> textarray=rp.get_vector<std::string>("textarray");\n
 * where rp is an instance of the class readparameters.\n
 */ 

// --------------------------------------------------------- 



#ifndef READPARAMETERS_HPP
#define READPARAMETERS_HPP

#include <vector>
#include <sstream>
#include <map>
#include <string>
#include <exception>

class readparameters_exception: public std::exception {

public:
    readparameters_exception() :
	what_("Unknown exception.") {
    }
    readparameters_exception(std::string c) :
	what_(c) {
    }
    virtual ~readparameters_exception() throw() {
    }
    virtual const char* what() const throw();
private:
    std::string what_;
};


class readparameters {

private:
    /** \name Constructors and destructors */ 
    /* @{ */ 
    
    /** 
     * The default constructor is not ment to be used (yet) and therefore private. 
     */ 
    readparameters() {
    } 
public:
    /** 
     * Constructor taking the name of a parameter file as argument. 
     */ 
    explicit readparameters(const char* filename);
    /** 
     * The default destructor. 
     */ 
    ~readparameters();
    
    /* @} */ 

    /** \name Member functions */ 
    /* @{ */ 
    
    /** 
     * @return Check whether 'name' exists in parameter file.
     */  
    bool exists(std::string name) const {
	return tabular.find(name)!=tabular.end();
    }
    /** 
     * @return Get value of variable of type T called 'name'.
     */     
    template <class T>
    T get(std::string name) const {
	T ret;
	if (exists(name)) {
	    std::stringstream ss;
	    ss << tabular[name][0];
	    ss >> ret;
	}
	else {
	    readparameters_exception excep(std::string("Parameter '")+name+std::string("' does not exist in ")+std::string(filename_)+std::string("."));
	    throw excep;
	}
	return ret;	  
    }
    /** 
     * @return Get value of variable of type std::vector<T> called 'name'.
     */    
    template <class T>
    std::vector<T> get_vector(std::string name) const {
	std::vector<T> ret;
	if (exists(name)) {
	    std::stringstream ss;
	    for (std::vector<std::string>::const_iterator it=tabular[name].begin();it!=tabular[name].end();++it) {
		std::stringstream ss;
		T val;
		ss<<(*it);
		ss>>val;
		ret.push_back(val);
	    }
	}
	else {
	    readparameters_exception excep(std::string("Parameter '")+name+std::string("' does not exist in ")+std::string(filename_)+std::string("."));
	    throw excep;
	}
	return ret;
    }
    
     /** 
     * @return Print content of private member tabular. Only for debugging purposes.
     */ 
    void print_tabular();

     /* @} */
    
private:
    mutable std::map<std::string,std::vector<std::string> > tabular;
    const char* filename_;
};


std::string int2string(int i);

#endif
