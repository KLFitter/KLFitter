#include "InterfaceRoot.h" 

#include <TDirectory.h>
#include <TFile.h>

#include <iostream> 

// --------------------------------------------------------- 
KLFitter::InterfaceRoot::InterfaceRoot() : fRootFile(0) 
{
}

// --------------------------------------------------------- 
KLFitter::InterfaceRoot::~InterfaceRoot()
{
  if (fRootFile)
    delete fRootFile; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceRoot::OpenRootFile(const char * filename, Option_t * opt)
{
  // delete old file 
  if (fRootFile)
    delete fRootFile; 

  // get current directory
  TDirectory * dir = gDirectory; 

  // create new file 
  fRootFile = new TFile(filename, opt); 

  // check if file is open
  if (!fRootFile->IsOpen()) {
    std::cout << "KLFitter::InterfaceRoot::OpenFile(). Could not open Root file." << std::endl; 
    delete fRootFile; 
    return 0; 
  }

  // change directory 
  gDirectory = dir; 

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceRoot::CloseRootFile()
{
  // check if file exists 
  if (!fRootFile) {
    std::cout << "KLFitter::InterfaceRoot::CloseRootFile(). No Root file defined." << std::endl; 
    return 0; 
  } 

  // check if file is open 
  if (!fRootFile->IsOpen()) { 
    std::cout << "KLFitter::InterfaceRoot::CloseRootFile(). Root file not open."<< std::endl; 
    return 0; 
  }
  
  // close file
  fRootFile->Close(); 

  // no error 
  return 1; 
}

// --------------------------------------------------------- 

std::vector<std::string> KLFitter::InterfaceRoot::ReadInputFiles(const char * filename)
{
  // define input file 
  std::ifstream inputfile; 

  // open file 
  inputfile.open(filename); 

  // check if file is open 
  if (!inputfile.is_open())
    {
      std::cout << "KLFitter::InterfaceRoot::ReadInputFiles(). File \"" << filename << "\" not found." << std::endl;              
    }

  // reset parameters 
  std::vector<std::string> filenameVec; 
  std::string name = "empty";
	
	// Read in, in case input files are separeted with "," 
  // read a string via file since long string causes memory error in CINT when it is read via stdin
  std::string argStr;
  std::ifstream ifs(filename);
  std::getline(ifs,argStr);

	// split by ','
  for (size_t i=0,n; i <= argStr.length(); i=n+1)
    {
      n = argStr.find_first_of(',',i);
      if (n == std::string::npos && filenameVec.size()!=0) //end of file
				n = argStr.length();
			else if(n == std::string::npos) //either only one entry or separated by line
        break;
      std::string tmp = argStr.substr(i,n-i);
      filenameVec.push_back(tmp);
    }
	
	//If there are no "," read in line by line
	if(filenameVec.size() ==0){
  	while(!inputfile.eof())
    	{
      	name = ""; 
      	inputfile >> name;
				if(name.size()==0)
					break;
      	filenameVec.push_back(name); 
    	}
	}	
	// close file 
  inputfile.close();    

  // no error 
  return filenameVec; 
}

// -----------------------------------------------------------------------------------------------------------


