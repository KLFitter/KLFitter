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
	 
  // read parameters 
  //for (int i = 0; i < nfiles; ++i)
  while(!inputfile.eof())
    {
      name = ""; 
      inputfile >> name;
			if(name.size()==0)
				break;
      filenameVec.push_back(name); 
    }
	
	// close file 
  inputfile.close();    

  // no error 
  return filenameVec; 
}

// -----------------------------------------------------------------------------------------------------------


