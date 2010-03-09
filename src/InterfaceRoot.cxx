#include "InterfaceRoot.h" 
#include "TDirectory.h"
#include <iostream> 

// --------------------------------------------------------- 
KLFitter::InterfaceRoot::InterfaceRoot()
{
	fRootFile = 0; 
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
	if (!fRootFile -> IsOpen())
		{
			std::cout << "KLFitter::InterfaceRoot::OpenFile(). Could not open Root file." << std::endl; 
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
	if (!fRootFile)
		{
			std::cout << "KLFitter::InterfaceRoot::CloseRootFile(). No Root file defined." << std::endl; 
			return 0; 
		} 

	// check if file is open 
	if (!fRootFile -> IsOpen())
		{ 
			std::cout << "KLFitter::InterfaceRoot::CloseRootFile(). Root file not open."<< std::endl; 
			return 0; 
		}
	
	fRootFile -> Close(); 

	// no error 
	return 1; 
}

// --------------------------------------------------------- 

