/*!
 * \class KLFitter::InterfaceRoot
 * \brief A class for interfacing a Root file. 
 * \author Kevin Kr&ouml;ninger
 * \version 1.3
 * \date 03.12.2009
 *
 * This class interfaces to a Root file.
 */

// --------------------------------------------------------- 

#ifndef INTERFACEROOT
#define INTERFACEROOT

// --------------------------------------------------------- 

#include "PREPROC.h"
#include "InterfaceBase.h" 

#include <TROOT.h>
#include <fstream>

// --------------------------------------------------------- 

class TFile;

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{

  class InterfaceRoot : public InterfaceBase
  {
                
  public: 
                
    /** \name Constructors and destructors */ 
    /* @{ */ 
                
    /** 
     * The default constructor. 
     */ 
    InterfaceRoot(); 
                
    /**
     * The default destructor.
     */
    virtual ~InterfaceRoot(); 

    /* @} */
    /** \name Member functions (Get)  */
    /* @{ */

    /**
     * Fill measured particles with data from tree. 
     * @param index The event number. 
     * @return An error code.
     */
    virtual int Event(int KLFITTER_UNUSED(index))
    { return 0; };  

    /**
     * Return the number of events. 
     * @return The number of events
     */ 
    virtual int NEvents()
    { return 0; }; 

    /**
     * Return the measured missing transverse energy. 
     * @return The missing ET
     */ 
    virtual double ET_miss()
    { return 0; }; 

    /**
     * Return the measured missing transverse energy (x component). 
     * @return The missing ET (x component)
     */ 
    virtual double ET_miss_x()
    { return 0; }; 

    /**
     * Return the measured missing transverse energy (y component). 
     * @return The missing ET (y component)
     */ 
    virtual double ET_miss_y()
    { return 0; }; 

    /* @} */
    /** \name Member functions (Set)  */
    /* @{ */
    /**
     * Set a flag. If flag is true the input is Signal MC.
     * @param flag If true, truth particle container is filled.
     */ 
    void SetFlagIsSignalMC(bool flag)
    { fFlagIsSignalMC = flag; };

    /**
     * Set a flag. Needed for filling the truth particle container.
     * @param flag If flag is true the input is Herwig MC.
     */ 
    void SetFlagIsHerwigMC(bool flag)
    { fFlagIsHerwigMC = flag; };

    /* @} */
    /** \name Member functions (misc)  */
    /* @{ */

    /**
     * Open Root file containing tree.
     * @param filename The filename. 
     * @param opt Options.
     * @return An error code.
     */ 
    virtual int OpenRootFile(const char * filename, Option_t * opt = "READ"); 
		
		 /**
     * Open Root files containing tree.
     * @param filenames The filenames. 
     * @param opt Options.
     * @return An error code.
     */ 
		virtual int OpenRootFiles(std::vector<std::string> filenames, Option_t * opt= "READ")
		{ return 1; }; 
		
    /**
     * Close Root file. 
     * @return An error code.
     */ 
    virtual int CloseRootFile(); 
		/**
     * Read input files from txt-file.
     * @param filename The name of the txt-file.
     * @return An error code. 
     */     
    virtual std::vector<std::string> ReadInputFiles(const char * filename);	

    /* @} */

  protected: 
    /**
     * A flag for using Signal MC as input.
     */ 
    bool fFlagIsSignalMC;

    /**
     * A flag for using Signal MC as input.
     */ 
    bool fFlagIsHerwigMC;

    /**
     * The Root file. 
     */ 
    TFile * fRootFile; 
  }; 

} // namespace KLFitter 

// --------------------------------------------------------- 

#endif 

