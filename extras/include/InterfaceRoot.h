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

    /** \name Enumerators */
    /* @{ */

    /**
     * Enumerate for signal MC generator
     */
    enum MCGenerator {kHerwig, kAcer};

    /* @} */ 
                
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
     * Return the transverse mass of the W boson
     * @return The MWT
     */ 
    virtual double MWT()
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

    /**
     * Return the measured total scalar transverse energy. 
     * @return The total scalar ET
     */ 
    virtual double Sum_ET()
    { return 0; }; 

    /* @} */
    /** \name Member functions (Set)  */
    /* @{ */

    /**
     * Set FlagWriteSignalMCTruth and chose the corresponding signal MC generator with default kHerwig.
     * @param flag If true, truth particle container is filled.
     * @param mcgen The current signal MC generator for TruthMapping.
     * @return An error code.
     */ 
    int WriteSignalMCTruth(bool flag, KLFitter::InterfaceRoot::MCGenerator mcgen = KLFitter::InterfaceRoot::kHerwig)
    { fFlagWriteSignalMCTruth = flag; fSignalMCGen = mcgen; return 1; };

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
	

    /**
     * Set the proper (global) b-tagging (SV0) cut-value, efficiencies and rejections.
     */
    void SetBtaggingInformation(double cut, double eff, double rej) { fBtagCut = cut; fBtagEff = eff; fBtagRej = rej; };

    /* @} */

  protected: 
    /**
     * A flag for writing Signal MC truth to output.
     */ 
    bool fFlagWriteSignalMCTruth;
    
    /**
     * The current Signal MC generator for TruthMapping.
     */ 
    KLFitter::InterfaceRoot::MCGenerator fSignalMCGen;
    /**
     * The Root file. 
     */ 
    TFile * fRootFile; 

    /**
     * B-tagging stuff
     */
    double fBtagCut;
    double fBtagEff;
    double fBtagRej;
  }; 

} // namespace KLFitter 

// --------------------------------------------------------- 

#endif 

