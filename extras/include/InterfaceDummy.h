/*!
 * \class KLFitter::InterfaceRoot
 * \brief A class for interfacing a Root file. 
 * \author Kevin Kr&ouml;ninger
 * \version 1.3
 * \date 03.12.2009
 *
 * This class interfaces to a Root file which contains events. The
 * events are stored in a Root tree of a certain structure. 
 */

// --------------------------------------------------------- 

#ifndef INTERFACEDUMMY
#define INTERFACEDUMMY

// --------------------------------------------------------- 

#include "InterfaceRoot.h" 
#include "Particles.h" 
#include "TFile.h"
#include "TTree.h" 

// --------------------------------------------------------- 

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{

  class InterfaceDummy : public InterfaceRoot
  {
                
  public: 
                
    /** \name Constructors and destructors */ 
    /* @{ */ 
                
    /** 
     * The default constructor. 
     */ 
    InterfaceDummy(); 
                
    /**
     * The default destructor.
     */
    virtual ~InterfaceDummy(); 

    /* @} */
    /** \name Member functions (Get)  */
    /* @{ */

    /**
     * Return the number of events. 
     * @return The number of events
     */ 
    int NEvents(); 

    /**
     * Fill measured particles with data from tree. 
     * @param index The event number. 
     * @return An error code.
     */
    virtual int Event(int index);

    /**
     * Return the measured missing transverse energy. 
     * @return The missing ET
     */ 
    double ET_miss()
    { return MET_Et; }; 

    /**
     * Return the measured missing transverse energy (x component). 
     * @return The missing ET (x component)
     */ 
    double ET_miss_x()
    { return MET_Etx; }; 

    /**
     * Return the measured missing transverse energy (y component). 
     * @return The missing ET (y component)
     */ 
    double ET_miss_y()
    { return MET_Ety; }; 

    /* @} */
        
    /** \name Member functions (Set)  */
    /* @{ */

    /* @} */
    /** \name Member functions (misc)  */
    /* @{ */

    /**
     * Open Root file containing tree.
     * @param filename The filename. 
     * @param opt options.
     * @return An error code.
     */ 
    int OpenRootFile(const char * filename, Option_t * opt = "READ"); 

    /**
     * Get a tree from Root file and set branch addresses. 
     * @param treename The name of the tree.
     * @return An error code. 
     */ 
    int ConnectTree(const char * treename); 

    /**
     * Fill list of particles.
     * @return An error code. 
     */ 
    int FillParticles(); 

    /* @} */

  protected: 

  private: 

    /**
     * The Root tree. 
     */ 
    TTree * fTree; 
                
    /** \name Tree variables  */
    /* @{ */

    double fBhad_E; 
    double fBhad_px; 
    double fBhad_py; 
    double fBhad_pz; 

    double fQup_E; 
    double fQup_px; 
    double fQup_py; 
    double fQup_pz; 

    double fQdown_E; 
    double fQdown_px; 
    double fQdown_py; 
    double fQdown_pz; 

    double fBlep_E; 
    double fBlep_px; 
    double fBlep_py; 
    double fBlep_pz; 

    double fLepton_E; 
    double fLepton_px; 
    double fLepton_py; 
    double fLepton_pz; 

    double fPhoton_E; 
    double fPhoton_px; 
    double fPhoton_py; 
    double fPhoton_pz; 

    double MET_Et; 
    double MET_Etx; 
    double MET_Ety; 
    double MET_Phi; 

    double fTrue_Bhad_E; 
    double fTrue_Bhad_px; 
    double fTrue_Bhad_py; 
    double fTrue_Bhad_pz; 

    double fTrue_Qup_E; 
    double fTrue_Qup_px; 
    double fTrue_Qup_py; 
    double fTrue_Qup_pz; 

    double fTrue_Qdown_E; 
    double fTrue_Qdown_px; 
    double fTrue_Qdown_py; 
    double fTrue_Qdown_pz; 

    double fTrue_Blep_E; 
    double fTrue_Blep_px; 
    double fTrue_Blep_py; 
    double fTrue_Blep_pz; 

    double fTrue_Lepton_E; 
    double fTrue_Lepton_px; 
    double fTrue_Lepton_py; 
    double fTrue_Lepton_pz; 

    double fTrue_Photon_E; 
    double fTrue_Photon_px; 
    double fTrue_Photon_py; 
    double fTrue_Photon_pz; 

    double fTrue_Neutrino_E; 
    double fTrue_Neutrino_px; 
    double fTrue_Neutrino_py; 
    double fTrue_Neutrino_pz; 
    /* @} */
  }; 

} // namespace KLFitter 

// --------------------------------------------------------- 

#endif 

