/*!
 * \class KLFitter::InterfaceStructNtup
 * \brief A class for interfacing a Root file. 
 * \author Kevin Kr&ouml;ninger
 * \version 1.3
 * \date 03.12.2009
 *
 * This class interfaces to a Root file which contains events. The
 * events are stored in a Root tree of a certain structure. 
 */

// --------------------------------------------------------- 

#ifndef INTERFACESTRUCTNTUP
#define INTERFACESTRUCTNTUP

// --------------------------------------------------------- 

#include "InterfaceRoot.h" 

// --------------------------------------------------------- 

class TTree;

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{

  class InterfaceStructNtup : public InterfaceRoot
  {
                
  public: 
                
    /** \name Constructors and destructors */ 
    /* @{ */ 
                
    /** 
     * The default constructor. 
     */ 
    InterfaceStructNtup(); 
                
    /**
     * The default destructor.
     */
    virtual ~InterfaceStructNtup(); 

    /* @} */
    /** \name Member functions (Get)  */
    /* @{ */
                
    /**
     * Return the number of events. 
     * @return The number of events
     */ 
    int NEvents(); 

    /**
     * Return the measured missing transverse energy. 
     * @return The missing ET
     */ 
    double ET_miss()
    { return EtMiss / 1000.; }; 

    /**
     * Return the measured missing transverse energy (x component). 
     * @return The missing ET (x component)
     */ 
    double ET_miss_x()
    { return PxMiss / 1000.; }; 

    /**
     * Return the measured missing transverse energy (y component). 
     * @return The missing ET (y component)
     */ 
    double ET_miss_y()
    { return PyMiss / 1000.; }; 

    /* @} */
    /** \name Member functions (Set)  */
    /* @{ */

    /* @} */
    /** \name Member functions (misc)  */
    /* @{ */

    /**
     * Open Root file containing tree.
     * @param filename The filename. 
     * @param opt Options.
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
     * Get event from Root tree. 
     * @param index The event index.
     * @return An error code.
     */ 
    int Event(int index); 

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

    int event_entry; 

    int GoodMuonVec_N;  
    double GoodMuonVec_E[13];  
    double GoodMuonVec_px[13];  
    double GoodMuonVec_py[13];  
    double GoodMuonVec_pz[13];  
    double GoodMuonVec_p_T[13];  
    double GoodMuonVec_eta[13];  
    double GoodMuonVec_phi[13];  
    
    int GoodElectronVec_N;  
    double GoodElectronVec_E[13];  
    double GoodElectronVec_px[13];  
    double GoodElectronVec_py[13];  
    double GoodElectronVec_pz[13];  
    double GoodElectronVec_p_T[13];  
    double GoodElectronVec_eta[13];  
    double GoodElectronVec_phi[13];  
    

    int GoodJetVec_N;  
    double GoodJetVec_E[13];  
    double GoodJetVec_px[13];  
    double GoodJetVec_py[13];  
    double GoodJetVec_pz[13];  
    double GoodJetVec_p_T[13];  
    double GoodJetVec_eta[13];  
    double GoodJetVec_phi[13];  
    double GoodJetVec_weight_SV0[13]; 

    double EtMiss; 
    double PhiMiss; 
    double PxMiss; 
    double PyMiss; 

    
    /* @} */
  }; 

} // namespace KLFitter 

// --------------------------------------------------------- 

#endif 

