/*!
 * \class KLFitter::InterfaceD3PD
 * \brief A class for interfacing a Root file. 
 * \author Kevin Kr&ouml;ninger
 * \version 1.3
 * \date 03.12.2009
 *
 * This class interfaces to a Root file which contains events. The
 * events are stored in a Root tree of a certain structure. 
 */

// --------------------------------------------------------- 

#ifndef INTERFACED3PD
#define INTERFACED3PD

// --------------------------------------------------------- 

#include "InterfaceRoot.h" 
#include <TChain.h>

// --------------------------------------------------------- 

class TTree;

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{

  class InterfaceD3PD : public InterfaceRoot
  {
                
  public: 
                
    /** \name Constructors and destructors */ 
    /* @{ */ 
                
    /** 
     * The default constructor. 
     */ 
    InterfaceD3PD(); 
                
    /**
     * The default destructor.
     */
    virtual ~InterfaceD3PD(); 

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
    { return topMET_et / 1000.; }; 

    /**
     * Return the measured missing transverse energy (x component). 
     * @return The missing ET (x component)
     */ 
    double ET_miss_x()
    { return topMET_etx / 1000.; }; 

    /**
     * Return the measured missing transverse energy (y component). 
     * @return The missing ET (y component)
     */ 
    double ET_miss_y()
    { return topMET_ety / 1000.; }; 

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
     * Open Root files containing trees.
     * @param filenames vector containing the filenames. 
     * @param opt Options.
     * @return An error code.
     */ 
    int OpenRootFiles(std::vector<std::string> filenames, Option_t * opt = "READ");		

    /**
     * Get a tree from Root file and set branch addresses. 
     * @param treename The name of the tree.
     * @return An error code. 
     */ 
    int ConnectTree(const char * treename);

    /**
     * Set branch addresses for a given tree directly, skip OpenRootFile. 
     * @param treename The name of the tree.
     * @return An error code. 
     */ 
    int ConnectTree(TTree * fTree); 
    
    /**
     * Set branch addresses for a given Chain directly, skip OpenRootFile. 
     * @param fChain The name of the Chain.
     * @return An error code. 
     */ 
    int ConnectChain(TChain * fChain); 

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
    /**
     * The Root chain. 
     */      
    TChain * fChain;            
    /** \name Tree variables  */
    /* @{ */

    UInt_t EventNumber; 
    std::vector<double> * mcevt_weight;  
    
    int topMu_n;
    std::vector<int> * topMu_index;
    std::vector<int> * topMu_use;
    std::vector<int> * topMu_inTrigger;   
    std::vector<float> * mu_E;  
    std::vector<float> * mu_px;  
    std::vector<float> * mu_py;  
    std::vector<float> * mu_pz;  
    std::vector<float> * mu_pt;  
    std::vector<float> * mu_eta;
    std::vector<float> * mu_dettheta;    
    std::vector<float> * mu_phi;  

    int topEl_n;
    std::vector<int> * topEl_index;
    std::vector<int> * topEl_use;
    std::vector<int> * topEl_inTrigger;     
    std::vector<float> * el_E;  
    std::vector<float> * el_px;  
    std::vector<float> * el_py;  
    std::vector<float> * el_pz;  
    std::vector<float> * el_pt;  
    std::vector<float> * el_eta;
    std::vector<float> * el_deteta;  
    std::vector<float> * el_phi;  

    int topJet_n;
    std::vector<int> * topJet_index;
    std::vector<int> * topJet_use;
    std::vector<int> * topJet_inTrigger;       
    std::vector<float> * jet_E;  
    //std::vector<float> * jet_px;  
    //std::vector<float> * jet_py;  
    //std::vector<float> * jet_pz;  
    std::vector<float> * jet_pt;  
    std::vector<float> * jet_eta;
    std::vector<float> * jet_deteta;  
    std::vector<float> * jet_phi;  
    std::vector<double> * jet_flavor_weight_SV0; 

    float topMET_et; 
    float topMET_phi; 
    float topMET_etx; 
    float topMET_ety; 

    
    /* @} */
  }; 

} // namespace KLFitter 

// --------------------------------------------------------- 

#endif 

