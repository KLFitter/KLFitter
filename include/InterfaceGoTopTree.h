/*!
 * \class KLFitter::InterfaceGoTopTree
 * \brief A class for interfacing a Root file. 
 * \author Kevin Kr&ouml;ninger
 * \version 1.3
 * \date 03.12.2009
 *
 * This class interfaces to a Root file which contains events. The
 * events are stored in a Root tree of a certain structure. 
 */

// --------------------------------------------------------- 

#ifndef INTERFACEGoTopTree
#define INTERFACEGoTopTree

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

  class InterfaceGoTopTree : public InterfaceRoot
  {
                
  public: 
                
    /** \name Constructors and destructors */ 
    /* @{ */ 
                
    /** 
     * The default constructor. 
     */ 
    InterfaceGoTopTree(); 
                
    /**
     * The default destructor.
     */
    virtual ~InterfaceGoTopTree(); 

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

    int Event_EventNumber; 

    int Muon_N;  
    std::vector<double> * Muon_E;  
    std::vector<double> * Muon_Px;  
    std::vector<double> * Muon_Py;  
    std::vector<double> * Muon_Pz;  
    std::vector<double> * Muon_Pt;  
    std::vector<double> * Muon_Eta;
    std::vector<double> * Muon_Phi;  
    std::vector<bool> * Muon_IsTopInputs;

    int Electron_N;  
    std::vector<double> * Electron_E;  
    std::vector<double> * Electron_Px;  
    std::vector<double> * Electron_Py;  
    std::vector<double> * Electron_Pz;  
    std::vector<double> * Electron_Pt;  
    std::vector<double> * Electron_Eta;
    std::vector<double> * Electron_DetEta;  
    std::vector<double> * Electron_Phi;  
    std::vector<bool> * Electron_IsTopInputs;

    int Jet_N;  
    std::vector<double> * Jet_E;  
    std::vector<double> * Jet_Px;  
    std::vector<double> * Jet_Py;  
    std::vector<double> * Jet_Pz;  
    std::vector<double> * Jet_Pt;  
    std::vector<double> * Jet_Eta;
    std::vector<double> * Jet_DetEta;  
    std::vector<double> * Jet_Phi;  
    std::vector<double> * Jet_SV0_Weight; 
    std::vector<bool>   * Jet_IsTopInputs;

    int Photon_N;  
    std::vector<double> * Photon_E;  
    std::vector<double> * Photon_Px;  
    std::vector<double> * Photon_Py;  
    std::vector<double> * Photon_Pz;  
    std::vector<double> * Photon_Pt;  
    std::vector<double> * Photon_Eta;  
    std::vector<double> * Photon_Phi;  
    std::vector<bool>   * Photon_IsTopInputs;

    double MET_Et; 
    double MET_Phi; 
    double MET_Etx; 
    double MET_Ety; 

    bool Truth_IsProperMCEvent; 

    int TruthPart_N; 
    std::vector<long> * TruthPart_PDG; 
    std::vector<int> * TruthPart_NParents; 
    std::vector< std::vector<int> > * TruthPart_ParentIdx;
    std::vector<int> * TruthPart_NChildren; 
    std::vector< std::vector<int> > * TruthPart_ChildIdx;
    std::vector<double> * TruthPart_E; 
    std::vector<double> * TruthPart_Px; 
    std::vector<double> * TruthPart_Py; 
    std::vector<double> * TruthPart_Pz; 
    std::vector<double> * TruthPart_Eta; 
    std::vector<double> * TruthPart_Phi; 
    std::vector<double> * TruthPart_Pt; 

    // debugKK
    //          bool Truth_WplusHad; 
    //          bool Truth_WminusHad; 
    char Truth_WplusHad; 
    char Truth_WminusHad; 

    int TruthIdx_b; 
    int TruthIdx_bbar; 
    int TruthIdx_Wminus; 
    int TruthIdx_Wplus; 
    int TruthIdx_lminus; 
    int TruthIdx_lplus; 
    int TruthIdx_n; 
    int TruthIdx_nbar; 
    int TruthIdx_t; 
    int TruthIdx_tbar; 
    int TruthIdx_QfromWminus; 
    int TruthIdx_QfromWplus; 
    int TruthIdx_QbarfromWminus; 
    int TruthIdx_QbarfromWplus; 
    int TruthIdx_photon; 

    /* @} */
  }; 

} // namespace KLFitter 

// --------------------------------------------------------- 

#endif 

