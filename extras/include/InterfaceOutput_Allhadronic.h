/*!
 * \class KLFitter::InterfaceOutput_Allhadronic
 * \brief An output interface 
 *
 * This class is an interface for output data. The data will be
 * written to a Root tree. 
 */

// --------------------------------------------------------- 

#ifndef INTERFACEOUTPUT_ALLHADRONIC
#define INTERFACEOUTPUT_ALLHADRONIC

// --------------------------------------------------------- 

#include "InterfaceRoot.h" 
#include "Particles.h" 

#include <string>

// --------------------------------------------------------- 

class TTree;

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{

  class Fitter;
  class SelectionTool;
  class MatchingTool;

  class InterfaceOutput_Allhadronic : public InterfaceRoot 
  {
                
  public: 
                
    /** \name Constructors and destructors */ 
    /* @{ */ 
                
    /** 
     * The default constructor. 
     */ 
    InterfaceOutput_Allhadronic(); 
                
    /**
     * The default destructor.
     */
    virtual ~InterfaceOutput_Allhadronic(); 

    /* @} */
    /** \name Member functions (Get)  */
    /* @{ */

    /* @} */
    /** \name Member functions (Set)  */
    /* @{ */

    /**
     * Set pointer to fitter.
     * @param fitter A fitter.
     * @return An error code. 
     */ 
    int SetFitter(KLFitter::Fitter * fitter); 

    /**
     * Set pointer to matching tool. 
     * @param matchingtool A matching tool.
     * @return An error code. 
     */ 
    int SetMatchingTool(KLFitter::MatchingTool * matchingtool); 

    /**
     * Set pointer to selection tool. 
     * @param selectiontool A selection tool.
     * @return An error code. 
     */ 
    int SetSelectionTool(KLFitter::SelectionTool * selectiontool); 

    /**
     * Set pointer to a set of truth particles.
     * @param pparticles A set of particles.
     * @return An error code. 
     */ 
    int SetParticlesTruth(KLFitter::Particles ** pparticles)
    { fParticlesTruth = pparticles; return 1; }; 

    /**
     * Set pointer to a set of measured particles.
     * @param pparticles A set of particles.
     * @return An error code. 
     */ 
    int SetParticlesMeasured(KLFitter::Particles ** pparticles)
    { fParticlesMeasured = pparticles; return 1; }; 

    /* @} */
    /** \name Member functions (misc)  */
    /* @{ */

    /**
     * Open Root file containing tree.
     * @param filename The filename. 
     * @param opt options.
     * @return An error code.
     */ 
    virtual int OpenRootFile(const char * filename, Option_t * opt = "RECREATE"); 

    /**
     * Close Root file. 
     * @return An error code.
     */ 
    virtual int CloseRootFile(); 

    /**
     * Creates trees for the truth, model and measured particles. 
     * @return An error code. 
     */ 
    int CreateTrees(); 

    /**
     * Creates a tree for the model particles. 
     * @return An error code. 
     */ 
    int CreateTreeModel(); 

    /**
     * Creates a tree for the measured particles. 
     * @return An error code. 
     */ 
    int CreateTreeMeasured(); 

    /**
     * Creates a tree for the selected particles. 
     * @return An error code. 
     */ 
    int CreateTreeSelected(); 

    /**
     * Creates a tree for the truth particles. 
     * @return An error code. 
     */ 
    int CreateTreeTruth(); 

    /**
     * Creates a tree for the matching information.
     * @return An error code. 
     */ 
    int CreateTreeMatching(); 

    /**
     * Creates a tree for the mapping information.
     * @return An error code. 
     */ 
    int CreateTreeMap(); 

    /**
     * Fills copies of the model particles of the current permutation
     * into the tree variables.
     * @return An error code. 
     */ 
    int FillTreeModelPermutation();

    /**
     * Deletes all but the n best permutation from TreeModel.
     * @param n Number of best permutations which should be kept.
     * @return An error code. 
     */ 
    int TreeModelDeleteAllButNBestPermutations(unsigned int n);

    /**
     * Fills copies of the measured particles into the tree variables.
     * @return An error code.
     */ 
    int FillTreeMeasured();

    /**
     * Fills copies of the measured particles into the tree variables.
     * @return An error code.
     */ 
    int FillTreeSelected();

    /**
     * Fills copies of the truth particles into the tree variables.
     * @return An error code.
     */ 
    int FillTreeTruth();

    /**
     * Fill copies of the matching information in the tree variables.
     * @return An error code.
     */ 
    int FillTreeMatching(); 

    /**
     * Fill copies of the mapping information in the tree variables.
     * @return An error code.
     */ 
    int FillTreeMap(); 

    /**
     * Fill the trees.
     * @return An error code.
     */ 
    int FillTrees(); 

    /**
     * Fill the event weight.
     * @return An error code.
     */
    int SetEventWeight(double weight);

    /**
     * Fill the pileup weight.
     * @return An error code.
     */
    int SetPileupWeight(double weight);          

    /**
     * Fill the photon type.
     * @return An error code.
     */
    int SetPhotonType(bool isNotClassified, bool isRadTopProd, bool isHadTopRadDecay, bool isLepTopRadDecay, bool isHadWRadDecay, bool isLepWRadDecay);

    /**
     * A helper class which removes blanks from a string. 
     * @param str A string. 
     * @return A string. 
     */ 
    std::string ModifyString(std::string str); 

    /* @} */

  protected: 

    /**
     * A pointer to a pointer to a set of truth particles. 
     */ 
    KLFitter::Particles ** fParticlesTruth; 

    /**
     * A pointer to a pointer to a set of model particles. 
     */ 
    KLFitter::Particles ** fParticlesModel; 
                
    /**
     * A pointer to a pointer to a set of measured particles. 
     */ 
    KLFitter::Particles ** fParticlesMeasured; 

    /**
     * A pointer to a pointer to a set of selected particles. 
     */ 
    KLFitter::Particles ** fParticlesSelected; 

    /**
     * A pointer to the matching tool
     */ 
    KLFitter::MatchingTool * fMatchingTool; 

    /**
     * A pointer to the selection tool
     */ 
    KLFitter::SelectionTool * fSelectionTool; 

    /**
     * A pointer to a fitter 
     */ 
    KLFitter::Fitter * fFitter; 

    /**
     * The Root tree for the truth particles
     */ 
    TTree * fTreeTruth; 

    /**
     * The Root tree for the model particles
     */ 
    TTree * fTreeModel; 

    /**
     * The Root tree for the measured particles
     */ 
    TTree * fTreeMeasured; 

    /**
     * The Root tree for the selected particles
     */ 
    TTree * fTreeSelected; 

    /**
     * The Root tree for the matching information
     */ 
    TTree * fTreeMatching; 

    /**
     * The Root tree for mapping the selected to the measured particles 
     */ 
    TTree * fTreeMap; 

  private: 

    /** 
     * A helper method for resizing vectors of vectors.
     * @param v The vector of vectors.
     * @param length The new length.
     */
    template<class type> void Resize(std::vector<std::vector<type> * > * v, unsigned int length);

    /** 
     * A helper method for creating the branch names for TreeMeasured
     * @param The type of the particle
     * @return Particle name - empty string for particles that cannot be measured directly
     */

    std::string TreeMeasuredBranchName(KLFitter::Particles::ParticleType pType);

    /** \name The tree variables */ 
    /* @{ */ 
                
    int fTreeVarEventNumber; 
    int fTreeVarNPermutations;
    int fTreeVarNBTags;
    std::vector<int> * fTreeVarBestPermutation; 
    std::vector<double> * fTreeVarLogLikelihood; 
//    std::vector<std::vector<double> > * fTreeVarLogLikelihoodComponents; 
    std::vector<double> * fTreeVarLogLikelihoodComp_TF_bhad1; 
    std::vector<double> * fTreeVarLogLikelihoodComp_TF_bhad2; 
    std::vector<double> * fTreeVarLogLikelihoodComp_TF_lq1; 
    std::vector<double> * fTreeVarLogLikelihoodComp_TF_lq2; 
    std::vector<double> * fTreeVarLogLikelihoodComp_TF_lq3;
    std::vector<double> * fTreeVarLogLikelihoodComp_TF_lq4;
    std::vector<double> * fTreeVarLogLikelihoodComp_BW_Whad1;
    std::vector<double> * fTreeVarLogLikelihoodComp_BW_Whad2;
    std::vector<double> * fTreeVarLogLikelihoodComp_BW_Thad1;
    std::vector<double> * fTreeVarLogLikelihoodComp_BW_Thad2;

    std::vector<double> * fTreeVarMinuitStatus; 
    std::vector<unsigned int> * fTreeVarConvergenceStatus; 
    std::vector<double> * fTreeVarIntegral;
    std::vector<double> * fTreeVarEventProbability; 

    int fTreeVarNPartonsMeasured; 
    int fTreeVarNElectronsMeasured; 
    int fTreeVarNMuonsMeasured; 
    int fTreeVarNPhotonsMeasured; 

    int fTreeVarNPartonsSelected; 
    int fTreeVarNElectronsSelected; 
    int fTreeVarNMuonsSelected; 
    int fTreeVarNPhotonsSelected; 

    std::vector <std::vector<double> *> * fTreeVarParameters; 
    std::vector <std::vector<double> *> * fTreeVarParameterErrors; 
    std::vector <std::vector<double> *> * fTreeVarModel;
    std::vector <std::vector<int> *> 		* fTreeIntVarModel; 
    std::vector <std::vector<double> *> * fTreeVarMeasured;
    std::vector <std::vector<double> *> * fTreeVarSelected;
    std::vector <std::vector<double> *> * fTreeVarTruth; 
                
    std::vector <int *> * fTreeVarNMatchedPartons; 
    std::vector <int *> * fTreeVarNMatchedElectrons; 
    std::vector <int *> * fTreeVarNMatchedMuons; 
    std::vector <int *> * fTreeVarNMatchedPhotons; 

    std::vector <std::vector<int> *> * fTreeVarMatchedPartons; 
    std::vector <std::vector<int> *> * fTreeVarMatchedElectrons; 
    std::vector <std::vector<int> *> * fTreeVarMatchedMuons; 
    std::vector <std::vector<int> *> * fTreeVarMatchedPhotons; 

    std::vector<int> * fTreeVarMapJets; 
    std::vector<int> * fTreeVarMapElectrons; 
    std::vector<int> * fTreeVarMapMuons; 
    std::vector<int> * fTreeVarMapPhotons; 

    double fEventWeight; 
    double fPileupWeight;  

    bool fIsNotClassified;
    bool fIsRadTopProd;
    bool fIsHadTopRadDecay;
    bool fIsLepTopRadDecay;
    bool fIsHadWRadDecay;
    bool fIsLepWRadDecay;

    /* @} */

  }; 

} // namespace KLFitter 

// --------------------------------------------------------- 

#endif 

