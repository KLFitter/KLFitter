/*!
 * \class KLFitter::LikelihoodTopLeptonJetsUDSep
 * \brief A class implementing a likelihood for the ttbar lepton+jets channel. 
 * \author Kevin Kr&ouml;ninger
 * \version 1.3
 * \date 03.12.2009
 *
 * This class represents a likelihood for the ttbar into lepton+jets.
 */

// --------------------------------------------------------- 

#ifndef LIKELIHOODTOPLEPTONJETSUDSEP
#define LIKELIHOODTOPLEPTONJETSUDSEP

// --------------------------------------------------------- 

#include "PREPROC.h"
#include "LikelihoodTopLeptonJets.h" 
#include "TLorentzVector.h"
#include "ResolutionBase.h"
#include "TH1F.h"
#include "TH2F.h"

#include <iostream>

// --------------------------------------------------------- 

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{

  class LikelihoodTopLeptonJetsUDSep : public KLFitter::LikelihoodTopLeptonJets
  {
                
  public: 

    /**
     * Enumerate for lJet reweighting methods
     */
    enum LJetSeparationMethod{
      kNone,
      kPermReweight,
      kPermReweight2D
    };
                
    /** \name Constructors and destructors */ 
    /* @{ */ 
                
    /** 
     * The default constructor. 
     */ 
    LikelihoodTopLeptonJetsUDSep(); 
                
    /**
     * The default destructor.
     */
    virtual ~LikelihoodTopLeptonJetsUDSep(); 

    /* @} */

    /**
     * Define the parameters of the fit. 
     */ 
    virtual void DefineParameters();

    /**
     * Return the log of the event probability fof the current
     * combination
     * @return The event probability 
     */ 
    double LogEventProbability(); 

    /**
     * Return the contribution from b tagging to the log of the 
     * event probability for the current combination
     * @return The event probability contribution
     */ 
    double LogEventProbabilityBTag(); 

    /**
     * Return the contribution from pT and b tag weight probability (by LJetSeparationMethod)
     * to the log of the event probability for the current combination
     * @return The event probability contribution
     */ 
    double LogEventProbabilityLJetReweight(); 

    /**
     * Returns the probability of a jet to have the pT of an up type jet.
     * @return The probability.
     */ 
    double UpJetPt(double pt); 

    /**
     * Returns the probability of a jet to have the pT of an down type jet.
     * @return The probability.
     */ 
    double DownJetPt(double pt); 

    /**
     * Returns the probability of a jet to have the pT of an b type jet.
     * @return The probability.
     */ 
    double BJetPt(double pt); 

    /**
     * Returns the probability of a jet to have the tag weight of an up type jet.
     * @return The probability.
     */ 
    double UpJetTagWeight(double tagweight); 

    /**
     * Returns the probability of a jet to have the tag weight of an down type jet.
     * @return The probability.
     */ 
    double DownJetTagWeight(double tagweight); 

    /**
     * Returns the probability of a jet to have the tag weight of an b type jet.
     * @return The probability.
     */ 
    double BJetTagWeight(double tagweight); 
//
    /**
     * Returns the probability of a jet to have the tag weight and pT of an up type jet.
     * @return The probability.
     */ 
    double UpJetProb(double tagweight, double pt); 

    /**
     * Returns the probability of a jet to have the tag weight and pT of an down type jet.
     * @return The probability.
     */ 
    double DownJetProb(double tagweight, double pt); 

    /**
     * Returns the probability of a jet to have the tag weight and pT of an b type jet.
     * @return The probability.
     */ 
    double BJetProb(double tagweight, double pt); 
//

    /**
     * Set histogram for pT distribution of up jets (reco level).
     * @param hist Pointer to histogram. 
     * @return An error flag.
     */ 
    int SetUpJetPtHisto(TH1F* hist) { fUpJetPtHisto = hist; return 1; } 

    /**
     * Set histogram for pT distribution of down jets (reco level).
     * @param hist Pointer to histogram. 
     * @return An error flag.
     */ 
    int SetDownJetPtHisto(TH1F* hist) { fDownJetPtHisto = hist; return 1; } 

    /**
     * Set histogram for pT distribution of b jets (reco level).
     * @param hist Pointer to histogram. 
     * @return An error flag.
     */ 
    int SetBJetPtHisto(TH1F* hist) { fBJetPtHisto = hist; return 1; } 

    /**
     * Set histogram for tag weight distribution of up type jets.
     * @param hist Pointer to histogram. 
     * @return An error flag.
     */ 
    int SetUpJetTagWeightHisto(TH1F* hist) { fUpJetTagWeightHisto = hist; return 1; } 

    /**
     * Set histogram for tag weight distribution of down type jets.
     * @param hist Pointer to histogram. 
     * @return An error flag.
     */ 
    int SetDownJetTagWeightHisto(TH1F* hist) { fDownJetTagWeightHisto = hist; return 1; } 

    /**
     * Set histogram for tag weight distribution of b jets.
     * @param hist Pointer to histogram. 
     * @return An error flag.
     */ 
    int SetBJetTagWeightHisto(TH1F* hist) { fBJetTagWeightHisto = hist; return 1; } 

    /**
     * Set a flag. If flag is true the permutations are reweighted with the pT and tag weight probabilities.
     * @param flag The flag. 
     */ 
    void SetLJetSeparationMethod(KLFitter::LikelihoodTopLeptonJetsUDSep::LJetSeparationMethod flag)
    { fLJetSeparationMethod = flag; }; 

    /**
     * Check if the permutation is LH invariant.
     * @return Permutation of the invariant partner, -1 if there is no one. 
     */
    int LHInvariantPermutationPartner(int iperm, int nperms, int &switchpar1, int &switchpar2);

    /**
     * Set histogram for tag weight distribution of up type jets.
     * @param hist Pointer to histogram. 
     * @return An error flag.
     */ 
    int SetUpJet2DWeightHisto(TH2F* hist) { fUpJet2DWeightHisto = hist; return 1; } 

    /**
     * Set histogram for tag weight distribution of down type jets.
     * @param hist Pointer to histogram. 
     * @return An error flag.
     */ 
    int SetDownJet2DWeightHisto(TH2F* hist) { fDownJet2DWeightHisto = hist; return 1; } 

    /**
     * Set histogram for tag weight distribution of b jets.
     * @param hist Pointer to histogram. 
     * @return An error flag.
     */ 
    int SetBJet2DWeightHisto(TH2F* hist) { fBJet2DWeightHisto = hist; return 1; }


  protected: 

    /** \name Member functions (misc)  */
    /* @{ */


    /**
     * Define the model particles
     * @return An error code.
     */ 
    virtual int DefineModelParticles();

    /**
     * Remove invariant particle permutations.
     * @return An error code. 
     */ 
    int RemoveInvariantParticlePermutations(); 

     /**
     * Remove forbidden particle permutations.
     * @return An error code. 
     */ 
    int RemoveForbiddenParticlePermutations()
    { return 1; }; 


    /**
     * A flag for using an additional reweighting of the permutations with the pT and tag weight probability (default: false);
     */ 
    KLFitter::LikelihoodTopLeptonJetsUDSep::LJetSeparationMethod fLJetSeparationMethod; 

    /**
     * A pointer to the histogram of the up type jet pT distribution. 
     */ 
    TH1F* fUpJetPtHisto; 

    /**
     * A pointer to the histogram of the down type jet pT distribution. 
     */ 
    TH1F* fDownJetPtHisto; 

    /**
     * A pointer to the histogram of the down b jet pT distribution. 
     */ 
    TH1F* fBJetPtHisto; 

    /**
     * A pointer to the histogram of the up quark tag weight distribution. 
     */ 
    TH1F* fUpJetTagWeightHisto; 

    /**
     * A pointer to the histogram of the down quark tag weight distribution. 
     */ 
    TH1F* fDownJetTagWeightHisto; 

    /**
     * A pointer to the histogram of the up b tag weight distribution. 
     */ 
    TH1F* fBJetTagWeightHisto; 

    /**
     * A pointer to the 2d histogram "tag weight vs. pT for upQuarks"
     */ 
    TH2F* fUpJet2DWeightHisto; 

    /**
     * A pointer to the 2d histogram "tag weight vs. pT for downQuarks"
     */ 
    TH2F* fDownJet2DWeightHisto; 

    /**
     * A pointer to the 2d histogram "tag weight vs. pT for bQuarks"
     */ 
    TH2F* fBJet2DWeightHisto; 


    /* @} */

  }; 

} // namespace KLFitter 

// --------------------------------------------------------- 

#endif 

