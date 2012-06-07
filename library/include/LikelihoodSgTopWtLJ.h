/*!
 * \class KLFitter::LikelihoodSgTopWtLJ
 * \brief A class implementing a likelihood for the SgTop Wt -> lepton+jets channel. 
 * \author Kevin Kr&ouml;ninger, Thomas Loddenk&ouml;tter
 * \version 1.3
 * \date 03.12.2009
 *
 * This class represents a likelihood for the single top Wt channel into lepton+jets.
 */

// --------------------------------------------------------- 

#ifndef __LIKELIHOODSGTOPWTLJ__H
#define __LIKELIHOODSGTOPWTLJ__H

// --------------------------------------------------------- 

#include "PREPROC.h"
#include "LikelihoodBase.h"
#include "ResolutionBase.h"

// --------------------------------------------------------- 

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{

  class LikelihoodSgTopWtLJ : public KLFitter::LikelihoodBase
  {
    
  public: 
    
    /** \name Constructors and destructors */ 
    /* @{ */ 
    
    /** 
     * The default constructor. 
     */ 
    LikelihoodSgTopWtLJ(); 
    
    /**
     * The default destructor.
     */
    virtual ~LikelihoodSgTopWtLJ(); 

    /* @} */
    /** \name Member functions (Get)  */

    /**
     * @brief Get a pointer to the lepton in a particles object, depending on value of fTypeLepton
     * @return Pointer to the lepton's TLorentzVector
     */
    TLorentzVector* GetLepton(KLFitter::Particles* particles) ;

    /**
     * @brief Return the lepton type
     * @return fTypeLepton
     */
    int GetLeptonType() ;

    /**
     * @brief Return the top decay hypothesis
     * @return fHadronicTop
     */
    bool GetHadronicTop()
    { return fHadronicTop; };

    /* @{ */

    /* @} */
    /** \name Member functions (Set)  */
    /* @{ */

     /**
     * Enumerator for the lepton type.
     */
    enum LeptonType { kElectron, kMuon };

    /**
     * Enumerator for the parameters.
     */
    enum Parameters { parBE, parLQ1E, parLQ2E, parLepE, parNuPx, parNuPy, parNuPz };

    /**
     * Set the values for the missing ET x and y components and the SumET.
     * @param etx missing ET x component.
     * @param ety missing ET y component.
     * @param sumet total scalar ET. 
     * @return An error flag.
     */
    int SetET_miss_XY_SumET(double etx, double ety, double sumet);

    /**
     * Associate the hadronic leg of the event to the top quark for likelihood calculation.
     */ 
    void SetHadronicTop() 
    { fHadronicTop = kTRUE; }; 

    /**
     * Associate the leptonic leg of the event to the top quark for likelihood calculation.
     */ 
    void SetLeptonicTop() 
    { fHadronicTop = kFALSE; }; 

    void SetFlagUseJetMass(bool flag)
    { fFlagUseJetMass = flag; }; 

    /**
     * Set the type of lepton 
     * @param leptontype The type of lepton: electron(1) or muon (2) 
     */ 
    void SetLeptonType(int leptontype); 

    /* @} */
    /** \name Member functions (misc)  */
    /* @{ */

    /* @} */
    /** \name Member functions (BAT)  */
    /* @{ */

    /**
     * Define the parameters of the fit. 
     */ 
    virtual void DefineParameters();

    /** 
     * The prior probability definition, overloaded from BCModel. 
     * @param parameters A vector of parameters (double values). 
     * @return The logarithm of the prior probability. 
     */
    virtual double LogAPrioriProbability(std::vector <double> KLFITTER_UNUSED(parameters))
    { return 0; }; 
    
    /** 
     * The posterior probability definition, overloaded from BCModel. 
     * @param parameters A vector of parameters (double values). 
     * @return The logarithm of the prior probability. 
     */
    virtual double LogLikelihood(const std::vector<double> & parameters);

    /**
     * Get initial values for the parameters.
     * @return vector of initial values.
     */
    virtual std::vector<double> GetInitialParameters();

    /**
     * Check if there are TF problems.
     * @return Return false if TF problem.
     */
    virtual bool NoTFProblem(const std::vector<double> & parameters);

    /**
     * Return the set of model particles.
     * @return A pointer to the particles.
     */
    virtual KLFitter::Particles* ParticlesModel() {
      BuildModelParticles();
      return fParticlesModel;
    };
    virtual KLFitter::Particles** PParticlesModel() {
      BuildModelParticles();
      return &fParticlesModel;
    };


    /* @} */

  protected: 

    /** \name Member functions (misc)  */
    /* @{ */

    /**
     * Update 4-vectors of model particles. 
     * @return An error flag.
     */ 
    virtual int CalculateLorentzVectors(std::vector <double> parameters);

    /**
     * Initialize the likelihood for the event
     */ 
    virtual int Initialize();

    /**
     * Adjust parameter ranges 
     */ 
    virtual int AdjustParameterRanges();

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
     * Build the model particles from the best fit parameters.
     * @return An error code.
     */
    int BuildModelParticles();

    /* @} */

  protected: 
    
    /**
     * A flag for associating either the hadronic (true) 
     * or the leponic (false) leg of the event to the top quark
     * for likelihood calculation.
     */ 
    bool fHadronicTop; 

    /**
     * A flag for using the measured jet masses (true) instead of
     * parton masses (false);
     */ 
    bool fFlagUseJetMass; 

    /**
     * Calculates the neutrino pz solutions from the measured values
     * and the W mass.
     * @return A vector with 0, 1 or 2 neutrino pz solutions.
     */
    virtual std::vector<double> GetNeutrinoPzSolutions();

     /**
     * Save permuted particles.
     */
    int SavePermutedParticles();

    /**
     * Save resolution functions.
     */
    int SaveResolutionFunctions();

    /**
     * Set model parton mass according to fFlagUseJetMass.
     * @param The jet mass.
     * @param The quark mass.
     * @param The parton px (will be modified, if necessary).
     * @param The parton py (will be modified, if necessary).
     * @param The parton pz (will be modified, if necessary).
     * @param The parton energy (not modified).
     * @return The parton mass.
     */
    inline double SetPartonMass(double jetmass, double quarkmass, double &px, double &py, double &pz, double e) {
      double mass(0.);
      if (fFlagUseJetMass)
        mass = jetmass > 0. ? jetmass : 0.;
      else
        mass = quarkmass;
      double p_orig = sqrt(px*px + py*py + pz*pz);
      double p_newmass = sqrt(e*e - mass*mass);
      double scale = p_newmass / p_orig;
      px *= scale;
      py *= scale;
      pz *= scale;
      return mass;
    }

    /**
     * The values of the x component of the missing ET.
     */
    double ETmiss_x;

    /**
     * The values of the y component of the missing ET.
     */
    double ETmiss_y;

    /**
     * The values of the total scalar ET.
     */
    double SumET;

    /**
     * An index deciding if the event is electron (1) or muon (2) plus
     * jets.
     */ 
    int fTypeLepton; 

    /**
     * Global variable for TF problems.
     */
    bool fTFgood;

    /**
     * Save resolution functions since the eta of the partons is not fitted.
     */
    ResolutionBase * fResEnergyB;
    ResolutionBase * fResEnergyLQ1;
    ResolutionBase * fResEnergyLQ2;
    ResolutionBase * fResLepton;
    ResolutionBase * fResMET;

    /**
     * Save measured particle values for frequent calls
     */
    double b_meas_e;
    double b_meas_p;
    double b_meas_m;
    double b_meas_deteta;
    double b_meas_eta;
    double b_meas_phi;
    double b_meas_px;
    double b_meas_py;
    double b_meas_pz;

    double lq1_meas_e;
    double lq1_meas_p;
    double lq1_meas_m;
    double lq1_meas_deteta;
    double lq1_meas_eta;
    double lq1_meas_phi;
    double lq1_meas_px;
    double lq1_meas_py;
    double lq1_meas_pz;

    double lq2_meas_e;
    double lq2_meas_p;
    double lq2_meas_m;
    double lq2_meas_deteta;
    double lq2_meas_eta;
    double lq2_meas_phi;
    double lq2_meas_px;
    double lq2_meas_py;
    double lq2_meas_pz;

    double lep_meas_e;
    double lep_meas_deteta;
    double lep_meas_sintheta;
    double lep_meas_pt;
    double lep_meas_px;
    double lep_meas_py;
    double lep_meas_pz;


    /**
     * Save fit particle values for frequent calls
    */
    double b_fit_e;
    double b_fit_px;
    double b_fit_py;
    double b_fit_pz;

    double lq1_fit_e;
    double lq1_fit_px;
    double lq1_fit_py;
    double lq1_fit_pz;

    double lq2_fit_e;
    double lq2_fit_px;
    double lq2_fit_py;
    double lq2_fit_pz;

    double lep_fit_e;
    double lep_fit_px;
    double lep_fit_py;
    double lep_fit_pz;

    double nu_fit_e;
    double nu_fit_px;
    double nu_fit_py;
    double nu_fit_pz;

    double whad_fit_m;
    double wlep_fit_m;
    double thad_fit_m;
    double tlep_fit_m;

  }; 

} // namespace KLFitter 

// --------------------------------------------------------- 

#endif // __LIKELIHOODSGTOPWTLJ__H

