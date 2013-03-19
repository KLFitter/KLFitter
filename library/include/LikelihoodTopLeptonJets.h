/*!
 * \class KLFitter::LikelihoodTopLeptonJets
 * \brief A class implementing a likelihood for the ttbar lepton+jets channel. 
 * \author Kevin Kr&ouml;ninger
 * \version 1.3
 * \date 03.12.2009
 *
 * This class represents a likelihood for the ttbar into lepton+jets.
 */

// --------------------------------------------------------- 

#ifndef LIKELIHOODTOPLEPTONJETS
#define LIKELIHOODTOPLEPTONJETS

// --------------------------------------------------------- 

#include "PREPROC.h"
#include "LikelihoodBase.h" 
#include "TLorentzVector.h"
#include "ResolutionBase.h"

#include <iostream>

// --------------------------------------------------------- 

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{

  class LikelihoodTopLeptonJets : public KLFitter::LikelihoodBase
  {
                
  public: 
                
    /** \name Constructors and destructors */ 
    /* @{ */ 
                
    /** 
     * The default constructor. 
     */ 
    LikelihoodTopLeptonJets(); 
                
    /**
     * The default destructor.
     */
    virtual ~LikelihoodTopLeptonJets(); 

    /* @} */
    /** \name Member functions (Get)  */
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
    enum Parameters { parBhadE, parBlepE, parLQ1E, parLQ2E, parLepE, parNuPx, parNuPy, parNuPz, parTopM };

    /**
     * Set the values for the missing ET x and y components and the SumET.
     * @param etx missing ET x component.
     * @param ety missing ET y component.
     * @param sumet total scalar ET.
     * @return An error flag.
     */
    int SetET_miss_XY_SumET(double etx, double ety, double sumet);

    /**
     * Set a flag. If flag is true the invariant top quark mass is
     * fixed to the pole mass. 
     * @param flag The flag. 
     */ 
    void SetFlagTopMassFixed(bool flag)
    { fFlagTopMassFixed = flag; }; 

    void SetFlagUseJetMass(bool flag)
    { fFlagUseJetMass = flag; }; 

    void SetFlagGetParSigmasFromTFs(bool flag)
    { fFlagGetParSigmasFromTFs = flag; }; 

    /**
     * Set the type of lepton 
     * @param leptontype The type of lepton: kElectron or kMuon
     */ 
    void SetLeptonType(LeptonType leptontype); 

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
    virtual double LogAPrioriProbability(const std::vector <double> & KLFITTER_UNUSED(parameters))
    { return 0; }; 
                
    /** 
     * The posterior probability definition, overloaded from BCModel. 
     * @param parameters A vector of parameters (double values). 
     * @return The logarithm of the prior probability. 
     */
    virtual double LogLikelihood(const std::vector <double> & parameters); 

    /** 
     * The posterior probability definition, overloaded from BCModel. Split up into several subcomponents 
     * @param parameters A vector of parameters (double values). 
     * @return A vector with the components of the logarithm of the prior probability. Its components are:
     * 0:  TF_bhad
     * 1:  TF_blep
     * 2:  TF_lq1
     * 3:  TF_lq2
     * 4:  TF_lep
     * 5:  TF_METx
     * 6:  TF_METy
     * 7:  BW_Whad
     * 8:  BW_Wlep
     * 9:  BW_Thad
     * 10: BW_Tlep
     */
    virtual std::vector<double> LogLikelihoodComponents(std::vector <double> parameters); 

    /**
     * Get initial values for the parameters.
     * @return vector of initial values.
     */
    virtual std::vector<double> GetInitialParameters();

    /**
     * Get initial values for the parameters with a dummy of "0.0" for the neutrino pz.
     * The decision on the initial value for the neutrino pz then needs to be done in
     * GetInitialParameters().
     * @return vector of initial values.
     */
    virtual std::vector<double> GetInitialParametersWoNeutrinoPz();

    /**
     * Check if there are TF problems.
     * @return Return false if TF problem.
     */
    virtual bool NoTFProblem(std::vector<double> parameters);

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
    virtual int CalculateLorentzVectors(std::vector <double> const& parameters); 

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
     * Remove forbidden particle permutations.
     * @return An error code. 
     */ 
    int RemoveForbiddenParticlePermutations();      

    /**
     * Build the model particles from the best fit parameters.
     * @return An error code.
     */
    int BuildModelParticles();

    /* @} */

  protected: 
                
    /**
     * A flag for using a fixed top mass (true) or not (false).
     */ 
    bool fFlagTopMassFixed; 

    /**
     * A flag for using the measured jet masses (true) instead of
     * parton masses (false);
     */ 
    bool fFlagUseJetMass; 

    /**
     *  Flag for using ResolutionBase::GetSigma() to retrieve the parameter ranges
     */
    bool fFlagGetParSigmasFromTFs;

    /**
     * Return the neutrino pz solutions from the measured values
     * and the W mass.
     * @return A vector with 0, 1 or 2 neutrino pz solutions.
     */
    virtual std::vector<double> GetNeutrinoPzSolutions();

    /**
     * Calculates the neutrino pz solutions from the measured values
     * and the W mass. An additional particle to be added to the 
     * charged lepton may be specified, for example a photon
     * in ttbargamma, which is radiated from the leptonic W
     * or the charged lepton;
     * @param Pointer to a 4-vector of a particle which is added to the charged lepton in the calculation
     * @return A vector with 0, 1 or 2 neutrino pz solutions.
     */
    std::vector<double> CalculateNeutrinoPzSolutions(TLorentzVector * additionalParticle = 0x0);

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
    LeptonType fTypeLepton; 

    /**
     * Global variable for TF problems.
     */
    bool fTFgood;

    /**
     * Save resolution functions since the eta of the partons is not fitted.
     */
    ResolutionBase * fResEnergyBhad;
    ResolutionBase * fResEnergyBlep;
    ResolutionBase * fResEnergyLQ1;
    ResolutionBase * fResEnergyLQ2;
    ResolutionBase * fResLepton;
    ResolutionBase * fResMET;

    /**
     * Save measured particle values for frequent calls
     */
    double bhad_meas_e;
    double bhad_meas_p;
    double bhad_meas_m;
    double bhad_meas_deteta;
    double bhad_meas_eta;
    double bhad_meas_phi;
    double bhad_meas_px;
    double bhad_meas_py;
    double bhad_meas_pz;

    double blep_meas_e;
    double blep_meas_p;
    double blep_meas_m;
    double blep_meas_deteta;
    double blep_meas_eta;
    double blep_meas_phi;
    double blep_meas_px;
    double blep_meas_py;
    double blep_meas_pz;

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
    double bhad_fit_e;
    double bhad_fit_px;
    double bhad_fit_py;
    double bhad_fit_pz;

    double blep_fit_e;
    double blep_fit_px;
    double blep_fit_py;
    double blep_fit_pz;

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

#endif 

