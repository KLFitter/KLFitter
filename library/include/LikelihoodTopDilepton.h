/*!
 * \class KLFitter::LikelihoodTopDilepton
 * \brief A class implementing a likelihood for the ttbar dilepton channel. 
 * \author Tamara Vazquez Schröder
 *
 * This class represents a likelihood for the ttbar dilepton channel.
 */

// --------------------------------------------------------- 

#ifndef LIKELIHOODTOPDILEPTON
#define LIKELIHOODTOPDILEPTON

// --------------------------------------------------------- 

#include "PREPROC.h"
#include "LikelihoodBase.h" 
#include "TLorentzVector.h"
#include "ResolutionBase.h"
#include <assert.h>
#include <iostream>
#include <cmath>

#include "BAT/BCModel.h"
#include "BAT/BCH1D.h"

// --------------------------------------------------------- 

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{
  // Neutrino Solution Set
  class NuSolutions {
  public:
    TLorentzVector nu1, nu2;
    int NSolutions;
    NuSolutions():NSolutions(0) {};
      ~NuSolutions() {};
  };

  class LikelihoodTopDilepton : public KLFitter::LikelihoodBase
  {
                
  public: 
                
    /** \name Constructors and destructors */ 
    /* @{ */ 
                
    /** 
     * The default constructor. 
     */ 
    LikelihoodTopDilepton(); 
                
    /**
     * The default destructor.
     */
    virtual ~LikelihoodTopDilepton(); 

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
    //enum Parameters { parTopM=0, parB1E, parB2E, parLep1E, parLep2E, parNuEta, parAntiNuEta /*parNuPx, parNuPy, parNuPz,parTopM*/ };
    enum Parameters { parTopM=0, parB1E, parB2E, parLep1E, parLep2E, parAntiNuEta, parNuEta };

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

    /**
     * Set the neutrino pseudorapidity sigma linear dependency on mtop 
     * according to SM expectations
     */ 
    void SetEtaNuParams(std::vector<double> etanuparam)
    { nueta_params = etanuparam; };

    /**
     * Set the type of lepton 
     * @param leptontype The type of lepton: kElectron or kMuon
     */ 
    void SetLeptonType(LeptonType leptontype_1, LeptonType leptontype_2); 

    /**
     * Set the type of lepton 
     * @param leptontype The type of lepton: electron(1) or muon (2) 
     */ 
    void SetLeptonType(int leptontype_1, int leptontype_2); 

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
     * Define sharp gauss prior for mtop par if mtop fixed
     */ 
    virtual void DefinePrior();
    
    /**
     * Define BCH1D and TH1D histograms to be filled 
     * in MCMCIterationInterface
     */
    void DefineHistograms();
    
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
     * 0:  NuWT
     * 1:  TF_b1
     * 2:  TF_b2
     * 3:  TF_lep1
     * 4:  TF_lep2
     * 5:  AntiNu_Eta
     * 6:  Nu_Eta
     * 7:  Minv(lep,jet)
     */
    virtual std::vector<double> LogLikelihoodComponents(std::vector <double> parameters); 

    /**
     * Get initial values for the parameters.
     * @return vector of initial values.
     */
    virtual std::vector<double> GetInitialParameters();

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

    /**
     * Return Gaussian term for neutrino 
     * pseudorapidity. 
     * @return A double. 
     */ 
    double GaussNuEta(std::vector<double> parameters);
    /**
     * Return Gaussian term for antineutrino 
     * pseudorapidity. 
     * @return A double. 
     */ 
    double GaussAntiNuEta(std::vector<double> parameters);
    /**
     * Return NuWT weight
     * @return A double. 
     */ 
    double CalculateWeight(const std::vector<double> & parameters);
    /**
     * Return NuWT weight for a set of jet1, jet2, lep1, lep2
     * @return A double. 
     */ 
    double CalculateWeightPerm(TLorentzVector * l1, TLorentzVector * l2, TLorentzVector * j1, TLorentzVector * j2, const std::vector<double> & parameters);
    /**
     * Return set of neutrino/antineutrino kinematic solutions 
     * (up to 2)
     * @return A KLFitter::NuSolutions object. 
     */ 
    KLFitter::NuSolutions SolveForNuMom(TLorentzVector * l, TLorentzVector * b, double mtop, double nueta);
    /**
     * Return neutrino weight for a given nu solution and antinu solution
     * @return A double. 
     */
    double neutrino_weight(TLorentzVector nu,TLorentzVector nubar);
    /**
     * Return sum of invariant masses of each (lep,jet) pair, 
     * including a tuning factor alpha.
     * @return A double. 
     */
    double CalculateMLepJet(const std::vector<double> & parameters);
    
    /**
     * Set a flag. If flag is true the sumloglikelihood
     * option is used, instead of the default best-permutation
     * @param flag The flag. 
     */ 
    void SetDoSumLogLik(bool flag)
    { doSumloglik = flag; }; 

    /**
     * TH1D histograms to be filled
     * with functions of interest, e.g.: mttbar,
     * costheta*,etc. for each MCNCiteration
     */ 
    TH1D *  hist_mttbar;
    TH1D *  hist_costheta;
   
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

    /**
     * Calculate other variables out of the KLFitter parameters for each MCMCiteration
     * 
     */
    void MCMCIterationInterface();

    /**
     * Get BAT BCH1D histograms of Mttbar
     * @return BCH1D histograms
     */
    BCH1D * GetHistMttbar() { return fHistMttbar; }

    /**
     * Get BAT BCH1D histograms of CosTheta
     * @return BCH1D histograms
     */
    BCH1D * GetHistCosTheta() {	return fHistCosTheta;  }
    
    /**
     * calculate cos(theta*) for both top and antitop
     */
    std::pair<float, float> CalculateCosTheta(std::vector <TLorentzVector> *particles);

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
     * Lepton 1 Type (electron or muon)
     */ 
    LeptonType fTypeLepton_1; 
    /**
     * Lepton 2 Type (electron or muon) 
     */ 
    LeptonType fTypeLepton_2; 

    /**
     * vector including nu pseudorapidity sigma
     * dependency on mtop ( if sigma=a + b*mtop => 
     * vector[0]=a, vector[1]=b)
     */
    std::vector<double> nueta_params;

    /**
     * A flag for using sumloglikelihood option
     */
    bool doSumloglik;
    
    /**
     * Global variable for TF problems.
     */
    bool fTFgood;

    /**
     * BAT BCH1D Histogram for mttbar
     */
    BCH1D * fHistMttbar;

    /**
     * BAT BCH1D Histogram cos(theta*)
     */
    BCH1D * fHistCosTheta;

    /**
     * BAT BCH1D Histogram for dR(truth top, fit top)
     */
    BCH1D * fHistdRTop;
    /**
     * BAT BCH1D Histogram for dR(truth antitop, fit antitop)
     */
    BCH1D * fHistdRAntiTop;
    /**
     * BAT BCH1D Histogram for dR(truth nu, fit nu)
     */
    BCH1D * fHistdRNu;
    /**
     * BAT BCH1D Histogram for dR(truth antinu, fit antinu)
     */
    BCH1D * fHistdRAntiNu;


    /**
     * Save resolution functions since the eta of the partons is not fitted.
     */
    ResolutionBase * fResEnergyB1;
    ResolutionBase * fResEnergyB2;
    ResolutionBase * fResLepton1;
    ResolutionBase * fResLepton2;
    ResolutionBase * fResMET;

    /**
     * Save measured particle values for frequent calls
     */

    double b1_meas_e;
    double b1_meas_p;
    double b1_meas_m;
    double b1_meas_deteta;
    double b1_meas_eta;
    double b1_meas_phi;
    double b1_meas_px;
    double b1_meas_py;
    double b1_meas_pz;

    double b2_meas_e;
    double b2_meas_p;
    double b2_meas_m;
    double b2_meas_deteta;
    double b2_meas_eta;
    double b2_meas_phi;
    double b2_meas_px;
    double b2_meas_py;
    double b2_meas_pz;

    double lep1_meas_e;
    double lep1_meas_deteta;
    float  lep1_meas_charge;
    double lep1_meas_sintheta;
    double lep1_meas_pt;
    double lep1_meas_px;
    double lep1_meas_py;
    double lep1_meas_pz;

    double lep2_meas_e;
    double lep2_meas_deteta;
    float  lep2_meas_charge;
    double lep2_meas_sintheta;
    double lep2_meas_pt;
    double lep2_meas_px;
    double lep2_meas_py;
    double lep2_meas_pz;

    /**
     * Save fit particle values for frequent calls
     */

    double lep1_fit_e;
    double lep1_fit_px;
    double lep1_fit_py;
    double lep1_fit_pz;

    double lep2_fit_e;
    double lep2_fit_px;
    double lep2_fit_py;
    double lep2_fit_pz;

    double b1_fit_e;
    double b1_fit_px;
    double b1_fit_py;
    double b1_fit_pz;

    double b2_fit_e;
    double b2_fit_px;
    double b2_fit_py;
    double b2_fit_pz;

  }; 

} // namespace KLFitter 

// --------------------------------------------------------- 

#endif 

