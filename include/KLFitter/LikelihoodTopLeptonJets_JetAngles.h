/*
 * Copyright (c) 2009--2018, the KLFitter developer team
 *
 * This file is part of KLFitter.
 *
 * KLFitter is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at
 * your option) any later version.
 *
 * KLFitter is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with KLFitter. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef KLFITTER_LIKELIHOODTOPLEPTONJETS_JETANGLES_H_
#define KLFITTER_LIKELIHOODTOPLEPTONJETS_JETANGLES_H_

#include <iostream>
#include <vector>

class TLorentzVector;

namespace KLFitter {
  class ResolutionBase;
}

#include "KLFitter/LikelihoodBase.h"

// ---------------------------------------------------------

/**
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter {
/**
  * \class KLFitter::LikelihoodTopLeptonJets_JetAngles
  * \brief Add brief description here
  *
  * Add a detailed description here.
  */
class LikelihoodTopLeptonJets_JetAngles : public KLFitter::LikelihoodBase {
 public:
  /** \name Constructors and destructors */
  /* @{ */

  /**
    * The default constructor.
    */
  LikelihoodTopLeptonJets_JetAngles();

  /**
    * The (defaulted) destructor.
    */
  ~LikelihoodTopLeptonJets_JetAngles();

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
  enum Parameters { parBhadE, parBlepE, parLQ1E, parLQ2E, parLepE, parNuPx, parNuPy, parNuPz, parBhadEta, parBlepEta, parLQ1Eta, parLQ2Eta, parBhadPhi, parBlepPhi, parLQ1Phi, parLQ2Phi, parTopM };

  /**
    * Set the values for the missing ET x and y components and the SumET.
    * @param etx missing ET x component.
    * @param ety missing ET y component.
    * @param sumet total scalar ET.
    * @return An error flag.
    */
  int SetET_miss_XY_SumET(double etx, double ety, double sumet) override;

  /**
    * Set a flag. If flag is true the invariant top quark mass is
    * fixed to the pole mass.
    * @param flag The flag.
    */
  void SetFlagTopMassFixed(bool flag) { fFlagTopMassFixed = flag; }

  void SetFlagGetParSigmasFromTFs(bool flag) { fFlagGetParSigmasFromTFs = flag; }

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

  /**
   * Set the detector. This needs a reimplementation, because for
   * this likelihood we need to check whether the chosen detector
   * includes eta and phi transfer functions.
   * @param detector A pointer to a pointer of the detector.
   * @return An error flag
   */
  int SetDetector(KLFitter::DetectorBase** detector) override;

  /* @} */
  /** \name Member functions (misc)  */
  /* @{ */

  /* @} */
  /** \name Member functions (BAT)  */
  /* @{ */

  /**
    * Define the parameters of the fit.
    */
  void DefineParameters() override;

  /**
    * The posterior probability definition, overloaded from BCModel.
    * @param parameters A vector of parameters (double values).
    * @return The logarithm of the prior probability.
    */
  double LogLikelihood(const std::vector <double> & parameters) override;

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
    * 7:  TFeta_bhad
    * 8:  TFeta_blep
    * 9:  TFeta_lq1
    *10:  TFeta_lq2
    *11:  TFphi_bhad
    *12:  TFphi_blep
    *13:  TFphi_lq1
    *14:  TFphi_lq2
    *15:  BW_Whad
    *16:  BW_Wlep
    *17:  BW_Thad
    *18:  BW_Tlep
    */
  std::vector<double> LogLikelihoodComponents(std::vector <double> parameters) override;

  /**
    * Get initial values for the parameters.
    * @return vector of initial values.
    */
  std::vector<double> GetInitialParameters() override;

  /**
    * Get initial values for the parameters with a dummy of "0.0" for the neutrino pz.
    * The decision on the initial value for the neutrino pz then needs to be done in
    * GetInitialParameters().
    * @return vector of initial values.
    */
  std::vector<double> GetInitialParametersWoNeutrinoPz();

  /* @} */

 protected:
  /** \name Member functions (misc)  */
  /* @{ */

  /**
    * Update 4-vectors of model particles.
    * @return An error flag.
    */
  int CalculateLorentzVectors(std::vector <double> const& parameters) override;

  /**
    * Adjust parameter ranges
    */
  int AdjustParameterRanges() override;

  /**
    * Define the model particles
    * @return An error code.
    */
  int DefineModelParticles() override;

  /**
    * Remove invariant particle permutations.
    * @return An error code.
    */
  int RemoveInvariantParticlePermutations() override;

  /**
    * Build the model particles from the best fit parameters.
    * @return An error code.
    */
  int BuildModelParticles() override;

  /* @} */

 protected:
  /**
    * A flag for using a fixed top mass (true) or not (false).
    */
  bool fFlagTopMassFixed;

  /**
    *  Flag for using ResolutionBase::GetSigma() to retrieve the parameter ranges
    */
  bool fFlagGetParSigmasFromTFs;

  /**
    * Return the neutrino pz solutions from the measured values
    * and the W mass.
    * @return A vector with 0, 1 or 2 neutrino pz solutions.
    */
  std::vector<double> GetNeutrinoPzSolutions();

  /**
    * Calculates the neutrino pz solutions from the measured values
    * and the W mass. An additional particle to be added to the
    * charged lepton may be specified, for example a photon
    * in ttbargamma, which is radiated from the leptonic W
    * or the charged lepton;
    * @param additionalParticle Pointer to a 4-vector of a particle which is
    * added to the charged lepton in the calculation
    * @return A vector with 0, 1 or 2 neutrino pz solutions.
    */
  std::vector<double> CalculateNeutrinoPzSolutions(TLorentzVector * additionalParticle = 0x0);

  /**
    * A helper for calculating delta(phi1, phi2)
    */
  double diffPhi(double phi1, double phi2);

  /**
    * Save permuted particles.
    */
  int SavePermutedParticles() override;

  /**
    * Save resolution functions.
    */
  int SaveResolutionFunctions() override;

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
    * Save resolution functions for all particles where the eta is not fitted.
    */
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
}  // namespace KLFitter

#endif  // KLFITTER_LIKELIHOODTOPLEPTONJETS_JETANGLES_H_
