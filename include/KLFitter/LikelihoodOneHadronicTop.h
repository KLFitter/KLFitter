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

#ifndef KLFITTER_LIKELIHOODSINGLETOPALLHADRONIC_H_
#define KLFITTER_LIKELIHOODSINGLETOPALLHADRONIC_H_

#include <iostream>
#include <vector>

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
  * \class KLFitter::LikelihoodOneHadronicTop
  * \brief A class implementing a likelihood for the ttbar allhadronic channel.
  *
  * This class represents a likelihood for the ttbar allhadronic channel.
  */
class LikelihoodOneHadronicTop : public KLFitter::LikelihoodBase {
 public:
  /** \name Constructors and destructors */
  /* @{ */

  /**
    * The default constructor.
    */
  LikelihoodOneHadronicTop();

  /**
    * The (defaulted) destructor.
    */
  ~LikelihoodOneHadronicTop();

  /* @} */
  /** \name Member functions (Get)  */
  /* @{ */

  /* @} */
  /** \name Member functions (Set)  */
  /* @{ */

  /**
    * Enumerator for the parameters.
    */
  enum Parameters { parBhad1E, parLQ1E, parLQ2E, parTopM };

  /**
   * Set the values for the missing ET x and y components and the SumET.
   * Reimplemented with dummy implementation to overwrite purely virtual
   * function in base class.
   * @return An error flag.
   */
  int SetET_miss_XY_SumET(double /*etx*/, double /*ety*/, double /*sumet*/) override { return 1; }

  /// Request the necessary resolution functions from the detector.
  void RequestResolutionFunctions() override;

  /**
    * Set a flag. If flag is true the invariant top quark mass is
    * fixed to the pole mass.
    * @param flag The flag.
    */
  void SetFlagTopMassFixed(bool flag) { fFlagTopMassFixed = flag; }

  void SetFlagGetParSigmasFromTFs(bool flag) { fFlagGetParSigmasFromTFs = flag; }

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
  double LogLikelihood(const std::vector <double> &  parameters) override;

  /**
    * The posterior probability definition, overloaded from BCModel. Split up into several subcomponents
    * @param parameters A vector of parameters (double values).
    * @return A vector with the components of the logarithm of the prior probability. Its components are:
    * 0:  TF_bhad1
    * 1:  TF_lq1
    * 2:  TF_lq2
    * 3:  BW_Whad1
    * 4:  BW_Thad1
    */
  std::vector<double> LogLikelihoodComponents(std::vector <double> parameters) override;

  /**
    * Get initial values for the parameters.
    * @return vector of initial values.
    */
  std::vector<double> GetInitialParameters() override;

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
    * Remove forbidden particle permutations - forcing b-jets on the position of a b parton.
    * @return An error code.
    */
  //    int RemoveForbiddenBParticlePermutations();

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
    * Save permuted particles.
    */
  int SavePermutedParticles() override;

  /**
    * Save resolution functions.
    */
  int SaveResolutionFunctions() override;

  /**
    * Save resolution functions since the eta of the partons is not fitted.
    */
  ResolutionBase * fResEnergyBhad1;
  ResolutionBase * fResEnergyLQ1;
  ResolutionBase * fResEnergyLQ2;

  /**
    * Save measured particle values for frequent calls
    */
  double bhad1_meas_e;
  double bhad1_meas_p;
  double bhad1_meas_m;
  double bhad1_meas_deteta;
  double bhad1_meas_eta;
  double bhad1_meas_phi;
  double bhad1_meas_px;
  double bhad1_meas_py;
  double bhad1_meas_pz;

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

  /**
    * Save fit particle values for frequent calls
    */
  double bhad1_fit_e;
  double bhad1_fit_px;
  double bhad1_fit_py;
  double bhad1_fit_pz;

  double lq1_fit_e;
  double lq1_fit_px;
  double lq1_fit_py;
  double lq1_fit_pz;

  double lq2_fit_e;
  double lq2_fit_px;
  double lq2_fit_py;
  double lq2_fit_pz;

  double whad1_fit_m;
  double thad1_fit_m;
};
}  // namespace KLFitter

#endif  // KLFITTER_LIKELIHOODSINGLETOPALLHADRONIC_H_
