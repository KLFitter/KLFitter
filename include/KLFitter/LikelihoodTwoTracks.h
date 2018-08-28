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

#ifndef KLFITTER_LIKELIHOODTWOTRACKS_H_
#define KLFITTER_LIKELIHOODTWOTRACKS_H_

#include <iostream>

#include "KLFitter/LikelihoodBase.h"

class TLorentzVector;

// ---------------------------------------------------------

namespace KLFitter {
class ResolutionBase;

/**
 * This class implements a likelihood for the reconstruction of a
 * particle decaying to two (oppositely) charged particles
 * producing tracks in the detector. Per track, it uses the
 * four-vector (usually constructed from phi, theta, and q/p
 * assuming a charge of +/- 1 and pion/kaon mass) as well as the
 * covariance matrix of the variables phi, theta, and q/p as
 * parameters.
 */
class LikelihoodTwoTracks : public LikelihoodBase {
 public:
  /// The default constructor.
  LikelihoodTwoTracks();

  /// The (defaulted) destructor.
  ~LikelihoodTwoTracks();

  /** Calculate 3D Gaussian.
    * @param x0 First variable point of evalulation
    * @param x1 Second variable point of evaluation
    * @param x2 Third variable point of evaluation
    * @param mean0 First variable mean of distribution
    * @param mean1 Second variable mean of distribution
    * @param mean2 Third variable mean of distribution
    * @param sigma00 First variable variance squared
    * @param sigma10 Covariance of first and second variable
    * @param sigma11 Second variable variance squared
    * @param sigma20 Covariance of first and third variable
    * @param sigma21 Covariance of second and third variable
    * @param sigma22 Third variable variance squared
    * @return Evaluated value of the 3D Gaussian
    */
  double Log3DGaus(double x0, double x1, double x2,
                    double mean0, double mean1, double mean2,
                    double sigma00, double sigma10, double sigma11,
                    double sigma20, double sigma21, double sigma22);

  /// Enumerator for the parameters.
  enum Parameters { parPiPlusPhi,     ///< First (positive) track phi
                    parPiPlusTheta,   ///< First (positive) track theta
                    parPiPlusP,       ///< First (positive) track momentum
                    parPiMinusPhi,    ///< Second (negative) track phi
                    parPiMinusTheta,  ///< Second (negative) track theta
                    parPiMinusP,      ///< Second (negative) track momentum
                    parKShortM        ///< Mass of track particle
  };

  /** \name Member functions (BAT)  */
  /** @{ */

  /// Define the parameters of the fit.
  void DefineParameters() override;

  /**
    * The posterior probability definition, overloaded from BCModel.
    * @param parameters A vector of parameters (double values).
    * @return The logarithm of the prior probability.
    */
  double LogLikelihood(const std::vector<double>& parameters) override;

  /**
    * The posterior probability definition, overloaded from
    * BCModel. Split up into several subcomponents.
    * @param parameters A vector of parameters (double values).
    * @return A vector with the components of the logarithm of the
    * prior probability. Its components are: \n
    *   0) TF_bhad \n
    *   1) TF_blep \n
    *   2) TF_lq1 \n
    *   3) TF_lq2 \n
    *   4) TF_lep \n
    *   5) TF_METx \n
    *   6) TF_METy \n
    *   7) BW_Whad \n
    *   8) BW_Wlep \n
    *   9) BW_Thad \n
    *   10) BW_Tlep
    */
  std::vector<double> LogLikelihoodComponents(std::vector<double> parameters) override;

  /**
    * Get initial values for the parameters.
    * @return vector of initial values.
    */
  std::vector<double> GetInitialParameters() override;

  /// Dummy reimplementation of the base class function.
  int SetET_miss_XY_SumET(double /*etx*/, double /*ety*/, double /*sumet*/) override { return 0; }

  /// Dummy reimplementation from the base class.
  void RequestResolutionFunctions() override { /* do nothing */ }

  /// Dummy reimplementation of the base class function.
  int AdjustParameterRanges() override { return 0; }

  /// Dummy reimplementation of the base class function.
  int SaveResolutionFunctions() override { return 0; }

  /** @} */

 protected:
  /**
    * Update 4-vectors of model particles.
    * @return An error flag.
    */
  int CalculateLorentzVectors(const std::vector<double>& parameters) override;

  /// Initialize the likelihood for the event.
  int Initialize() override;

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

  /// Save permuted particles.
  int SavePermutedParticles() override;

  /** \name Member attributes */
  /** @{ */

  /// Mass of the pion.
  const double m_pion_mass;

  /// Mass of k_short.
  const double m_kshort_mass;

  /// Decay width of k_short.
  const double m_kshort_width;

  /** @} */
  /** \name Member attributes (measured parameters) */
  /** @{ */

  double m_t1_meas_phi;
  double m_t1_meas_theta;
  double m_t1_meas_p;

  double m_t1_meas_sigma00;
  double m_t1_meas_sigma10;
  double m_t1_meas_sigma11;
  double m_t1_meas_sigma20;
  double m_t1_meas_sigma21;
  double m_t1_meas_sigma22;

  double m_t2_meas_phi;
  double m_t2_meas_theta;
  double m_t2_meas_p;

  double m_t2_meas_sigma00;
  double m_t2_meas_sigma10;
  double m_t2_meas_sigma11;
  double m_t2_meas_sigma20;
  double m_t2_meas_sigma21;
  double m_t2_meas_sigma22;

  /** @} */
  /** \name Member attributes (fitted parameters) */
  /** @{ */

  double m_t1_fit_phi;
  double m_t1_fit_theta;
  double m_t1_fit_p;
  double m_t1_fit_m;

  double m_t2_fit_phi;
  double m_t2_fit_theta;
  double m_t2_fit_p;
  double m_t2_fit_m;

  double m_ks_fit_m;

  /** @} */
};
}  // namespace KLFitter

#endif  // KLFITTER_LIKELIHOODTWOTRACKS_H_
