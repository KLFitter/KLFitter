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

#ifndef KLFITTER_RESDOUBLEGAUSSPT_H_
#define KLFITTER_RESDOUBLEGAUSSPT_H_

#include <vector>

#include "KLFitter/Resolutions/ResDoubleGaussBase.h"

// ---------------------------------------------------------

namespace KLFitter {
namespace ResDoubleGauss {
/**
 * \class KLFitter::ResDoubleGaussPt
 * \brief A class describing a resolution parameterized with a double Gaussian.
 *
 * This class offers a simple parameterization of a resolution. The
 * parameterization is a double Gaussian with energy dependent
 * parameters.
 */
class ResDoubleGaussPt : public ResDoubleGaussBase {
 public:
  /// The default constructor.
  explicit ResDoubleGaussPt(const char * filename);

  /**
   * A constructor.
   * @param parameters The parameters of the parameterization.
   */
  explicit ResDoubleGaussPt(std::vector<double> const& parameters);

  /// The (defaulted) destructor.
  ~ResDoubleGaussPt();

  /** \name Member functions (Get)  */
  /** @{ */

  /**
   * Calculate the mean of the first Gaussian from the TF parameters and the value of x.
   * @param x The value of x.
   * @return The width.
   */
  double GetMean1(double x) override;

  /**
   * Calculate the width of the first Gaussian from the TF parameters and the value of x.
   * @param x The value of x.
   * @return The width.
   */
  double GetSigma1(double x) override;

  /**
   * Calculate the amplitude of the second Gaussian from the TF parameters and the value of x.
   * @param x The value of x.
   * @return The width.
   */
  double GetAmplitude2(double x) override;

  /**
   * Calculate the mean of the second Gaussian from the TF parameters and the value of x.
   * @param x The value of x.
   * @return The width.
   */
  double GetMean2(double x) override;

  /**
   * Calculate the width of the sedcond Gaussian from the TF parameters and the value of x.
   * @param x The value of x.
   * @return The width.
   */
  double GetSigma2(double x) override;

  /** @} */
};
}  // namespace ResDoubleGauss
}  // namespace KLFitter

#endif  // KLFITTER_RESDOUBLEGAUSSPT_H_
