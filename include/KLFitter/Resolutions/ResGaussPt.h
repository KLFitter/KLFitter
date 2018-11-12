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

#ifndef KLFITTER_RESGAUSSPT_H_
#define KLFITTER_RESGAUSSPT_H_

#include <vector>
#include "KLFitter/Resolutions/ResGaussBase.h"

// ---------------------------------------------------------

namespace KLFitter {
namespace ResGauss {
/**
 * \class KLFitter::ResGaussPt
 * \brief A class describing a Gaussian resolution.
 *
 * This class offers a simple parameterization of a resolution. The
 * parameterization is a Gaussian with a width parametrized
 * as a function of the true parameter.
 */
class ResGaussPt : public ResGaussBase {
 public:
  /// The default constructor.
  explicit ResGaussPt(const char * filename);

  /**
   * A constructor.
   * @param sigma The width of the Gaussian.
   */
  explicit ResGaussPt(double sigma);

  /**
   * A constructor.
   * @param sigma The width of the Gaussian.
   */
  explicit ResGaussPt(std::vector<double> const& parameters);

  /// The (defaulted) destructor.
  ~ResGaussPt();

  /** \name Member functions (Get)  */
  /** @{ */

  /**
   * Return the width of the TF depending on the value of energy x.
   * Used to adjust the range of the fit parameter that correspond to the TF.
   * @param par true energy as parameter of the TF.
   * @return The width.
   */
  double GetSigma(double par) override;

  /** @} */
};
}  // namespace ResGauss
}  // namespace KLFitter

#endif  // KLFITTER_ResGaussPt_H_
