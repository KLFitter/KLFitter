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

#ifndef KLFITTER_RESOLUTIONS_GAUSS_GAUSS_MET_H_
#define KLFITTER_RESOLUTIONS_GAUSS_GAUSS_MET_H_

#include <vector>

#include "KLFitter/Resolutions/Gauss/GaussBase.h"

// ---------------------------------------------------------

namespace KLFitter {
namespace Resolutions {
namespace Gauss {
/**
 * \class KLFitter::ResGauss_MET
 * \brief A class describing a Gaussian resolution, parametrized for MET.
 *
 * This class offers a simple parameterization of a resolution. The
 * parameterization is a Gaussian with a width of a constant times the
 * square root of the true parameter.
 */
class Gauss_MET : public GaussBase {
 public:
  /// The default constructor.
  explicit Gauss_MET(const char * filename);

  /**
   * A constructor.
   * @param parameters The parameters of the parameterization.
   */
  explicit Gauss_MET(std::vector<double> const& parameters);

  /// The (defaulted) destructor.
  ~Gauss_MET();

  /** \name Member functions (Get)  */
  /** @{ */

  /**
   * Return the width of the TF depending on the value of sumET.
   * Used to adjust the range of the fit parameter that correspond to the TF.
   * @param sumet SumET as parameter for the MET TF.
   * @return The width.
   */
  double GetSigma(double sumet) override;

  /** @} */
};
}  // namespace Gauss
}  // namespace Resolutions
}  // namespace KLFitter

#endif  // KLFITTER_RESOLUTIONS_GAUSS_GAUSS_MET_H_
