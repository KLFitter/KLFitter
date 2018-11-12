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

#ifndef KLFITTER_RESOLUTIONS_GAUSS_GAUSS_LINEARMET_H_
#define KLFITTER_RESOLUTIONS_GAUSS_GAUSS_LINEARMET_H_

#include <iostream>
#include <vector>

#include "KLFitter/Resolutions/Gauss/GaussBase.h"

// ---------------------------------------------------------

/**
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter {
namespace Resolutions {
namespace Gauss {
/**
  * This class offers a simple parameterization of a resolution. The
  * parameterization is a single Gaussian with energy dependent
  * parameters.
  */
class Gauss_linearMET : public GaussBase {
 public:
  /** \name Constructors and destructors */
  /* @{ */

  /**
    * The default constructor.
    */
  explicit Gauss_linearMET(const char * filename);

  /**
    * A constructor that takes parameters directly, unlike the default 
    * constructor that takes a path to the file with TFs.
    * @param parameters The parameters of the parameterization.
    */
  explicit Gauss_linearMET(std::vector<double> const& parameters);

  ///The (defaulted) destructor.
  ~Gauss_linearMET();

  /* @} */
  /** \name Member functions (Get)  */
  /* @{ */

  /**
    * Calculate the width of the first Gaussian from the TF parameters and the value of x.
    * @param x The value of x.
    * @return The width.
    */
  double GetSigma(double x) override;

  /* @} */
};
}  // namespace Gauss
}  // namespace Resolutions
}  // namespace KLFitter

#endif  // KLFITTER_RESOLUTIONS_GAUSS_GAUSS_LINEARMET_H_
