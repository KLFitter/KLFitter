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

#ifndef KLFITTER_RESOLUTIONS_SINGLEGAUSS_SINGLEGAUSSPT_H_
#define KLFITTER_RESOLUTIONS_SINGLEGAUSS_SINGLEGAUSSPT_H_

#include <iostream>
#include <vector>

#include "KLFitter/Resolutions/SingleGauss/SingleGaussBase.h"

// ---------------------------------------------------------

namespace KLFitter {
namespace Resolutions {
namespace SingleGauss {
/**
  * This class offers a simple parameterization of a resolution. The
  * parameterization is a single Gaussian with pT dependent
  * parameters.
  */
class SingleGaussPt : public SingleGaussBase {
 public:
  /**
    * The default constructor.
    */
  explicit SingleGaussPt(const char * filename);

  /**
    * A constructor that takes parameters directly, unlike the default 
    * constructor that takes a path to the file with TFs.
    * @param parameters The parameters of the parameterization.
    */
  explicit SingleGaussPt(std::vector<double> const& parameters);

  ///The (defaulted) destructor.
  ~SingleGaussPt();

  /**
    * Calculate the mean of the Gaussian from the TF parameters and the value of x.
    * @param x The value of x.
    * @return The mean.
    */
  double GetMean(double x) override;

  /**
    * Calculate the width of the first Gaussian from the TF parameters and the value of x.
    * @param x The value of x.
    * @return The width.
    */
  double GetSigma(double x) override;

};
}  // namespace SingleGauss
}  // namespace Resolutions
}  // namespace KLFitter

#endif  // KLFITTER_RESOLUTIONS_SINGLEGAUSS_SINGLEGAUSSPT_H_
