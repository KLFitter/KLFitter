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

#ifndef KLFITTER_RESCRYSTALBALLJETS_H_
#define KLFITTER_RESCRYSTALBALLJETS_H_

#include <iostream>
#include <vector>

#include "KLFitter/Resolutions/ResCrystalBallBase.h"

// ---------------------------------------------------------

namespace KLFitter {
namespace ResCrystalBall {
/**
  * This class offers a simple parameterization of a resolution. The
  * parameterization is a crystal ball with energy dependent
  * parameters.
  */
class ResCrystalBallJets : public ResCrystalBallBase {
 public:
  /**
    * The default constructor.
    */
  explicit ResCrystalBallJets(const char * filename);

  /**
    * A constructor that takes parameters directly, unlike the default 
    * constructor that takes a path to the file with TFs.
    * @param parameters The parameters of the parameterization.
    */
  explicit ResCrystalBallJets(std::vector<double> const& parameters);

  /**
    * The (defaulted) destructor.
    */
  ~ResCrystalBallJets();

  /**
    * Calculate the alpha from the TF parameters and the value of x.
    * @param x The value of x.
    * @return The alpha.
    */
  double GetAlpha(double x) override;

  /**
    * Calculate the n parameter from the TF parameters and the value of x.
    * @param x The value of x.
    * @return The n parameter.
    */
  double GetN(double x) override;

  /**
    * Calculate the sigma from the TF parameters and the value of x.
    * @param x The value of x.
    * @return The sigma.
    */
  double GetSigma(double x) override;

  /**
    * Calculate the mean from the TF parameters and the value of x.
    * @param x The value of x.
    * @return The mean.
    */
  double GetMean(double x) override;

};
}  // namespace ResCrystalBall
}  // namespace KLFitter

#endif  // KLFITTER_RESCRYSTALBALLJETS_H_
