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

#ifndef KLFITTER_RESDOUBLEGAUSSE_4_H_
#define KLFITTER_RESDOUBLEGAUSSE_4_H_

#include <vector>

#include "KLFitter/ResDoubleGaussBase.h"

// ---------------------------------------------------------

/**
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter {
/**
  * \class KLFitter::ResDoubleGaussE_4
  * \brief A class describing a resolution parameterized with a double Gaussian.
  *
  * This class offers a simple parameterization of a resolution. The
  * parameterization is a double Gaussian with energy dependent
  * parameters.
  */
class ResDoubleGaussE_4 : public ResDoubleGaussBase {
 public:
  /** \name Constructors and destructors */
  /* @{ */

  /**
    * The default constructor.
    */
  explicit ResDoubleGaussE_4(const char * filename);

  /**
    * A constructor.
    * @param parameters The parameters of the parameterization.
    */
  explicit ResDoubleGaussE_4(std::vector<double> const& parameters);

  /**
    * The default destructor.
    */
  virtual ~ResDoubleGaussE_4();

  /* @} */
  /** \name Member functions (Get)  */
  /* @{ */

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

  /* @} */
};
}  // namespace KLFitter

#endif  // KLFITTER_RESDOUBLEGAUSSE_4_H_
