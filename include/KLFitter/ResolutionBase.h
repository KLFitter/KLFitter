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

#ifndef KLFITTER_RESOLUTIONBASE_H_
#define KLFITTER_RESOLUTIONBASE_H_

#include <vector>

// ---------------------------------------------------------

/**
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter {
/**
  * \class KLFitter::ResolutionBase
  * \brief A base class for describing resolutions.
  *
  * This base class can be used to decribe resolutions.
  */
class ResolutionBase {
 public:
  /** \name Constructors and destructors */
  /* @{ */

  /**
    * The default constructor.
    * @param npar The number of parameters needed for the parameterization.
    */
  explicit ResolutionBase(int npar);

  /**
    * A constructor.
    * @param parameters A vector with the parameters.
    */
  explicit ResolutionBase(std::vector<double> parameters);

  /**
    * The (defaulted) destructor.
    */
  virtual ~ResolutionBase();

  /* @} */
  /** \name Member functions (Get)  */
  /* @{ */

  /**
    * Return the (approximate) width of the TF depending on the relevant parameter
    * (mostly xmeas, but sumET in case of ResGauss_MET).
    * Use to adjust the range of the fit parameter that correspond to the TF.
    * @param par Parameter on which the width depends
    * @return The width.
    */
  virtual double GetSigma(double par) = 0;

  /**
    * Return the probability of the true value of x given the
    * measured value, xmeas.
    * @param x The true value of x.
    * @param xmeas The measured value of x.
    * @param good False if problem with TF.
    * @param par Optional additional parameter (SumET in case of MET TF).
    * @return Logarithm of the probability.
    */
  virtual double logp(double /*x*/, double /*xmeas*/, bool *good, double /*par*/ = 0) { *good = true; return 0; }

  /**
    * Return a parameter of the parameterization.
    * @param index The parameter index.
    * @param par The parameter value.
    * @return An error flag.
    */
  int Par(int index, double *par);

  /* @} */
  /** \name Member functions (Set)  */
  /* @{ */

  /**
    * Set a parameter value.
    * @param index The parameter index.
    * @param value The parameter value.
    * @return An error code.
    */
  int SetPar(int index, double value);

  /**
    * Set all parameter values.
    * @param parameters A vector of parameters.
    * @return An error code.
    */
  int SetPar(std::vector<double> parameters);

  /* @} */
  /** \name Member functions (misc)  */
  /* @{ */

  /**
    * Read parameter values from ASCII file.
    * @param filename The name of the file.
    * @param nparameters The number of parameters.
    * @return An error code.
    */
  int ReadParameters(const char * filename, int nparameters);

  /**
    * Return a status code.
    * @return A status code (1: ok, 0: error).
    */
  int Status() { return fStatus; }

  /* @} */

 protected:
  /// The number of parameters.
  int fNParameters;

  /// The parameter values.
  std::vector <double> fParameters;

  /// The status of this class (1: ok, 0: error).
  int fStatus;
};
}  // namespace KLFitter

#endif  // KLFITTER_RESOLUTIONBASE_H_
