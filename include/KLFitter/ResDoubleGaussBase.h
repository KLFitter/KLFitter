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

#ifndef RESDOUBLEGAUSSBASE
#define RESDOUBLEGAUSSBASE

#include <vector>
#include "ResolutionBase.h" 
#include <iostream>

// --------------------------------------------------------- 

/**
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{

  /**
   * \class KLFitter::ResDoubleGaussBase
   * \brief A class describing a resolution parameterized with a double Gaussian. 
   * \author Kevin Kr&ouml;ninger
   *
   * This class offers a simple parameterization of a resolution. The
   * parameterization is a double Gaussian with energy dependent
   * parameters.
   */
  class ResDoubleGaussBase : public ResolutionBase
  {
                
  public: 
                
    /** \name Constructors and destructors */ 
    /* @{ */ 
                
    /** 
     * The default constructor. 
     */ 
    ResDoubleGaussBase(const char * filename); 

    /**
     * A constructor.
     * @param parameters The parameters of the parameterization. 
     */
    ResDoubleGaussBase(std::vector<double> const& parameters);

    /**
     * The default destructor.
     */
    virtual ~ResDoubleGaussBase(); 

    /* @} */
    /** \name Member functions (Get)  */
    /* @{ */

    /**
     * Calculate the mean of the first Gaussian from the TF parameters and the value of x.
     * @param x The value of x.
     * @return The width. 
     */ 
    virtual double GetMean1(double x) = 0;

    /**
     * Calculate the width of the first Gaussian from the TF parameters and the value of x.
     * @param x The value of x.
     * @return The width. 
     */ 
    virtual double GetSigma1(double x) = 0;

    /**
     * Calculate the amplitude of the second Gaussian from the TF parameters and the value of x.
     * @param x The value of x.
     * @return The width. 
     */ 
    virtual double GetAmplitude2(double x) = 0;

    /**
     * Calculate the mean of the second Gaussian from the TF parameters and the value of x.
     * @param x The value of x.
     * @return The width. 
     */ 
    virtual double GetMean2(double x) = 0;

    /**
     * Calculate the width of the sedcond Gaussian from the TF parameters and the value of x.
     * @param x The value of x.
     * @return The width. 
     */ 
    virtual double GetSigma2(double x) = 0;

    /**
     * Return the approximate width of the TF depending on the measured value of x.
     * Used to adjust the range of the fit parameter that correspond to the TF.
     * @param xmeas The measured value of x.
     * @return The width. 
     */ 
    virtual double GetSigma(double xmeas);

    /**
     * Return the probability of the true value of x given the
     * measured value, xmeas.
     * @param x The true value of x.
     * @param xmeas The measured value of x.
     * @param good False if problem with TF.
     * @return The probability. 
     */ 
    virtual double p(double x, double xmeas, bool &good); 

    /**
     * Return the probability of the true value of x given the
     * measured value, xmeas.
     * @param x The true value of x.
     * @param xmeas The measured value of x.
     * @param good False if problem with TF.
     * @param par Optional additional parameter (SumET in case of MET TF).
     * @return The probability. 
     */ 
    virtual double p(double x, double xmeas, bool &good, double par)
    { good = true; return 0; } 

    /* @} */
    /** \name Member functions (Set)  */
    /* @{ */
        
    /* @} */
    /** \name Member functions (misc)  */
    /* @{ */
                
    /* @} */

    /**
     * Sanity check for double gaussian parameters p2, p3 and p5 (1st sigma, scale and 2nd sigma).
     * @param sigma1 (the 1st sigma).
     * @param amplitude2 (the scale parameter).
     * @param sigma2 (the 2nd sigma).
     * @return False if problem with TF.
     */
    inline static bool CheckDoubleGaussianSanity(double &sigma1, double &amplitude2, double &sigma2) {
      if (amplitude2 < 0.) amplitude2 = 0.;      
      if (sigma1 < 0.) {
//        std::cout << "KLFitter::ResDoubleGauss::CheckDoubleGaussianSanity() ERROR IN TRANSFERFUNCTIONS the sigma of the 1st Gaussian is < 0  -  FIT RESULT MAY NOT BE RELIABLE" << std::endl;
//        std::cout << "--> Fitter is in a bad condition! Please check your input files (E out of validation scope?)." << std::endl;
        sigma1 = 0.00000001;
        return false;
	}
      if (sigma2 < 0.) {
//        std::cout << "KLFitter::ResDoubleGauss::CheckDoubleGaussianSanity() ERROR IN TRANSFERFUNCTIONS the sigma of the 2nd Gaussian is < 0  -  FIT RESULT MAY NOT BE RELIABLE" << std::endl;
//        std::cout << "--> Fitter is in a bad condition! Please check your input files (E out of validation scope?)." << std::endl;
        sigma2 = 0.000000001;
        return false;
	}

      return true;
    }

  private: 

  }; 
        
} // namespace KLFitter 

// --------------------------------------------------------- 

#endif 

