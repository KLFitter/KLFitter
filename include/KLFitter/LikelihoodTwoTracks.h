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


// --------------------------------------------------------- 

#ifndef KLFITTER_LIKELIHOODTWOTRACKS_H_
#define KLFITTER_LIKELIHOODTWOTRACKS_H_

#include <iostream>

#include "KLFitter/LikelihoodBase.h"

class TLorentzVector;
// --------------------------------------------------------- 

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{
  
  class ResolutionBase;

  class LikelihoodTwoTracks : public KLFitter::LikelihoodBase
  {
                
  public: 
                
    /** \name Constructors and destructors */ 
    /* @{ */ 
                
    /** 
     * The default constructor. 
     */ 
    LikelihoodTwoTracks(); 
                
    /**
    * The (defaulted) destructor.
    */
    ~LikelihoodTwoTracks(); 

    /** Calculate 3D Gaussian.
    * @param x0 First variable point of evalulation
    * @param x1 Second variable point of evaluation
    * @param x3 Third variable point of evaluation
    * @param x0 First variable mean of distribution
    * @param x1 Second variable mean of distribution
    * @param x3 Third variable mean of distribution  
    * @param sigma00 First variable variance squared
    * @param sigma10 Covariance of first and second variable
    * @param sigma11 Second variable variance squared
    * @param sigma20 Covariance of first and third variable
    * @param sigma21 Covariance of second and third variable
    * @param sigma22 Third variable variance squared
    * @return Evaluated value of the 3D Gaussian
    */
    double Log3DGaus(double x0, double x1, double x2, double mean0, double mean1, double mean2, double sigma00, double sigma10, double sigma11, double sigma20, double sigma21, double sigma22);


    /* @} */
    /** \name Member functions (Get)  */
    /* @{ */

    /* @} */
    /** \name Member functions (Set)  */
    /* @{ */

    /**
     * Enumerator for the parameters.
     */
    enum Parameters { parPiPlusPhi,     ///< First (positive) track phi
                      parPiPlusTheta,   ///< First (positive) track theta
                      parPiPlusP,       ///< First (positive) track momentum
                      parPiMinusPhi,    ///< Second (negative) track phi
                      parPiMinusTheta,  ///< Second (negative) track theta
                      parPiMinusP,      ///< Second (negative) track momentum
                      parKShortM        ///< Mass of track particle 
    };

    
    /* @} */
    /** \name Member functions (misc)  */
    /* @{ */

    /* @} */
    /** \name Member functions (BAT)  */
    /* @{ */

    /**
     * Define the parameters of the fit. 
     */ 
    void DefineParameters() override;

    /** 
     * The posterior probability definition, overloaded from BCModel. 
     * @param parameters A vector of parameters (double values). 
     * @return The logarithm of the prior probability. 
     */
    double LogLikelihood(const std::vector <double> & parameters) override; 

    /** 
     * The posterior probability definition, overloaded from BCModel. Split up into several subcomponents 
     * @param parameters A vector of parameters (double values). 
     * @return A vector with the components of the logarithm of the prior probability. Its components are:
     * 0:  TF_bhad
     * 1:  TF_blep
     * 2:  TF_lq1
     * 3:  TF_lq2
     * 4:  TF_lep
     * 5:  TF_METx
     * 6:  TF_METy
     * 7:  BW_Whad
     * 8:  BW_Wlep
     * 9:  BW_Thad
     * 10: BW_Tlep
     */
    std::vector<double> LogLikelihoodComponents(std::vector <double> parameters) override; 

    /**
     * Get initial values for the parameters.
     * @return vector of initial values.
     */
    std::vector<double> GetInitialParameters() override;

    /**
     * Dummy function because of KLFitter structure.
     */
    int SetET_miss_XY_SumET(double etx, double ety, double sumet) override { return 0;};
    int AdjustParameterRanges() override {return 0;};
    int SaveResolutionFunctions() override {return 0;};


  protected: 

    /** \name Member functions (misc)  */
    /* @{ */

    /**
     * Update 4-vectors of model particles. 
     * @return An error flag.
     */ 
    int CalculateLorentzVectors(const std::vector<double>& parameters) override; 

    /**
     * Initialize the likelihood for the event
     */ 
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


    /* @} */

  protected: 

    /**
     * Save permuted particles.
     */
    int SavePermutedParticles() override;

    const double m_pion_mass;
    const double m_kshort_mass;
    const double m_kshort_width;

    /**
     * Save measured particle values for frequent calls
     */


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

    /**
     * Save fit particle values for frequent calls
     */

    double m_t1_fit_phi;
    double m_t1_fit_theta;
    double m_t1_fit_p;
    double m_t1_fit_m;

    double m_t2_fit_phi;
    double m_t2_fit_theta;
    double m_t2_fit_p;
    double m_t2_fit_m;

    double m_ks_fit_m;

  }; 

} // namespace KLFitter 

// --------------------------------------------------------- 

#endif 

