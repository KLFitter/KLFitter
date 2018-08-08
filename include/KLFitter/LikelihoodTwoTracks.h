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

#ifndef LIKELIHOODTWOTRACKS
#define LIKELIHOODTWOTRACKS

// --------------------------------------------------------- 

#include "LikelihoodBase.h" 
#include "TLorentzVector.h"
#include "ResolutionBase.h"

#include <iostream>

// --------------------------------------------------------- 

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{

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
     * The default destructor.
     */
    ~LikelihoodTwoTracks(); 

    /**
     * Calculate 3D Gaussian.
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
    enum Parameters { parPiPlusPhi, parPiPlusTheta, parPiPlusP, parPiMinusPhi, parPiMinusTheta, parPiMinusP, parKShortM };

    
    /* @} */
    /** \name Member functions (misc)  */
    /* @{ */

    /* @} */
    /** \name Member functions (BAT)  */
    /* @{ */

    /**
     * Define the parameters of the fit. 
     */ 
    virtual void DefineParameters();

    /** 
     * The prior probability definition, overloaded from BCModel. 
     * @param parameters A vector of parameters (double values). 
     * @return The logarithm of the prior probability. 
     */
    virtual double LogAPrioriProbability(const std::vector <double> & /*parameters*/)
    { return 0; }; 
                
    /** 
     * The posterior probability definition, overloaded from BCModel. 
     * @param parameters A vector of parameters (double values). 
     * @return The logarithm of the prior probability. 
     */
    virtual double LogLikelihood(const std::vector <double> & parameters); 

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
    virtual std::vector<double> LogLikelihoodComponents(std::vector <double> parameters); 

    /**
     * Get initial values for the parameters.
     * @return vector of initial values.
     */
    virtual std::vector<double> GetInitialParameters();

    /**
     * Dummy function because of KLFitter structure.
     */
    int SetET_miss_XY_SumET(double etx, double ety, double sumet){ return 0;};
    int AdjustParameterRanges(){return 0;};
    int SaveResolutionFunctions(){return 0;};


  protected: 

    /** \name Member functions (misc)  */
    /* @{ */

    /**
     * Update 4-vectors of model particles. 
     * @return An error flag.
     */ 
    virtual int CalculateLorentzVectors(std::vector <double> const& parameters); 

    /**
     * Initialize the likelihood for the event
     */ 
    virtual int Initialize(); 

    /**
     * Define the model particles
     * @return An error code.
     */ 
    virtual int DefineModelParticles();

    /**
     * Remove invariant particle permutations.
     * @return An error code. 
     */ 
    int RemoveInvariantParticlePermutations(); 
    
    /**
     * Build the model particles from the best fit parameters.
     * @return An error code.
     */
    int BuildModelParticles();


    /* @} */

  protected: 
                

    /**
     * Calculates the neutrino pz solutions from the measured values
     * and the W mass. An additional particle to be added to the 
     * charged lepton may be specified, for example a photon
     * in ttbargamma, which is radiated from the leptonic W
     * or the charged lepton;
     * @param Pointer to a 4-vector of a particle which is added to the charged lepton in the calculation
     * @return A vector with 0, 1 or 2 neutrino pz solutions.
     */
    std::vector<double> CalculateNeutrinoPzSolutions(TLorentzVector * additionalParticle = 0x0);

    /**
     * Save permuted particles.
     */
    int SavePermutedParticles();

    double pion_mass;
    double kshort_mass;
    double kshort_width;

    /**
     * Save measured particle values for frequent calls
     */


    double t1_meas_phi;
    double t1_meas_theta;
    double t1_meas_p;

    double t1_meas_sigma00;
    double t1_meas_sigma10;
    double t1_meas_sigma11;
    double t1_meas_sigma20;
    double t1_meas_sigma21;
    double t1_meas_sigma22;

    double t2_meas_phi;
    double t2_meas_theta;
    double t2_meas_p;

    double t2_meas_sigma00;
    double t2_meas_sigma10;
    double t2_meas_sigma11;
    double t2_meas_sigma20;
    double t2_meas_sigma21;
    double t2_meas_sigma22;

    /**
     * Save fit particle values for frequent calls
     */

    double t1_fit_phi;
    double t1_fit_theta;
    double t1_fit_p;
    double t1_fit_m;

    double t2_fit_phi;
    double t2_fit_theta;
    double t2_fit_p;
    double t2_fit_m;

    double ks_fit_m;

  }; 

} // namespace KLFitter 

// --------------------------------------------------------- 

#endif 

