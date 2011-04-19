/*!
 * \class KLFitter::PhotonType
 * \brief A class describing truth photon classification by invariant mass cuts
 * \author Johannes Erdmann
 * \version 1.3
 * \date 01.03.2010
 *
 * This class implements a possibility to classify truth photons by cutting on
 * the invariant masses of the final state particles of hadronic and leptonic
 * top quarks and W bosons including the photon or not in the 4-vector sum.
 */

// --------------------------------------------------------- 

#ifndef PHOTONTYPE
#define PHOTONTYPE

#include "PhysicsConstants.h"
#include <cmath>
#include <TLorentzVector.h>

// --------------------------------------------------------- 

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{

  class Particles;

  class PhotonType
  {
                
  public: 
                
    /** \name Constructors and destructors */ 
    /* @{ */ 
                
    /** 
     * The default constructor. 
     */ 
    PhotonType(); 

    /**
     * The default destructor.
     */
    ~PhotonType() { ; };

    /* @} */
    /** \name Member function (Classification)  */
    /* @{ */

    /**
     * Classify truth photons
     * @return An error flag
     */
    int Classify();

    /* @} */
    /** \name Member functions (Set)  */
    /* @{ */

    /**
     * Set the top mass acceptance
     */
    void SetDeltaTopMass(double dM) { deltaTopMass = dM; };

    /**
     * Set the W mass acceptance
     */
    void SetDeltaWMass(double dM) { deltaWMass = dM; };

    /**
     * Set the truth particles. 
     * @param particles A pointer to a set of particles. 
     * @return An error flag. 
     */ 
    int SetTruthParticles(KLFitter::Particles ** particles) { fTruthParticles = particles; return 1; };

    /**
     * Set the table of physics constants.
     * @return A pointer to the physics constants.
     */ 
    void SetPhysicsConstants(KLFitter::PhysicsConstants * physicsConstants)
    { fPhysicsConstants = physicsConstants; };

    /* @} */
    /** \name Member functions (Get)  */
    /* @{ */

    /**
     * Is not classified
     * @return A boolean
     */
    bool IsNotClassified() { return isNotClassified; };

    /**
     * Is radiative top production
     * @return A boolean
     */
    bool IsRadTopProd() { return isRadTopProd; };

    /**
     * Is radiative decay of hadronic top
     * @return A boolean
     */
    bool IsHadTopRadDecay() { return isHadTopRadDecay; };

    /**
     * Is radiative decay of lepronic top
     * @return A boolean
     */
    bool IsLepTopRadDecay() { return isLepTopRadDecay; };

    /**
     * Is radiative decay of hadronic w
     * @return A boolean
     */
    bool IsHadWRadDecay() { return isHadWRadDecay; };

    /**
     * Is radiative decay of lepronic w
     * @return A boolean
     */
    bool IsLepWRadDecay() { return isLepWRadDecay; };

    /* @} */

  private: 

    /**
     * Make composite particles out of the decay particles
     */
    void MakeCompositeParticles();

    /**
     * The composite particles
     */
    TLorentzVector htop;
    TLorentzVector htop_g;
    TLorentzVector ltop;
    TLorentzVector ltop_g;
    TLorentzVector hadW;
    TLorentzVector hadW_g;
    TLorentzVector lepW;
    TLorentzVector lepW_g;
                
    /**
     * Methods for the actual classification
     * @return A boolean
     */
    bool ClassifyRadTopProd();
    bool ClassifyHadTopRadDecay();
    bool ClassifyLepTopRadDecay();
    bool ClassifyHadWRadDecay();
    bool ClassifyLepWRadDecay();
    /**
     * The top mass acceptance.
     */
    double deltaTopMass;

    /**
     * The W mass acceptance.
     */
    double deltaWMass;

    /**
     * Flag for non-classified events
     */
    bool isNotClassified;

    /**
     * Flag for radiative top production
     */
    bool isRadTopProd;

    /**
     * Flag for radiative decay of hadronic top
     */
    bool isHadTopRadDecay;

    /**
     * Flag for radiative decay of lepronic top
     */
    bool isLepTopRadDecay;

    /**
     * Flag for radiative decay of hadronic w
     */
    bool isHadWRadDecay;

    /**
     * Flag for radiative decay of lepronic w
     */
    bool isLepWRadDecay;

    /**
     * The truth particles
     */
    KLFitter::Particles ** fTruthParticles;

    /**
     * The pointer to the table with the physics constants.
     */
    KLFitter::PhysicsConstants * fPhysicsConstants;

    /**
     * Helper methods for "equal", "lower" and "bigger"
     * w.r.t. to the true particle mass and the delta(M)
     * acceptance.
     */
    bool eqTopM(double m, double acceptance) { return fabs(m - fPhysicsConstants->MassTop()) <  acceptance; }
    bool btTopM(double m, double acceptance) { return m - fPhysicsConstants->MassTop()       >  acceptance; }
    bool ltTopM(double m, double acceptance) { return m - fPhysicsConstants->MassTop()       < -acceptance; }
    bool eqWM(double m, double acceptance)   { return fabs(m - fPhysicsConstants->MassW())   <  acceptance; }
    bool btWM(double m, double acceptance)   { return m - fPhysicsConstants->MassW()         >  acceptance; }
    bool ltWM(double m, double acceptance)   { return m - fPhysicsConstants->MassW()         < -acceptance; }

  }; 
        
} // namespace KLFitter 

// --------------------------------------------------------- 

#endif 
