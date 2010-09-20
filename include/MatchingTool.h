/*
 * \class KLFitter::MatchingTool
 * \brief A class for matching truth to measured Quantities. 
 * \author Stefan Guindon
 * \version 1.3
 * \date 03.12.2009
 *
 * This class matches truth to measured sets of particles.
 * 
 */


// --------------------------------------------------------- 

#ifndef MATCHINGTOOL
#define MATCHINGTOOL

// --------------------------------------------------------- 

#include "Particles.h" 
 
// --------------------------------------------------------- 

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{

  class MatchingTool  
  {

  public: 

    /** \name Constructors and destructors */ 
    /* @{ */ 
        
    /** 
     * The default constructor. 
     * @param particles A pointer to a set of particles 
     * @param particlestruth A pointer to a set of truth particles 
     */ 
    MatchingTool(KLFitter::Particles ** particles, KLFitter::Particles ** particlestruth);

    /**
     * The default destructor.
     */
    virtual ~MatchingTool();

    /* @} */
    /** \name Member functions (Get)  */
    /* @{ */

    /**
     * Return the matching status of a truth and a reco particle. 
     * @param indextruth The index of the truth particle. 
     * @param indexreco The index of the reco particle
     * @param ptype The type of the particle. 
     * @return The matching status of a truth and a reec particle. 
     */ 
    int MatchingStatus(int indextruth, int indexreco, KLFitter::Particles::ParticleType ptype); 

    /**
     * Return the number of matches of a truth particle.
     * @param index  The index of the particle.
     * @param ptype The type of the particle. 
     * @return The number of matches of a truth particle. 
     */
    int NMatchedTruth(int index, KLFitter::Particles::ParticleType ptype);

    /**
     * Return the number of matches of a truth particle.
     * @param name The name of the particle. 
     * @return The number of matches of a truth particle. 
     */
    int NMatchedTruth(const char * name);

    /**
     * Return the list of matches of a truth particle. The length of
     * the vector equals the number of reco objects of that type. The
     * entries are the matching status (-1: unknown, 0: no match; 1:
     * match).
     * @param index  The index of the particle.
     * @param ptype The type of the particle. 
     * @return The list of matches of a truth particle. 
     */
    std::vector<int> ListMatchedTruth(int index,  KLFitter::Particles::ParticleType ptype);

    /**
     * Return the list of matches of a truth particle. The length of
     * the vector equals the number of reco objects of that type. The
     * entries are the matching status (-1: unknown, 0: no match; 1:
     * match).
     * @param name The name of the particle. 
     * @return The list of matches of a truth particle. 
     */
    std::vector<int> ListMatchedTruth(const char * name);

    /* @} */
    /** \name Member functions (Set)  */
    /* @{ */

    /**
     * Set delta R matching criterion for partons.
     * @param deltaR The maximum deltaR to match two particles.
     */
    void SetDeltaRPartons(double deltaR);

    /**
     * Set delta R matching criterion for electrons.
     * @param deltaR The maximum deltaR to match two particles.
     */
    void SetDeltaRElectrons(double deltaR);

    /**
     * Set delta R matching criterion for muons.
     * @param deltaR The maximum deltaR to match two particles.
     */
    void SetDeltaRMuons(double deltaR);

    /**
     * Set delta R matching criterion for photons.
     * @param deltaR The maximum deltaR to match two particles.
     */
    void SetDeltaRPhotons(double deltaR);

    /**
     * Set the matching status.
     * @param indextruth The truth particle. 
     * @param indexreco The reco particle. 
     * @param ptype The type of the particle. 
     * @param stat The status. 
     */ 
    int SetMatchingStatus(int indextruth, int indexreco, KLFitter::Particles::ParticleType ptype, int stat); 

    /* @} */
    /** \name Member functions (misc)  */
    /* @{ */

    /**
     * Initialize vectors of a certain type.
     * @param ptype The type of the particle. 
     * @return An error code.
     */
    int Initialize(KLFitter::Particles::ParticleType ptype);

    /**
     * Match all truth particles of a certain type to reco particles. 
     * @param ptype The particle type.
     * @return An error code. 
     */
    int MatchTruthAll(KLFitter::Particles::ParticleType ptype);

    /**
     * Match a truth particle of a certain type to all reco particles. 
     * @param index  The index of the particle.
     * @param ptype The type of the particle. 
     * @return Matching status flag (1: at least one match, 0: no matches).
     */
    int MatchTruth(int index, KLFitter::Particles::ParticleType ptype);

    /**
     * Check if two objects of a certain type are matched. 
     * @param index1 Index of measured particle.
     * @param index2 Index for truth particle. 
     * @param ptype The particle type. 
     * @return Matching status flag (1: matched, 0: not matched). 
     */
    virtual int MatchVectors(int index1, int index2, KLFitter::Particles::ParticleType ptype);

    /**
     * Perform matching of two Lorentz vectors based on a deltaR
     * criterion.
     * @param vect1 A pointer to a particle. 
     * @param vect2 A pointer to a particle. 
     * @param dR The delteR criterion.
     * @return Matching status flag (1: matched, 0: not matched). 
     */
    int MatchDeltaR(TLorentzVector * vect1, TLorentzVector * vect2, double dR);

    /* @} */

  protected:

    /**
     * A pointer to the original particles. 
     */ 
    KLFitter::Particles ** fParticles; 

    /**
     * A pointer to the original truth particles. 
     */ 
    KLFitter::Particles ** fParticlesTruth; 

    /**
     * DeltaR criterion for partons.
     */
    double fDeltaRPartons;

    /**
     * DeltaR criterion for electrons.
     */
    double fDeltaRElectrons;

    /**
     * DeltaR criterion for muons.
     */
    double fDeltaRMuons;

    /**
     * DeltaR criterion for photons.
     */
    double fDeltaRPhotons;

    /**
     * A pointer to the vector of vectors matching truth partons
     */
    std::vector< std::vector<int> > * fMatchedPartons;

    /**
     * A pointer to the vector of vectors matching truth electrons
     */
    std::vector< std::vector<int> > * fMatchedElectrons;

    /**
     * A pointer to the vector of vectors matching truth muons
     */
    std::vector< std::vector<int> > * fMatchedMuons; 

    /**
     * A pointer to the vector of vectors matching truth photons
     */
    std::vector< std::vector<int> > * fMatchedPhotons; 
  };

}// namespace KLFitter 

#endif 
