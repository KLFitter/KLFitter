/*!
 * \class KLFitter::Particles
 * \brief A class describing particles. 
 * \author Kevin Kr&ouml;ninger
 * \version 1.3
 * \date 03.12.2009
 *
 * This class contains sets of TLorentzVectors for quarks, leptons,
 * etc.
 */

// --------------------------------------------------------- 

#ifndef PARTICLES
#define PARTICLES

// --------------------------------------------------------- 

#include <TLorentzVector.h>
#include <vector> 
#include <string>

// --------------------------------------------------------- 

class TLorentzVector;

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{

  class Particles
  {
                
  public: 
                
    /** \name Enumerators */
    /* @{ */

    /**
     * An enumerator for the particle types 
     */
    enum ParticleType { kParton, kElectron, kMuon, kTau, kNeutrino, kBoson, kPhoton };

    /* @} */
    /** \name Constructors and destructors */ 
    /* @{ */ 
                
    /** 
     * The default constructor. 
     */ 
    Particles(); 
                
    /**
     * The default destructor.
     */
    virtual ~Particles(); 

    /* @} */
    /** \name Member functions (Get)  */
    /* @{ */

    /**
     * Return the number of partons. 
     * @return The number of partons. 
     */ 
    int NPartons()
    { return int(fPartons -> size()); }; 

    /**
     * Return the number of electrons. 
     * @return The number of electrons. 
     */ 
    int NElectrons()
    { return int (fElectrons -> size()); }; 

    /**
     * Return the number of muons. 
     * @return The number of muons. 
     */ 
    int NMuons()
    { return int (fMuons -> size()); }; 

    /**
     * Return the number of taus. 
     * @return The number of taus. 
     */ 
    int NTaus()
    { return int (fTaus -> size()); }; 

    /**
     * Return the number of neutrinos. 
     * @return The number of neutrinos. 
     */ 
    int NNeutrinos()
    { return int (fNeutrinos -> size()); }; 

    /**
     * Return the number of bosons. 
     * @return The number of bosons. 
     */ 
    int NBosons()
    { return int (fBosons -> size()); }; 

    /**
     * Return the number of photons. 
     * @return The number of photons. 
     */ 
    int NPhotons()
    { return int (fPhotons -> size()); }; 

    /**
     * Return the particle with a certain name
     * @param name The name of the particle.  
     * @return A pointer to the TLorentzVector of the particle. 
     */ 
    TLorentzVector * Particle(std::string name); 

    /**
     * Return a particle with some index and type.
     * @param index The index of the particle.
     * @param ptype The type of the particle. 
     * @return A pointer to the TLorentzVector of the particle. 
     */ 
    TLorentzVector * Particle(int index, KLFitter::Particles::ParticleType ptype); 

    /**
     * Find a particle by name. 
     * @param name The name of the particle. 
     * @param particle A reference to the pointer to the 4-vector. 
     * @param index The reference of the index. 
     * @param ptype The reference of the particle type. 
     * @return A flag (1: found, 0: not found). 
     */ 
    int FindParticle(std::string name, TLorentzVector * &particle, int &index, KLFitter::Particles::ParticleType &ptype);

    /**
     * Return the parton at some index. 
     * @param index The parton index 
     * @return A pointer to the TLorentzVector of the parton. 
     */ 
    TLorentzVector * Parton(int index); 

    /**
     * Return the electron at some index. 
     * @param index The electron index 
     * @return A pointer to the TLorentzVector of the electron. 
     */ 
    TLorentzVector * Electron(int index); 

    /**
     * Return the muon at some index. 
     * @param index The muon index 
     * @return A pointer to the TLorentzVector of the muon. 
     */ 
    TLorentzVector * Muon(int index); 

    /**
     * Return the tau at some index. 
     * @param index The tau index 
     * @return A pointer to the TLorentzVector of the tau. 
     */ 
    TLorentzVector * Tau(int index); 

    /**
     * Return the boson at some index. 
     * @param index The boson index 
     * @return A pointer to the TLorentzVector of the boson. 
     */ 
    TLorentzVector * Boson(int index); 

    /**
     * Return the neutrino at some index. 
     * @param index The neutrino index 
     * @return A pointer to the TLorentzVector of the neutrino. 
     */ 
    TLorentzVector * Neutrino(int index); 

    /**
     * Return the photon at some index. 
     * @param index The photon index 
     * @return A pointer to the TLorentzVector of the photon. 
     */ 
    TLorentzVector * Photon(int index); 

    /**
     * Return the number of particles. 
     * @return The number of particles. 
     */ 
    int NParticles()
    { return int(fPartons -> size() + fElectrons -> size() + fMuons -> size() + fTaus -> size() + fNeutrinos -> size() + fBosons -> size() + fPhotons -> size()); }; 
                
    /**
     * Return the number of particles of a certain type. 
     * @param ptype The particle type. 
     * @return The number of particles. 
     */ 
    int NParticles(KLFitter::Particles::ParticleType ptype); 

    /** 
     * Return the name of a particle. 
     * @param index The index of the particle.
     * @param ptype The type of the particle. 
     * @return The name of the particle. 
     */ 
    std::string NameParticle(int index, KLFitter::Particles::ParticleType ptype); 

    /** 
     * Return the name of a electron. 
     * @param index The index of the electron. 
     * @return The name of the electron. 
     */ 
    std::string NameElectron(int index); 

    /** 
     * Return the name of a muon. 
     * @param index The index of the muon. 
     * @return The name of the muon. 
     */ 
    std::string NameMuon(int index); 

    /** 
     * Return the name of a tau. 
     * @param index The index of the tau. 
     * @return The name of the tau. 
     */ 
    std::string NameTau(int index); 

    /** 
     * Return the name of a boson. 
     * @param index The index of the boson. 
     * @return The name of the boson. 
     */ 
    std::string NameBoson(int index); 

    /** 
     * Return the name of a neutrino. 
     * @param index The index of the neutrino. 
     * @return The name of the neutrino. 
     */ 
    std::string NameNeutrino(int index); 

    /** 
     * Return the name of a photon. 
     * @param index The index of the photon. 
     * @return The name of the photon. 
     */ 
    std::string NamePhoton(int index); 

    /**
     * Return the corresponding measured particle. 
     * @param index The index of the model particle. 
     * @param index The index of the measured particle. 
     */ 
    int JetIndex(int index);

    int ElectronIndex(int index);

    int MuonIndex(int index);

    int PhotonIndex(int index);

    /**
     * Return the b-tagging probability of a parton. 
     * measured particles: measured b-tagging probability (or discriminant)\n
     * model particles: true flavor (0: light, 1:b, 2:other).
     * @param index The parton index
     * @return The b-tagging probability 
     */ 
    double BTaggingProbability(int index); 

    /**
     * Return the experimental flavor tag of a parton. 
     * 0: tagged as light jet\n
     * 1: tagged as b-jet
     * @param index The parton index
     * @return The flavor tag
     */ 
    double FlavorTag(int index); 

    /**
     * Return the number of b-tags.
     */
    int NBTags(); 

    /* @} */
    /** \name Member functions (Set)  */
    /* @{ */

    /**
     * Set the flavor tag
     * @param index The parton index. 
     * @param prob The flavor tag.
     * @return An error flag.
     */ 
    int SetFlavorTag(int index, double tag); 
                
    /**
     * Set the b-tagging probability.
     * @param index The parton index. 
     * @param prob The b-tagging probability. 
     * @return An error flag.
     */ 
    int SetBTaggingProbability(int index, double prob);

    /* @} */
    /** \name Member functions (misc)  */
    /* @{ */

    /**
     * Add a particle to a list of particles. 
     * @param particle A pointer to the particle. 
     * @param ptype The type of particle.
     * @param name The name of the particle.
     * @param btagprob The b-tagging probability. 
     * @param flavortag The flavor-tag
     * @param measuredindex The index of the associated measured particle. 
     * @return An error code.
     */ 
    int AddParticle(TLorentzVector * particle, KLFitter::Particles::ParticleType ptype, std::string name = "", double btagprob = 0, double flavortag = 0, int measuredindex = -1); 

    /**
     * Removes a particle from a list of particles. 
     * @param index The index of the particle.
     * @param ptype The type of the particle.
     * @return An error code.
     */ 
    int RemoveParticle(int index, KLFitter::Particles::ParticleType ptype); 

    /**
     * Removes a particle from a list of particles. 
     * @param name The name of the particle.
     * @return An error code.
     */ 
    int RemoveParticle(std::string name); 

    /**
     * Return the particle container of a type of particles
     * @param ptype The type of the particle. 
     * @return The particle container. 
     */ 
    std::vector <TLorentzVector *> * ParticleContainer(KLFitter::Particles::ParticleType ptype); 

    /**
     * Return the particle name container of a type of particles
     * @param ptype The type of the particle. 
     * @return The particle name container. 
     */ 
    std::vector <std::string> * ParticleNameContainer(KLFitter::Particles::ParticleType ptype);

    /**
     * Checks if the index is within range.
     * @param container The particle container. 
     * @param index The index of particle.
     * @return An error flag. 
     */ 
    int CheckIndex(std::vector <TLorentzVector *> * container, int index);
    
    /**
     * Compare function for std::sort with pT as comparison criteria.
     * @param p1 First particle. 
     * @param p2 Second particle.
     * @return bool. 
     */ 
    static bool PtOrder(const TLorentzVector *p1, const TLorentzVector * p2){return p1->Pt() > p2->Pt();}

    /* @} */

  private: 

  protected: 

  private: 

    /**
     * A set of quarks and gluons. 
     */ 
    std::vector <TLorentzVector *> * fPartons; 

    /**
     * A set of electrons. 
     */ 
    std::vector <TLorentzVector *> * fElectrons; 

    /**
     * A set of muons. 
     */ 
    std::vector <TLorentzVector *> * fMuons; 

    /**
     * A set of taus. 
     */ 
    std::vector <TLorentzVector *> * fTaus; 

    /**
     * A set of neutrinos. 
     */ 
    std::vector <TLorentzVector *> * fNeutrinos; 

    /**
     * A set of bosons. 
     */ 
    std::vector <TLorentzVector *> * fBosons; 

    /**
     * A set of photons. 
     */ 
    std::vector <TLorentzVector *> * fPhotons; 

    /**
     * The name of the partons.
     */ 
    std::vector <std::string> * fNamePartons; 

    /**
     * The name of the electrons.
     */ 
    std::vector <std::string> * fNameElectrons; 

    /**
     * The name of the muons.
     */ 
    std::vector <std::string> * fNameMuons; 

    /**
     * The name of the taus.
     */ 
    std::vector <std::string> * fNameTaus; 

    /**
     * The name of the neutrinos.
     */ 
    std::vector <std::string> * fNameNeutrinos; 

    /**
     * The name of the bosons.
     */ 
    std::vector <std::string> * fNameBosons; 

    /**
     * The name of the photons.
     */ 
    std::vector <std::string> * fNamePhotons; 

    /**
     * The index of the corresponding measured parton.
     */ 
    std::vector <int> * fJetIndex; 

    /**
     * The index of the corresponding measured electron.
     */ 
    std::vector <int> * fElectronIndex; 

    /**
     * The index of the corresponding measured muon.
     */ 
    std::vector <int> * fMuonIndex; 

    /**
     * The index of the corresponding measured photon.
     */ 
    std::vector <int> * fPhotonIndex; 

    /**
     * Vector containing the b-tagging probability.
     * measured particles: measured b-tagging probability (or discriminant)\n
     * model particles: true flavor (0: light, 1:b, 2:other).
     */ 
    std::vector <double> * fBTaggingProbability; 

    /**
     * Vector containing the experimental flavor tag. 
     * 0: tagged as light jet\n
     * 1: tagged as b-jet
     */ 
    std::vector <double> * fFlavorTag; 

  }; 

} // namespace KLFitter 

inline 
KLFitter::Particles::ParticleType &operator++(KLFitter::Particles::ParticleType &ptype)
{
  return ptype = KLFitter::Particles::ParticleType(ptype + 1); 
}

// --------------------------------------------------------- 

#endif 

