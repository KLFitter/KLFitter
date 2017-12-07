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

    /**
     * An enumerator for the true jet flavor.
     */
    enum TrueFlavorType { kLight, kB, kLightUp, kLightDown, kNone };

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
     * Return the true flavor of a parton.
     * @param index The parton index
     * @return The parton true flavor.
     */
    TrueFlavorType TrueFlavor(int index) { return (*fTrueFlavor)[index]; }

    /**
     * Return has the jet been b-tagged?
     * @param index The parton index
     * @return The parton b-tagging boolean.
     */
    bool IsBTagged(int index) { return (*fIsBTagged)[index]; }

    /**
     * Return the jet b-tagging efficiency.
     * @param index The parton index
     * @return The jet b-tagging efficiency.
     */
    double BTaggingEfficiency(int index) { return (*fBTaggingEfficiency)[index]; }

    /**
     * Return the jet b-tagging rejection.
     * @param index The parton index
     * @return The jet b-tagging rejection.
     */
    double BTaggingRejection(int index) { return (*fBTaggingRejection)[index]; }

    /**
     * Return the jet b-tagging weight.
     * @param index The parton index
     * @return The jet b-tagging weight.
     */
    double BTagWeight(int index) { return (*fBTagWeight)[index]; }

    /**
     * Return the bool of a set tagging weight.
     * @param index The parton index
     * @return The bool of a set tagging weight
     */
    bool BTagWeightSet(int index) { return (*fBTagWeightSet)[index]; }

    /**
     * Return the detector eta of a particle with some index and type.
     * @param index The index of the particle
     * @param ptype The particle type.
     * @return The detector eta of the particle
     */
    double DetEta(int index, KLFitter::Particles::ParticleType ptype);

     /**
     * Return the charge of the lepton with some index and type.
     * @param index The index of the particle
     * @param ptype The particle type.
     * @return The charge of the lepton.
     */
    float LeptonCharge(int index, KLFitter::Particles::ParticleType ptype);

    /**
     * Return the number of b-tags.
     */
    int NBTags();

    /* @} */
    /** \name Member functions (Set)  */
    /* @{ */

    /**
     * Set has the jet been b-tagged?
     * @param index The parton index
     * @param isBtagged The parton b-tagging boolean.
     * @return An error flag.
     */
    int SetIsBTagged(int index, bool isBTagged);

    /**
     * Set the jet b-tagging efficiency.
     * @param index The parton index
     * @param btagEff The b-tagging efficiency.
     * @return An error flag.
     */
    int SetBTaggingEfficiency(int index, double btagEff);

    /**
     * Set the jet b-tagging rejection.
     * @param index The parton index
     * @param btagRej The b-tagging probability.
     * @return The jet b-tagging rejection.
     */
    int SetBTaggingRejection(int index, double btagRej);

    /**
     * Set the jet b-tagging weight.
     * @param index The parton index
     * @param btagRej The b-tagging weight.
     * @return An error flag.
     */
    int SetBTagWeight(int index, double btagweight);

    /**
     * Set bool for set jet b-tagging weight.
     * @param index The parton index
     * @param btagRej The b-tagging probability.
     * @return An error flag.
     */
    int SetBTagWeightSet(int index, bool btagweightset);

    /* @} */
    /** \name Member functions (misc)  */
    /* @{ */

    /**
     * Add a particle to a list of particles.
     * @param particle A pointer to the particle.
     * @param DetEta The Detector Eta of the particle.
     * @param LepEta The Charge of the particle.
     * @param ptype The type of particle.
     * @param name The name of the particle.
     * @param measuredindex The index of the associated measured particle.
     * @param isBtagged Has the particle been b-tagged?
     * @param bTagEff B-tagging efficiency of the particle.
     * @param bTagRej B-tagging rejection of the particle.
     * @param trueflav The true flavor (only for model particles).
     * @param btagweight The b tagger weight).
     * @return An error code.
     */
    int AddParticle(TLorentzVector * particle, double DetEta, float LepCharge, KLFitter::Particles::ParticleType ptype, std::string name = "", int measuredindex = -1);

    /**
     * Add a particle to a list of particles.
     * @param particle A pointer to the particle.
     * @param DetEta The Detector Eta of the particle.
     * @param ptype The type of particle.
     * @param name The name of the particle.
     * @param measuredindex The index of the associated measured particle.
     * @param isBtagged Has the particle been b-tagged?
     * @param bTagEff B-tagging efficiency of the particle.
     * @param bTagRej B-tagging rejection of the particle.
     * @param trueflav The true flavor (only for model particles).
     * @param btagweight The b tagger weight).
     * @return An error code.
     */
    int AddParticle(TLorentzVector * particle, double DetEta, KLFitter::Particles::ParticleType ptype, std::string name = "", int measuredindex = -1, bool isBtagged = false, double bTagEff = -1., double bTagRej = -1., TrueFlavorType trueflav = kNone, double btagweight = 999);


    /**
     * Add a particle to a list of particles.
     * @param particle A pointer to the particle.
     * @param ptype The type of particle.
     * @param name The name of the particle.
     * @param isBtagged Has the particle been b-tagged?
     * @param bTagEff B-tagging efficiency of the particle.
     * @param bTagRej B-tagging rejection of the particle.
     * @param measuredindex The index of the associated measured particle.
     * @param trueflav The true flavor (only for model particles).
     * @param btagweight The b tagger weight).
     * @return An error code.
     */
    int AddParticle(TLorentzVector * particle, KLFitter::Particles::ParticleType ptype, std::string name = "", int measuredindex = -1, bool isBtagged = false, double bTagEff = -1., double bTagRej = -1., TrueFlavorType trueflav = kNone, double btagweight = 999);

    /**
     * Add a particle to a list of particles (especially for model particles).
     * @param particle A pointer to the particle.
     * @param ptype The type of particle.
     * @param name The name of the particle.
     * @param measuredindex The index of the associated measured particle.
     * @param trueflav The true flavor (only for model particles).
     * @param btagweight The b tagger weight).
     * @return An error code.
     */
    int AddParticle(TLorentzVector * particle, KLFitter::Particles::ParticleType ptype, std::string name, int measuredindex, TrueFlavorType trueflav, double btagweight = 999);

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
     * Vector containing the true flavor.
     */
    std::vector<TrueFlavorType> * fTrueFlavor;

    /**
     * Vector containing a boolean for the b-tagging.
     */
    std::vector<bool> * fIsBTagged;

    /**
     * Vector containing the b-tagging efficiencies for the jets.
     */
    std::vector<double> * fBTaggingEfficiency;

    /**
     * Vector containing the b-tagging rejection for the jets.
     */
    std::vector<double> * fBTaggingRejection;

    /**
     * Vector containing the b-tagging weights for the jets.
     */
    std::vector<double> * fBTagWeight;

    /**
     * Vector containing the bool if b-tagging weights for the jets were set.
     */
    std::vector<bool> * fBTagWeightSet;

    /**
     * Vector containing the detector eta of electrons.
     */
    std::vector <double> * fElectronDetEta;
    /**
     * Vector containing the detector eta of muons.
     */
    std::vector <double> * fMuonDetEta;
    /**
     * Vector containing the detector eta of jets.
     */
    std::vector <double> * fJetDetEta;
    /**
     * Vector containing the detector eta of photons.
     */
    std::vector <double> * fPhotonDetEta;
   /**
     * Vector containing the charge of electrons.
     */
    std::vector <float> * fElectronCharge;
    /**
     * Vector containing the charge of muons.
     */
    std::vector <float> * fMuonCharge;

  };

} // namespace KLFitter

inline
KLFitter::Particles::ParticleType &operator++(KLFitter::Particles::ParticleType &ptype)
{
  return ptype = KLFitter::Particles::ParticleType(ptype + 1);
}

// ---------------------------------------------------------

#endif

