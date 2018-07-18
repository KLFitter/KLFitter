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

#ifndef KLFITTER_PARTICLES_H_
#define KLFITTER_PARTICLES_H_

#include <string>
#include <vector>

#include "TLorentzVector.h"

// ---------------------------------------------------------

/**
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter {
/**
  * \class KLFitter::Particles
  * \brief A class describing particles.
  *
  * This class contains sets of TLorentzVectors for quarks, leptons,
  * etc.
  */
class Particles final {
 public:
  /// Enumerator for particle types.
  enum ParticleType { kParton,    ///< Particle type parton
                      kElectron,  ///< Particle type electron
                      kMuon,      ///< Particle type muon
                      kTau,       ///< Particle type tau
                      kNeutrino,  ///< Particle type neutrino
                      kBoson,     ///< Particle type boson
                      kPhoton     ///< Particle type photon
  };

  /// An enumerator for the true jet flavor.
  enum TrueFlavorType { kLight,     ///< Light quark
                        kB,         ///< B-quark
                        kLightUp,   ///< Up-type light quark
                        kLightDown, ///< Down-type light quark
                        kNone       ///< Not specified
  };

  /** \name Constructors and destructors */
  /* @{ */

  /// The default constructor.
  Particles();

  /// The copy constructor.
  explicit Particles(const Particles& o);

  /// The (defaulted) destructor.
  ~Particles();

  /// The assignment operator.
  Particles& operator=(const Particles& o);

  /* @} */
  /** \name Member functions (Get)  */
  /* @{ */

  /**
   * Return the number of partons.
   * @return The number of partons.
   */
  int NPartons() const { return static_cast<int>(m_partons.size()); }

  /**
   * Return the number of electrons.
   * @return The number of electrons.
   */
  int NElectrons() const { return int (m_electrons.size()); }

  /**
   * Return the number of muons.
   * @return The number of muons.
   */
  int NMuons() const { return int (m_muons.size()); }

  /**
   * Return the number of taus.
   * @return The number of taus.
   */
  int NTaus() const { return int (m_taus.size()); }

  /**
   * Return the number of neutrinos.
   * @return The number of neutrinos.
   */
  int NNeutrinos() const { return int (m_neutrinos.size()); }

  /**
   * Return the number of bosons.
   * @return The number of bosons.
   */
  int NBosons() const { return int (m_bosons.size()); }

  /**
   * Return the number of photons.
   * @return The number of photons.
   */
  int NPhotons() const { return int (m_photons.size()); }

  /**
   * Return the particle with a certain name
   * @param name The name of the particle.
   * @return A pointer to the TLorentzVector of the particle.
   */
  TLorentzVector* Particle(const std::string& name);

  /**
   * Return a particle with some index and type.
   * @param index The index of the particle.
   * @param ptype The type of the particle.
   * @return A pointer to the TLorentzVector of the particle.
   */
  TLorentzVector* Particle(int index, KLFitter::Particles::ParticleType ptype);

  /**
   * Find a particle by name.
   * @param name The name of the particle.
   * @param particle A reference to the pointer to the 4-vector.
   * @param index The pointer to the index.
   * @param ptype The pointer to the particle type.
   * @return A flag (1: found, 0: not found).
   */
  int FindParticle(const std::string& name, TLorentzVector* &particle, int* index, KLFitter::Particles::ParticleType* ptype);

  /**
   * Return the parton at some index.
   * @param index The parton index
   * @return A pointer to the TLorentzVector of the parton.
   */
  TLorentzVector* Parton(int index);

  /**
   * Return the electron at some index.
   * @param index The electron index
   * @return A pointer to the TLorentzVector of the electron.
   */
  TLorentzVector* Electron(int index);

  /**
   * Return the muon at some index.
   * @param index The muon index
   * @return A pointer to the TLorentzVector of the muon.
   */
  TLorentzVector* Muon(int index);

  /**
   * Return the tau at some index.
   * @param index The tau index
   * @return A pointer to the TLorentzVector of the tau.
   */
  TLorentzVector* Tau(int index);

  /**
   * Return the boson at some index.
   * @param index The boson index
   * @return A pointer to the TLorentzVector of the boson.
   */
  TLorentzVector* Boson(int index);

  /**
   * Return the neutrino at some index.
   * @param index The neutrino index
   * @return A pointer to the TLorentzVector of the neutrino.
   */
  TLorentzVector* Neutrino(int index);

  /**
   * Return the photon at some index.
   * @param index The photon index
   * @return A pointer to the TLorentzVector of the photon.
   */
  TLorentzVector* Photon(int index);

  /**
   * Return the number of particles.
   * @return The number of particles.
   */
  int NParticles() const { return static_cast<int>(m_partons.size() + m_electrons.size() + m_muons.size() + m_taus.size() + m_neutrinos.size() + m_bosons.size() + m_photons.size()); }

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
  std::string NameElectron(int index) const;

  /**
   * Return the name of a muon.
   * @param index The index of the muon.
   * @return The name of the muon.
   */
  std::string NameMuon(int index) const;

  /**
   * Return the name of a tau.
   * @param index The index of the tau.
   * @return The name of the tau.
   */
  std::string NameTau(int index) const;

  /**
   * Return the name of a boson.
   * @param index The index of the boson.
   * @return The name of the boson.
   */
  std::string NameBoson(int index) const;

  /**
   * Return the name of a neutrino.
   * @param index The index of the neutrino.
   * @return The name of the neutrino.
   */
  std::string NameNeutrino(int index) const;

  /**
   * Return the name of a photon.
   * @param index The index of the photon.
   * @return The name of the photon.
   */
  std::string NamePhoton(int index) const;

  /**
   * Return the corresponding measured particle.
   * @param index The index of the model particle.
   * @param index The index of the measured particle.
   */
  int JetIndex(int index) const;

  int ElectronIndex(int index) const;

  int MuonIndex(int index) const;

  int PhotonIndex(int index) const;

  /**
   * Return the true flavor of a parton.
   * @param index The parton index
   * @return The parton true flavor.
   */
  TrueFlavorType TrueFlavor(int index) const { return m_true_flavors[index]; }

  /**
   * Return has the jet been b-tagged?
   * @param index The parton index
   * @return The parton b-tagging boolean.
   */
  bool IsBTagged(int index) const { return m_jet_btagged_bools[index]; }

  /**
   * Return the jet b-tagging efficiency.
   * @param index The parton index
   * @return The jet b-tagging efficiency.
   */
  double BTaggingEfficiency(int index) const { return m_btag_efficiencies[index]; }

  /**
   * Return the jet b-tagging rejection.
   * @param index The parton index
   * @return The jet b-tagging rejection.
   */
  double BTaggingRejection(int index) const { return m_btag_rejections[index]; }

  /**
   * Return the jet b-tagging weight.
   * @param index The parton index
   * @return The jet b-tagging weight.
   */
  double BTagWeight(int index) const { return m_btag_weights[index]; }

  /**
   * Return the bool of a set tagging weight.
   * @param index The parton index
   * @return The bool of a set tagging weight
   */
  bool BTagWeightSet(int index) const { return m_btag_weights_set[index]; }

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

  /// Return the number of b-tags.
  int NBTags() const;

  /* @} */
  /** \name Member functions (Set)  */
  /* @{ */

  /**
   * Set has the jet been b-tagged?
   * @param index The parton index
   * @param isBTagged The parton b-tagging boolean.
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
   * @param btagweight The b-tagging weight.
   * @return An error flag.
   */
  int SetBTagWeight(int index, double btagweight);

  /**
   * Set bool for set jet b-tagging weight.
   * @param index The parton index
   * @param btagweightset The b-tagging probability.
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
   * @param LepCharge The Charge of the particle.
   * @param ptype The type of particle.
   * @param name The name of the particle.
   * @param measuredindex The index of the associated measured particle.
   * @return An error code.
   */
  int AddParticle(const TLorentzVector& particle, double DetEta, float LepCharge, KLFitter::Particles::ParticleType ptype, std::string name = "", int measuredindex = -1);

  /**
   * DEPRECATED FUNCTION. This is an overloaded implementation of
   * the previous function using TLorentzVector pointers instead
   * of const ref. The usage of this function is deprecated and
   * it will be removed in the next major release. Please switch
   * to the above implementation.
   */
  int AddParticle(const TLorentzVector* const particle, double DetEta, float LepCharge, KLFitter::Particles::ParticleType ptype, std::string name = "", int measuredindex = -1);

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
  int AddParticle(const TLorentzVector& particle, double DetEta, KLFitter::Particles::ParticleType ptype, std::string name = "", int measuredindex = -1, bool isBtagged = false, double bTagEff = -1., double bTagRej = -1., TrueFlavorType trueflav = kNone, double btagweight = 999);

  /**
   * DEPRECATED FUNCTION. This is an overloaded implementation of
   * the previous function using TLorentzVector pointers instead
   * of const ref. The usage of this function is deprecated and
   * it will be removed in the next major release. Please switch
   * to the above implementation.
   */
  int AddParticle(const TLorentzVector* const particle, double DetEta, KLFitter::Particles::ParticleType ptype, std::string name = "", int measuredindex = -1, bool isBtagged = false, double bTagEff = -1., double bTagRej = -1., TrueFlavorType trueflav = kNone, double btagweight = 999);

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
  int AddParticle(const TLorentzVector& particle, KLFitter::Particles::ParticleType ptype, std::string name = "", int measuredindex = -1, bool isBtagged = false, double bTagEff = -1., double bTagRej = -1., TrueFlavorType trueflav = kNone, double btagweight = 999);

  /**
   * DEPRECATED FUNCTION. This is an overloaded implementation of
   * the previous function using TLorentzVector pointers instead
   * of const ref. The usage of this function is deprecated and
   * it will be removed in the next major release. Please switch
   * to the above implementation.
   */
  int AddParticle(const TLorentzVector* const particle, KLFitter::Particles::ParticleType ptype, std::string name = "", int measuredindex = -1, bool isBtagged = false, double bTagEff = -1., double bTagRej = -1., TrueFlavorType trueflav = kNone, double btagweight = 999);

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
  int AddParticle(const TLorentzVector& particle, KLFitter::Particles::ParticleType ptype, std::string name, int measuredindex, TrueFlavorType trueflav, double btagweight = 999);

  /**
   * DEPRECATED FUNCTION. This is an overloaded implementation of
   * the previous function using TLorentzVector pointers instead
   * of const ref. The usage of this function is deprecated and
   * it will be removed in the next major release. Please switch
   * to the above implementation.
   */
  int AddParticle(const TLorentzVector* const particle, KLFitter::Particles::ParticleType ptype, std::string name, int measuredindex, TrueFlavorType trueflav, double btagweight = 999);

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
  int RemoveParticle(const std::string& name);

  /**
   * Return the particle container of a type of particles
   * @param ptype The type of the particle.
   * @return The particle container.
   */
  std::vector<TLorentzVector>* ParticleContainer(KLFitter::Particles::ParticleType ptype);

  /**
   * Return the particle name container of a type of particles
   * @param ptype The type of the particle.
   * @return The particle name container.
   */
  std::vector<std::string>* ParticleNameContainer(KLFitter::Particles::ParticleType ptype);

  /**
   * Checks if the index is within range.
   * @param container The particle container.
   * @param index The index of particle.
   * @return An error flag.
   */
  int CheckIndex(const std::vector<TLorentzVector>& container, int index) const;

  /* @} */

 private:
  /**
   * A set of quarks and gluons.
   */
  std::vector<TLorentzVector> m_partons;

  /**
   * A set of electrons.
   */
  std::vector<TLorentzVector> m_electrons;

  /**
   * A set of muons.
   */
  std::vector<TLorentzVector> m_muons;

  /**
   * A set of taus.
   */
  std::vector<TLorentzVector> m_taus;

  /**
   * A set of neutrinos.
   */
  std::vector<TLorentzVector> m_neutrinos;

  /**
   * A set of bosons.
   */
  std::vector<TLorentzVector> m_bosons;

  /**
   * A set of photons.
   */
  std::vector<TLorentzVector> m_photons;

  /**
   * The name of the partons.
   */
  std::vector<std::string> m_parton_names;

  /**
   * The name of the electrons.
   */
  std::vector<std::string> m_electron_names;

  /**
   * The name of the muons.
   */
  std::vector<std::string> m_muon_names;

  /**
   * The name of the taus.
   */
  std::vector<std::string> m_tau_names;

  /**
   * The name of the neutrinos.
   */
  std::vector<std::string> m_neutrino_names;

  /**
   * The name of the bosons.
   */
  std::vector<std::string> m_boson_names;

  /**
   * The name of the photons.
   */
  std::vector<std::string> m_photon_names;

  /**
   * The index of the corresponding measured parton.
   */
  std::vector<int> m_jet_indices;

  /**
   * The index of the corresponding measured electron.
   */
  std::vector<int> m_electron_indices;

  /**
   * The index of the corresponding measured muon.
   */
  std::vector<int> m_muon_indices;

  /**
   * The index of the corresponding measured photon.
   */
  std::vector<int> m_photon_indices;

  /**
   * Vector containing the true flavor.
   */
  std::vector<TrueFlavorType> m_true_flavors;

  /**
   * Vector containing a boolean for the b-tagging.
   */
  std::vector<bool> m_jet_btagged_bools;

  /**
   * Vector containing the b-tagging efficiencies for the jets.
   */
  std::vector<double> m_btag_efficiencies;

  /**
   * Vector containing the b-tagging rejection for the jets.
   */
  std::vector<double> m_btag_rejections;

  /**
   * Vector containing the b-tagging weights for the jets.
   */
  std::vector<double> m_btag_weights;

  /**
   * Vector containing the bool if b-tagging weights for the jets were set.
   */
  std::vector<bool> m_btag_weights_set;

  /**
   * Vector containing the detector eta of electrons.
   */
  std::vector<double> m_electron_det_etas;

  /**
   * Vector containing the detector eta of muons.
   */
  std::vector<double> m_muon_det_etas;

  /**
   * Vector containing the detector eta of jets.
   */
  std::vector<double> m_jet_det_etas;

  /**
   * Vector containing the detector eta of photons.
   */
  std::vector<double> m_photon_det_etas;

  /**
   * Vector containing the charge of electrons.
   */
  std::vector<float> m_electron_charges;

  /**
   * Vector containing the charge of muons.
   */
  std::vector<float> m_muon_charges;
};
}  // namespace KLFitter

inline KLFitter::Particles::ParticleType &operator++(KLFitter::Particles::ParticleType &ptype) {
  return ptype = KLFitter::Particles::ParticleType(ptype + 1);
}

#endif  // KLFITTER_PARTICLES_H_
