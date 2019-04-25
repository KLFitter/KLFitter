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

#include "KLFitter/PermutationHandler.h"

#include <algorithm>
#include <iostream>
#include <numeric>
#include <sstream>

namespace KLFitter {
// ---------------------------------------------------------
bool Permutation::next_permutation() {
  return std::next_permutation(partons.begin(), partons.end())
    || std::next_permutation(electrons.begin(), electrons.end())
    || std::next_permutation(muons.begin(), muons.end())
    || std::next_permutation(photons.begin(), photons.end())
    || std::next_permutation(tracks.begin(), tracks.end());
}

// ---------------------------------------------------------
std::string Permutation::print() {
  std::stringstream ss;
  ss << "Partons = ( ";
  for (const auto& parton : partons) {
    ss << parton << " ";
  }
  ss << ")";
  return ss.str();
}

// ---------------------------------------------------------
PermutationHandler::PermutationHandler(ParticleCollection ** p, ParticleCollection ** pp)
  : m_particles(p)
  , m_particles_permuted(pp)
  , m_permutation_index(-1) {
  // empty
}

// ---------------------------------------------------------
PermutationHandler::PermutationHandler(const PermutationHandler& o) = default;

// ---------------------------------------------------------
PermutationHandler::~PermutationHandler() = default;

// ---------------------------------------------------------
PermutationHandler& PermutationHandler::operator=(const PermutationHandler& obj) = default;

// ---------------------------------------------------------
bool PermutationHandler::IsVetoed(int index) {
  return m_permutation_list.at(index).vetoed;
}

// ---------------------------------------------------------
int PermutationHandler::SetPermutation(int index) {
  // check index
  if (index < 0 || index >= NPermutations()) {
    std::cout << "KLFitter::Permutations::SetPermutation(). Index out of range." << std::endl;
    return 0;
  }

  if (!m_particles_permuted) {
    std::cout << "KLFitter::Permutations::SetPermutation(). Pointer to permuted particles not available." << std::endl;
    return 0;
  }

  // set permutation
  (*m_particles_permuted) = &m_particles_table[index];

  // set permutation index
  m_permutation_index = index;

  // no error
  return 1;
}

// ---------------------------------------------------------
int PermutationHandler::CreatePermutations() {
  // reset existing particle and permuation tables
  Reset();

  // create new table of particles
  m_particles_table = std::vector<ParticleCollection>{};

  // check particles
  CheckParticles();

  // get number of objects per category
  size_t npartons   = (*m_particles)->partons.size();
  size_t nelectrons = (*m_particles)->electrons.size();
  size_t nmuons     = (*m_particles)->muons.size();
  size_t nphotons     = (*m_particles)->photons.size();
  size_t ntracks    = (*m_particles)->tracks.size();

  // Create initial permutation from which the permutations are calculated.
  Permutation initial_perm;
  initial_perm.partons.resize(npartons);
  std::iota(initial_perm.partons.begin(), initial_perm.partons.end(), 0);
  initial_perm.electrons.resize(nelectrons);
  std::iota(initial_perm.electrons.begin(), initial_perm.electrons.end(), 0);
  initial_perm.muons.resize(nmuons);
  std::iota(initial_perm.muons.begin(), initial_perm.muons.end(), 0);
  initial_perm.photons.resize(nphotons);
  std::iota(initial_perm.photons.begin(), initial_perm.photons.end(), 0);
  initial_perm.tracks.resize(ntracks);
  std::iota(initial_perm.tracks.begin(), initial_perm.tracks.end(), 0);

  // Store a list of all possible permutations.
  do {
    m_permutation_list.emplace_back(initial_perm);
  } while (initial_perm.next_permutation());

  for (const auto& perm : m_permutation_list) {
    ParticleCollection particles{};

    // loop over all partons
    for (const auto& parton : perm.partons) {
      particles.AddParticle((*m_particles)->partons.at(parton));
    }

    // loop over all electrons
    for (const auto& electron : perm.electrons) {
      particles.AddParticle((*m_particles)->electrons.at(electron));
    }

    // loop over all muons
    for (const auto& muon : perm.muons) {
      particles.AddParticle((*m_particles)->muons.at(muon));
    }

    // loop over all photons
    for (const auto& photon : perm.photons) {
      particles.AddParticle((*m_particles)->photons.at(photon));
    }

    // loop over all tracks
    for (const auto& track : perm.tracks) {
      particles.AddParticle((*m_particles)->tracks.at(track));
    }

    // add particles to table
    m_particles_table.emplace_back(particles);
  }

  // no error
  return 1;
}

// ---------------------------------------------------------
int PermutationHandler::Reset() {
  // Clear particle and permutation tables.
  m_particles_table.clear();
  m_permutation_list.clear();

  // no error
  return 1;
}

// ---------------------------------------------------------
int PermutationHandler::InvariantParticlePermutations(Particles::Type ptype, std::vector<int> indices) {
  // check if particles are defined
  if (!CheckParticles()) return 0;

  // no error
  int err = 1;

  if (indices.size() == 2) {
    // If we only need to compare two indices, things are easy:
    // we can go through the list of permutations, and veto all
    // permutations, where the position of the first index is
    // larger than that of the second index. This way we only
    // keep permutations with pos(index1) < pos(index2).
    auto index1 = indices.at(0);
    auto index2 = indices.at(1);

    for (auto& iperm : m_permutation_list) {
      if (ptype == Particles::Type::kParton && iperm.partons.at(index1) > iperm.partons.at(index2)) {
        iperm.vetoed = true;
      }
      if (ptype == Particles::Type::kElectron && iperm.electrons.at(index1) > iperm.electrons.at(index2)) {
        iperm.vetoed = true;
      }
      if (ptype == Particles::Type::kMuon && iperm.muons.at(index1) > iperm.muons.at(index2)) {
        iperm.vetoed = true;
      }
      if (ptype == Particles::Type::kPhoton && iperm.photons.at(index1) > iperm.photons.at(index2)) {
        iperm.vetoed = true;
      }
      if (ptype == Particles::Type::kTrack && iperm.tracks.at(index1) > iperm.tracks.at(index2)) {
        iperm.vetoed = true;
      }
    }
  } else {
    // In case we have more than two indices, take the last entry
    // in the index vector and compare it with all the others.
    // When compared, remove that last entry from the index
    // vector. Repeat until only two indices are left.
    while (indices.size() >= 2) {
      int index2 = indices.back();
      for (unsigned int i = 0, I = indices.size()-1; i < I; i++) {
        int index1 = indices.at(i);
        err *= InvariantParticlePermutations(ptype, std::vector<int>{index1, index2});
      }
      indices.erase(--indices.end());
    }
  }

  // return error code
  return err;
}

// ---------------------------------------------------------
int PermutationHandler::InvariantParticleGroupPermutations(Particles::Type ptype, std::vector<int> positions1,  std::vector<int> positions2) {
  // Check if particles are defined.
  if (!CheckParticles()) return 0;

  // Perform several sanity checks on the index vectors.
  {
    // Check if index vectors have the same size.
    if (positions1.size()!= positions2.size()) {
      std::cout << "KLFitter::PermutationHandler::InvariantParticleGroupPermutations(). ";
      std::cout << "Index Vectors need to have the same size." << std::endl;
      return 0;
    }

    // To check whether any position appears twice, append the
    // two position vectors to each other, sort the total vector
    // and use std::unique() to check for duplicates. If
    // duplicates are found, return with an error.
    auto sorted = positions1;
    sorted.insert(sorted.end(), positions2.begin(), positions2.end());
    std::sort(sorted.begin(), sorted.end());
    if (std::unique(sorted.begin(), sorted.end()) != sorted.end()) {
      std::cout << "KLFitter::PermutationHandler::InvariantParticleGroupPermutations(). ";
      std::cout << "Position vectors contain at least one duplicate entry" << std::endl;
      return 0;
    }
  }

  // Now go through all pairs of permutations and compare if all
  // positions stored in positions1 are swapped with the
  // positions stored in positions2 between the two permutations.
  // If yes, veto one of the two.
  for (auto perm_itr1 = m_permutation_list.begin(); perm_itr1 != m_permutation_list.end(); ++perm_itr1) {
    // Skip if permutation is already vetoed.
    if (perm_itr1->vetoed) continue;

    for (auto perm_itr2 = perm_itr1 + 1; perm_itr2 != m_permutation_list.end(); ++perm_itr2) {
      // Skip if permutation is already vetoed.
      if (perm_itr2->vetoed) continue;

      // For convenience, store const pointers to the two
      // containers that we need to compare against each other.
      const std::vector<int>* cont1;
      const std::vector<int>* cont2;

      if (ptype == Particles::Type::kParton) {
        cont1 = &perm_itr1->partons;
        cont2 = &perm_itr2->partons;
      } else if (ptype == Particles::Type::kElectron) {
        cont1 = &perm_itr1->electrons;
        cont2 = &perm_itr2->electrons;
      } else if (ptype == Particles::Type::kMuon) {
        cont1 = &perm_itr1->muons;
        cont2 = &perm_itr2->muons;
      } else if (ptype == Particles::Type::kPhoton) {
        cont1 = &perm_itr1->photons;
        cont2 = &perm_itr2->photons;
      } else if (ptype == Particles::Type::kTrack) {
        cont1 = &perm_itr1->tracks;
        cont2 = &perm_itr2->tracks;
      }

      // If we find as many matches as stored indices, it's a
      // full match and permutation 2 is to be vetoed.
      size_t matched{0};
      for (size_t i = 0; i < positions1.size(); ++i) {
        if (cont1->at(positions1[i]) == cont2->at(positions2[i])) ++matched;
      }
      if (matched == positions1.size()) perm_itr2->vetoed = true;
    }
  }

  // Return with no errors.
  return 1;
}

// ---------------------------------------------------------
int PermutationHandler::RemoveParticlePermutations(Particles::Type ptype, int index, int position) {
  // check if particles are defined
  if (!CheckParticles())
    return 0;

  // Veto all permutations, where 'index' is found in the
  // specified position of the particle container.
  for (auto& perm : m_permutation_list) {
    if (ptype == Particles::Type::kParton && perm.partons.at(position) == index) perm.vetoed = true;
    if (ptype == Particles::Type::kElectron && perm.electrons.at(position) == index) perm.vetoed = true;
    if (ptype == Particles::Type::kMuon && perm.muons.at(position) == index) perm.vetoed = true;
    if (ptype == Particles::Type::kPhoton && perm.photons.at(position) == index) perm.vetoed = true;
    if (ptype == Particles::Type::kTrack && perm.tracks.at(position) == index) perm.vetoed = true;
  }

  // return error code;
  return 1;
}

// ---------------------------------------------------------
int PermutationHandler::CheckParticles() {
  // check if particles are defined
  if (!m_particles) {
    std::cout << "KLFitter::PermutationHandler::CheckParticles(). Pointer to particles not defined yet." << std::endl;
    return 0;
  }

  if (!(*m_particles)) {
    std::cout << "KLFitter::PermutationHandler::CheckParticles(). Particles not defined yet." << std::endl;
    return 0;
  }

  // no error
  return 1;
}
}  // namespace KLFitter
