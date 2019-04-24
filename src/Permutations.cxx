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

#include "KLFitter/Permutations.h"

#include <algorithm>
#include <iostream>
#include <numeric>
#include <set>

namespace KLFitter {
// ---------------------------------------------------------
Permutations::Permutations(ParticleCollection ** p, ParticleCollection ** pp)
  : m_particles(p)
  , m_particles_permuted(pp)
  , m_permutation_index(-1) {
  // empty
}

// ---------------------------------------------------------
Permutations::Permutations(const Permutations& o) = default;

// ---------------------------------------------------------
Permutations::~Permutations() = default;

// ---------------------------------------------------------
Permutations& Permutations::operator=(const Permutations& obj) = default;

// ---------------------------------------------------------
int Permutations::SetPermutation(int index) {
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
int Permutations::CreatePermutations() {
  // reset existing particle and permuation tables
  Reset();

  // create new table of particles
  m_particles_table = std::vector<ParticleCollection>{};

  // create new table of permutations
  m_permutation_table = std::vector<std::vector<int> >{};

  // check particles
  CheckParticles();

  // get number of objects per category
  size_t npartons   = (*m_particles)->partons.size();
  size_t nelectrons = (*m_particles)->electrons.size();
  size_t nmuons     = (*m_particles)->muons.size();
  size_t nphotons     = (*m_particles)->photons.size();
  size_t ntracks    = (*m_particles)->tracks.size();

  int npermoverall = npartons + nelectrons + nmuons + nphotons + ntracks;

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
    // create new permutation
    std::vector<int> permutation(npermoverall);
    ParticleCollection particles{};

    // loop over all partons
    size_t index = 0;
    for (const auto& parton : perm.partons) {
      particles.AddParticle((*m_particles)->partons.at(parton));
      permutation.at(index) = parton;
      index++;
    }

    // loop over all electrons
    index = 0;
    for (const auto& electron : perm.electrons) {
      particles.AddParticle((*m_particles)->electrons.at(electron));
      permutation.at(npartons + index) = electron;
      index++;
    }

    // loop over all muons
    index = 0;
    for (const auto& muon : perm.muons) {
      particles.AddParticle((*m_particles)->muons.at(muon));
      permutation.at(npartons + nelectrons + index) = muon;
      index++;
    }

    // loop over all photons
    index = 0;
    for (const auto& photon : perm.photons) {
      particles.AddParticle((*m_particles)->photons.at(photon));
      permutation.at(npartons + nelectrons + nmuons + index) = photon;
      index++;
    }

    // loop over all tracks
    index = 0;
    for (const auto& track : perm.tracks) {
      particles.AddParticle((*m_particles)->tracks.at(track));
      permutation.at(npartons + nelectrons + nmuons + nphotons + index) = track;
      index++;
    }

    // add particles to table
    m_particles_table.emplace_back(particles);

    // add permutation to table
    m_permutation_table.emplace_back(permutation);
  }

  // no error
  return 1;
}

// ---------------------------------------------------------
int Permutations::Reset() {
  // Clear particle and permutation tables.
  m_particles_table.clear();
  m_permutation_table.clear();

  // no error
  return 1;
}

// ---------------------------------------------------------
int Permutations::InvariantParticlePermutations(Particles::Type ptype, std::vector<int> indexVector) {
  // check if particles are defined
  if (!CheckParticles())
    return 0;

  // check indices
  std::set<int> indexSet;
  for (unsigned int i = 0, I = indexVector.size(); i < I; i++) {
    if (indexSet.find(indexVector[i]) != indexSet.end()) {
      std::cout << "KLFitter::Permutations::InvariantParticlePermutations(). Indices have to be different." << std::endl;
      return 0;
    } else {
      indexSet.insert(indexVector[i]);
    }
  }

  for (unsigned int i = 0, I = indexVector.size(); i < I; i++) {
    int index = indexVector[i];
    if (index < 0 || static_cast<size_t>(index) >= (*m_particles)->NParticles(ptype)) {
      std::cout << "KLFitter::Permutations::InvariantParticlePermutations(). Index out of range." << std::endl;
      return 0;
    }
  }

  // swap indices
  indexVector.clear();
  std::set<int>::iterator it_indexSetBegin = indexSet.begin();
  std::set<int>::iterator it_indexSetEnd   = indexSet.end();
  for (; it_indexSetBegin != it_indexSetEnd; it_indexSetBegin++)
    indexVector.push_back(*it_indexSetBegin);

  // no error
  int err = 1;

  // loop over all permutations (if there are only 2 indices left)
  if (indexVector.size() == 2) {
    // get number of permutations
    int nperm = NPermutations();

    for (int iperm = nperm-1; iperm >= 0; --iperm) {
      int offset = 0;
      for (Particles::Type itype = Particles::Type::kParton; itype < ptype; ++itype)
        offset += (*m_particles)->NParticles(itype);
      int index1 = indexVector[0] + offset;
      int index2 = indexVector[1] + offset;

      // get permutation
      const std::vector<int>& permutation = m_permutation_table[iperm];

      // check indices
      if (permutation[index1] >= permutation[index2]) {
        m_permutation_table.erase(m_permutation_table.begin() + iperm);

        m_particles_table.erase(m_particles_table.begin() + iperm);
      }
    }
  } else {
    // repeat until there are only 2 indices left
    while (indexVector.size() >= 2) {
      int index2 = indexVector.back();
      for (unsigned int i = 0, I = indexVector.size()-1; i < I; i++) {
        int index1 = indexVector[i];
        std::vector<int> newIndexVector;
        newIndexVector.push_back(index1);
        newIndexVector.push_back(index2);
        err *= InvariantParticlePermutations(ptype, newIndexVector);
      }
      indexVector.erase(--indexVector.end());
    }
  }

  // return error code
  return err;
}

// ---------------------------------------------------------
int Permutations::InvariantParticleGroupPermutations(Particles::Type ptype, std::vector<int> indexVectorPosition1,  std::vector<int> indexVectorPosition2) {
  // check if particles are defined
  if (!CheckParticles())
    return 0;

  // check if index vectors have the same size
  if (indexVectorPosition1.size()!= indexVectorPosition2.size()) {
    std::cout << "KLFitter::Permutations::InvariantParticleGroupPermutations(). Index Vectors need to have the same size." << std::endl;
    return 0;
  }

  // check indices
  std::set<int> indexSetPosition1;
  std::set<int> indexSetPosition2;
  for (unsigned int i = 0, I = indexVectorPosition1.size(); i < I; i++) {
    if (indexSetPosition1.find(indexVectorPosition1[i]) != indexSetPosition1.end()) {
      std::cout << "KLFitter::Permutations::InvariantParticleGroupPermutations(). Indices within same index vector have to be different." << std::endl;
      return 0;
    } else if (indexSetPosition2.find(indexVectorPosition2[i]) != indexSetPosition2.end()) {
      std::cout << "KLFitter::Permutations::InvariantParticleGroupPermutations(). Indices within same index vector have to be different." << std::endl;
      return 0;
    } else if (indexVectorPosition1[i] == indexVectorPosition2[i]) {
      std::cout << "KLFitter::Permutations::InvariantParticleGroupPermutations(). Indices have to be different." << std::endl;
      return 0;
    } else {
      indexSetPosition1.insert(indexVectorPosition1[i]);
      indexSetPosition2.insert(indexVectorPosition2[i]);
    }
  }

  for (unsigned int i = 0, I = indexVectorPosition1.size(); i < I; i++) {
    int index1 = indexVectorPosition1[i];
    if (index1 < 0 || static_cast<size_t>(index1) >= (*m_particles)->NParticles(ptype)) {
      std::cout << "KLFitter::Permutations::InvariantParticleGroupPermutations(). Index out of range." << index1 << " " << (*m_particles)->NParticles(ptype) << std::endl;
      return 0;
    }
    int index2 = indexVectorPosition2[i];
    if (index2 < 0 || static_cast<size_t>(index2) >= (*m_particles)->NParticles(ptype)) {
      std::cout << "KLFitter::Permutations::InvariantParticleGroupPermutations(). Index out of range." << index2 << " " << (*m_particles)->NParticles(ptype) << std::endl;
      return 0;
    }
  }

  // swap indices
  indexVectorPosition1.clear();
  std::set<int>::iterator it_indexSetPosition1Begin = indexSetPosition1.begin();
  std::set<int>::iterator it_indexSetPosition1End   = indexSetPosition1.end();
  for (; it_indexSetPosition1Begin != it_indexSetPosition1End; it_indexSetPosition1Begin++)
    indexVectorPosition1.push_back(*it_indexSetPosition1Begin);

  // no error
  int err = 1;

  // loop over all permutations
  // get number of permutations
  int nperm = NPermutations();

  for (int iperm1 = nperm-1; iperm1 >= 1; --iperm1) {
    int offset = 0;
    for (Particles::Type itype = Particles::Type::kParton; itype < ptype; ++itype)
      offset += (*m_particles)->NParticles(itype);

    // get permutation
    const std::vector<int> permutation1 = m_permutation_table[iperm1];

    // Count numbers of removed permutations to adjust the index for the
    // for-loop over permutation 1.
    unsigned int removed_perms = 0;

    for (int iperm2 = iperm1-1; iperm2 >= 0; --iperm2) {
      // get second permutation
      const std::vector<int> permutation2 = m_permutation_table[iperm2];

      // loop over index vectors
      unsigned int numberOfInvariantMatches(0);

      for (unsigned int i = 0, I = indexVectorPosition1.size(); i < I; i++)  {
        int indexPosition1 = indexVectorPosition1[i] + offset;
        int indexPosition2 = indexVectorPosition2[i] + offset;

        // check indices
        if (permutation1[indexPosition1] == permutation2[indexPosition2] && permutation1[indexPosition2] == permutation2[indexPosition1])
          numberOfInvariantMatches++;
      }

      if (numberOfInvariantMatches == indexVectorPosition1.size()) {
        m_permutation_table.erase(m_permutation_table.begin() + iperm2);

        m_particles_table.erase(m_particles_table.begin() + iperm2);

        removed_perms++;
      }
    }  // second permutation

    // Decrement the first permutation index by the number of removed
    // permutations, otherwise we might go out of scope..
    iperm1 -= removed_perms;
  }  // first permutation

  // return error code
  return err;
}

// ---------------------------------------------------------
int Permutations::RemoveParticlePermutations(Particles::Type ptype, int index, int position) {
  // check if particles are defined
  if (!CheckParticles())
    return 0;

  // check index
  if (index < 0 || static_cast<size_t>(index) >= (*m_particles)->NParticles(ptype)) {
    std::cout << "KLFitter::Permutations::RemoveParticlePermutations(). Index out of range." << std::endl;
    return 0;
  }

  // get offset for the particle type
  int offset = 0;
  for (Particles::Type itype = Particles::Type::kParton; itype < ptype; ++itype)
    offset += (*m_particles)->NParticles(itype);
  position += offset;

  // loop over all permutations
  for (int iPerm(NPermutations()-1); iPerm >= 0; --iPerm) {
    const std::vector<int>& permutation = m_permutation_table[iPerm];

    if (permutation[position] == index) {
      m_permutation_table.erase(m_permutation_table.begin() + iPerm);

      m_particles_table.erase(m_particles_table.begin() + iPerm);
    }
  }

  // return error code;
  return 1;
}

// ---------------------------------------------------------
int Permutations::CheckParticles() {
  // check if particles are defined
  if (!m_particles) {
    std::cout << "KLFitter::Permutations::CheckParticles(). Pointer to particles not defined yet." << std::endl;
    return 0;
  }

  if (!(*m_particles)) {
    std::cout << "KLFitter::Permutations::CheckParticles(). Particles not defined yet." << std::endl;
    return 0;
  }

  // no error
  return 1;
}
}  // namespace KLFitter
