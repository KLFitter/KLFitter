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
bool Permutations::IsVetoed(int index) {
  return m_permutation_list.at(index).vetoed;
}

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
int Permutations::Reset() {
  // Clear particle and permutation tables.
  m_particles_table.clear();
  m_permutation_list.clear();

  // no error
  return 1;
}

// ---------------------------------------------------------
int Permutations::InvariantParticlePermutations(Particles::Type ptype, std::vector<int> indices) {
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

  // for (int iperm1 = nperm-1; iperm1 >= 1; --iperm1) {
  //   int offset = 0;
  //   for (Particles::Type itype = Particles::Type::kParton; itype < ptype; ++itype)
  //     offset += (*m_particles)->NParticles(itype);

  //   // get permutation
  //   const std::vector<int> permutation1 = m_permutation_table[iperm1];

  //   // Count numbers of removed permutations to adjust the index for the
  //   // for-loop over permutation 1.
  //   unsigned int removed_perms = 0;

  //   for (int iperm2 = iperm1-1; iperm2 >= 0; --iperm2) {
  //     // get second permutation
  //     const std::vector<int> permutation2 = m_permutation_table[iperm2];

  //     // loop over index vectors
  //     unsigned int numberOfInvariantMatches(0);

  //     for (unsigned int i = 0, I = indexVectorPosition1.size(); i < I; i++)  {
  //       int indexPosition1 = indexVectorPosition1[i] + offset;
  //       int indexPosition2 = indexVectorPosition2[i] + offset;

  //       // check indices
  //       if (permutation1[indexPosition1] == permutation2[indexPosition2] && permutation1[indexPosition2] == permutation2[indexPosition1])
  //         numberOfInvariantMatches++;
  //     }

  //     if (numberOfInvariantMatches == indexVectorPosition1.size()) {
  //       // m_permutation_table.erase(m_permutation_table.begin() + iperm2);

  //       // m_particles_table.erase(m_particles_table.begin() + iperm2);

  //       removed_perms++;
  //     }
  //   }  // second permutation

  //   // Decrement the first permutation index by the number of removed
  //   // permutations, otherwise we might go out of scope..
  //   iperm1 -= removed_perms;
  // }  // first permutation

  // return error code
  return err;
}

// ---------------------------------------------------------
int Permutations::RemoveParticlePermutations(Particles::Type ptype, int index, int position) {
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
