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
#include <set>

namespace {
// ---------------------------------------------------------
std::vector<int> get_int_vector(int i) {
  std::vector<int> vtmp;
  vtmp.push_back(i);
  return vtmp;
}

// ---------------------------------------------------------
std::vector<int> get_int_plus_vector(int i, std::vector<int> v) {
  std::vector<int> vtmp = get_int_vector(i);
  for (unsigned int j(0), jend(v.size()); j < jend; ++j) {
    vtmp.push_back(v[j]);
  }
  return vtmp;
}

// ---------------------------------------------------------
std::vector<std::vector<int> > get_m_from_n(unsigned int N, unsigned int M, unsigned int start = 0) {
  std::vector<std::vector<int> > v(0);
  for (unsigned int i(start); i < N-(M-1); ++i) {
    if (M == 1) {
      v.push_back(get_int_vector(i));
    } else {
      std::vector<std::vector<int> > vnext = get_m_from_n(N, M-1, i+1);
      for (unsigned int j(0), jend(vnext.size()); j < jend; ++j) {
        v.push_back(get_int_plus_vector(i, vnext[j]));
      }
    }
  }
  return v;
}

// ---------------------------------------------------------
void create_subtable(int Nobj, std::vector<std::vector<int> >* table, int Nmax = -1) {
  if (Nmax < 0) {
    std::vector<int> vidx;
    for (int i(0); i < Nobj; ++i) {
      vidx.push_back(i);
    }

    do {
      table->emplace_back(std::vector<int>(vidx));
    } while (std::next_permutation(vidx.begin(), vidx.end()));
  } else {
    std::vector<std::vector<int> > v = get_m_from_n(Nobj, Nmax);

    for (unsigned int i(0), n(v.size()); i < n; ++i) {
      std::vector<int> vidx = v[i];
      do {
        table->emplace_back(std::vector<int>(vidx));
      } while (std::next_permutation(vidx.begin(), vidx.end()));
    }
  }
}
}  // namespace


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

  bool isDilepton(false);

  if (nelectrons != 0 && (*m_particles)->electrons.at(0).GetCharge() != -9)
    isDilepton = true;

  if (nmuons != 0 && (*m_particles)->muons.at(0).GetCharge() != -9)
    isDilepton = true;

  // create table for parton, electron, muon, photon, and track's permutations
  m_table_partons = std::vector<std::vector<int> >{};
  create_subtable(npartons, &m_table_partons);

  m_table_electrons = std::vector<std::vector<int> >{};
  create_subtable(nelectrons, &m_table_electrons);

  m_table_muons = std::vector<std::vector<int> >{};
  create_subtable(nmuons, &m_table_muons);

  m_table_photons = std::vector<std::vector<int> >{};
  create_subtable(nphotons, &m_table_photons);

  m_table_tracks = std::vector<std::vector<int> >{};
  create_subtable(ntracks, &m_table_tracks);

  int npartonsPerm = npartons;

  // get number of possible permutations for each category
  int npermpartons   = m_table_partons.size() <= 0 ? 1 : m_table_partons.size();
  int npermelectrons = m_table_electrons.size() <= 0 ? 1 : m_table_electrons.size();
  int npermmuons     = m_table_muons.size() <= 0 ? 1 : m_table_muons.size();
  int npermphotons     = m_table_photons.size() <= 0 ? 1 : m_table_photons.size();
  int npermtracks    = m_table_tracks.size() <= 0 ? 1 : m_table_tracks.size();
  int npermoverall   = npartonsPerm + nelectrons + nmuons + nphotons + ntracks;

  // loop over all parton permutations
  for (int ipermparton = 0; ipermparton < npermpartons; ++ipermparton) {
    // loop over all electron permutations
    for (int ipermelectron = 0; ipermelectron < npermelectrons; ++ipermelectron) {
      // loop over all muon permutations
      for (int ipermmuon = 0; ipermmuon < npermmuons; ++ipermmuon) {
        // loop over all photon permutations
        for (int ipermphoton = 0; ipermphoton < npermphotons; ++ipermphoton) {
          // loop over all track permutations
          for (int ipermtrack = 0; ipermtrack < npermtracks; ++ipermtrack) {
            // create new particles
            ParticleCollection particles{};

            // create new permutation
            std::vector<int> permutation(npermoverall);

            // loop over all partons
            for (int i = 0; i < npartonsPerm; ++i) {
              // get index
              int index = m_table_partons[ipermparton][i];

              // add parton
              particles.AddParticle((*m_particles)->partons.at(index));

              // set permutation
              permutation[i] = index;
            }

            // loop over all electrons
            for (size_t i = 0; i < nelectrons; ++i) {
              // get index
              int index = m_table_electrons[ipermelectron][i];

              // if isDilepton include charge of the lepton
              if (isDilepton) {
                // add electron
                particles.AddParticle((*m_particles)->electrons.at(index));
              } else {
                // add electron
                particles.AddParticle((*m_particles)->electrons.at(index));
              }

              // set permutation
              permutation[npartonsPerm + i] = index;
            }

            // loop over all muons
            for (size_t i = 0; i < nmuons; ++i) {
              // get index
              int index = m_table_muons[ipermmuon][i];

              // if isDilepton include charge of the lepton
              if (isDilepton) {
                // add muon
                particles.AddParticle((*m_particles)->muons.at(index));
              } else {
                // add muon
                particles.AddParticle((*m_particles)->muons.at(index));
              }

              // set permutation
              permutation[npartonsPerm + nelectrons + i] = index;
            }

            // loop over all photons
            for (size_t i = 0; i < nphotons; ++i) {
              // get index
              int index = m_table_photons[ipermphoton][i];

              // add photon
              particles.AddParticle((*m_particles)->photons.at(index));

              // set permutation
              permutation[npartonsPerm + nelectrons + nmuons + i] = index;
            }

            // loop over all tracks
            for (size_t i = 0; i < ntracks; ++i) {
              // get index
              int index = m_table_tracks[ipermtrack][i];

              // add track
              particles.AddParticle((*m_particles)->tracks.at(index));

              // set permutation
              permutation[npartonsPerm + nelectrons + nmuons + nphotons + i] = index;
            }


            // add particles to table
            m_particles_table.emplace_back(particles);

            // add permutation to table
            m_permutation_table.emplace_back(permutation);
          }
        }
      }
    }
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
