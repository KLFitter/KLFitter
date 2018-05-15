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

// ---------------------------------------------------------
KLFitter::Permutations::Permutations(KLFitter::Particles ** p, KLFitter::Particles ** pp) :
    fParticles(p),
    fParticlesPermuted(pp),
    fPermutationIndex(-1) {
  // empty
}

// ---------------------------------------------------------
KLFitter::Permutations::~Permutations() = default;

// ---------------------------------------------------------
int KLFitter::Permutations::SetPermutation(int index) {
  // check if permutation table exists
  if (fParticlesTable.empty()) {
    std::cout << "KLFitter::Permutations::SetPermutation(). Table does not exist yet." << std::endl;
    return 0;
  }

  // check index
  if (index < 0 || index >= NPermutations()) {
    std::cout << "KLFitter::Permutations::SetPermutation(). Index out of range." << std::endl;
    return 0;
  }

  if (!fParticlesPermuted) {
    std::cout << "KLFitter::Permutations::SetPermutation(). Pointer to permuted particles not available." << std::endl;
    return 0;
  }

  // set permutation
  (*fParticlesPermuted) = &fParticlesTable[index];

  // set permutation index
  fPermutationIndex = index;

  // no error
  return 1;
}

// ---------------------------------------------------------
int KLFitter::Permutations::CreatePermutations(int nPartonsInPermutations) {
  // reset existing particle and permuation tables
  Reset();

  // create new table of particles
  fParticlesTable = std::vector<KLFitter::Particles>{};

  // create new table of permutations
  fPermutationTable = std::vector<std::vector<int> >{};

  // check particles
  CheckParticles();

  // get number of objects per category
  int npartons   = (*fParticles)->NPartons();
  int nelectrons = (*fParticles)->NElectrons();
  int nmuons     = (*fParticles)->NMuons();
  int nphotons     = (*fParticles)->NPhotons();

  bool isDilepton(false);

  if (nelectrons != 0 && (*fParticles)->LeptonCharge(0, KLFitter::Particles::kElectron) != -9)
    isDilepton = true;

  if (nmuons != 0 && (*fParticles)->LeptonCharge(0, KLFitter::Particles::kMuon) != -9)
    isDilepton = true;

  // create table for parton, electron, muon and photons permutations
  fTablePartons = std::vector<std::vector<int> >{};
  CreateSubTable(npartons, &fTablePartons, nPartonsInPermutations);

  fTableElectrons = std::vector<std::vector<int> >{};
  CreateSubTable(nelectrons, &fTableElectrons);

  fTableMuons = std::vector<std::vector<int> >{};
  CreateSubTable(nmuons, &fTableMuons);

  fTablePhotons = std::vector<std::vector<int> >{};
  CreateSubTable(nphotons, &fTablePhotons);

  int npartonsPerm = npartons;
  if (nPartonsInPermutations >= 0)
    npartonsPerm = nPartonsInPermutations;

  // get number of possible permutations for each category
  int npermpartons   = fTablePartons.size() <= 0 ? 1 : fTablePartons.size();
  int npermelectrons = fTableElectrons.size() <= 0 ? 1 : fTableElectrons.size();
  int npermmuons     = fTableMuons.size() <= 0 ? 1 : fTableMuons.size();
  int npermphotons     = fTablePhotons.size() <= 0 ? 1 : fTablePhotons.size();
  int npermoverall   = npartonsPerm + nelectrons + nmuons + nphotons;

  // loop over all parton permutations
  for (int ipermparton = 0; ipermparton < npermpartons; ++ipermparton) {
    // loop over all electron permutations
    for (int ipermelectron = 0; ipermelectron < npermelectrons; ++ipermelectron) {
      // loop over all muon permutations
      for (int ipermmuon = 0; ipermmuon < npermmuons; ++ipermmuon) {
        // loop over all photon permutations
        for (int ipermphoton = 0; ipermphoton < npermphotons; ++ipermphoton) {
          // create new particles
          KLFitter::Particles particles{};

          // create new permutation
          std::vector<int> permutation(npermoverall);

          // loop over all partons
          for (int i = 0; i < npartonsPerm; ++i) {
            // get index
            int index = fTablePartons[ipermparton][i];

            // add parton
            particles.AddParticle((*fParticles)->Parton(index),
                                  (*fParticles)->DetEta(index, KLFitter::Particles::kParton),
                                  KLFitter::Particles::kParton,
                                  (*fParticles)->NameParticle(index, KLFitter::Particles::kParton),
                                  (*fParticles)->JetIndex(index),
                                  (*fParticles)->IsBTagged(index),
                                  (*fParticles)->BTaggingEfficiency(index),
                                  (*fParticles)->BTaggingRejection(index),
                                  (*fParticles)->TrueFlavor(index),
                                  (*fParticles)->BTagWeight(index));

            // set permutation
            permutation[i] = index;
          }

          // loop over all electrons
          for (int i = 0; i < nelectrons; ++i) {
            // get index
            int index = fTableElectrons[ipermelectron][i];

            // if isDilepton include charge of the lepton
            if (isDilepton) {
              // add electron
              particles.AddParticle((*fParticles)->Electron(index),
                                    (*fParticles)->DetEta(index, KLFitter::Particles::kElectron),
                                    (*fParticles)->LeptonCharge(index, KLFitter::Particles::kElectron),
                                    KLFitter::Particles::kElectron,
                                    (*fParticles)->NameParticle(index, KLFitter::Particles::kElectron),
                                    (*fParticles)->ElectronIndex(index));
            } else {
              // add electron
              particles.AddParticle((*fParticles)->Electron(index),
                                    (*fParticles)->DetEta(index, KLFitter::Particles::kElectron),
                                    KLFitter::Particles::kElectron,
                                    (*fParticles)->NameParticle(index, KLFitter::Particles::kElectron),
                                    (*fParticles)->ElectronIndex(index));
            }

            // set permutation
            permutation[npartonsPerm + i] = index;
          }

          // loop over all muons
          for (int i = 0; i < nmuons; ++i) {
            // get index
            int index = fTableMuons[ipermmuon][i];

            // if isDilepton include charge of the lepton
            if (isDilepton) {
              // add muon
              particles.AddParticle((*fParticles)->Muon(index),
                                    (*fParticles)->DetEta(index, KLFitter::Particles::kMuon),
                                    (*fParticles)->LeptonCharge(index, KLFitter::Particles::kMuon),
                                    KLFitter::Particles::kMuon,
                                    (*fParticles)->NameParticle(index, KLFitter::Particles::kMuon),
                                    (*fParticles)->MuonIndex(index));
            } else {
              // add muon
              particles.AddParticle((*fParticles)->Muon(index),
                                    (*fParticles)->DetEta(index, KLFitter::Particles::kMuon),
                                    KLFitter::Particles::kMuon,
                                    (*fParticles)->NameParticle(index, KLFitter::Particles::kMuon),
                                    (*fParticles)->MuonIndex(index));
            }

            // set permutation
            permutation[npartonsPerm + nelectrons + i] = index;
          }

          // loop over all photons
          for (int i = 0; i < nphotons; ++i) {
            // get index
            int index = fTablePhotons[ipermphoton][i];

            // add photon
            particles.AddParticle((*fParticles)->Photon(index),
                                  (*fParticles)->DetEta(index, KLFitter::Particles::kPhoton),
                                  KLFitter::Particles::kPhoton,
                                  (*fParticles)->NameParticle(index, KLFitter::Particles::kPhoton),
                                  (*fParticles)->PhotonIndex(index));

            // set permutation
            permutation[npartonsPerm + nelectrons + nmuons + i] = index;
          }

          // add particles to table
          fParticlesTable.emplace_back(particles);

          // add permutation to table
          fPermutationTable.emplace_back(permutation);
        }
      }
    }
  }

  // no error
  return 1;
}

// ---------------------------------------------------------
int KLFitter::Permutations::Reset() {
  // Clear particle and permutation tables.
  fParticlesTable.clear();
  fPermutationTable.clear();

  // no error
  return 1;
}

// ---------------------------------------------------------
int KLFitter::Permutations::CreateSubTable(int Nobj, std::vector<std::vector<int> >* table, int Nmax) {
  if (Nmax < 0) {
    std::vector<int> vidx;
    for (int i(0); i < Nobj; ++i) {
      vidx.push_back(i);
    }

    do {
      table->emplace_back(std::vector<int>(vidx));
    } while (std::next_permutation(vidx.begin(), vidx.end()));
  } else {
    std::vector<std::vector<int> > v = Get_M_from_N(Nobj, Nmax);

    for (unsigned int i(0), n(v.size()); i < n; ++i) {
      std::vector<int> vidx = v[i];
      do {
        table->emplace_back(std::vector<int>(vidx));
      } while (std::next_permutation(vidx.begin(), vidx.end()));
    }
  }

  // no error
  return 1;
}

// ---------------------------------------------------------
int KLFitter::Permutations::InvariantParticlePermutations(KLFitter::Particles::ParticleType ptype, std::vector<int> indexVector) {
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
    if (index < 0 || index >= (*fParticles)->NParticles(ptype)) {
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

  // check particles table
  if (fParticlesTable.empty()) {
    std::cout << "KLFitter::Permutations::InvariantParticlePermutations(). Table does not exist yet." << std::endl;
    return 0;
  }

  // check permutation table
  if (fPermutationTable.empty()) {
    std::cout << "KLFitter::Permutations::InvariantParticlePermutations(). Table of parton permutations doesn ot exist." << std::endl;
    return 0;
  }

  // no error
  int err = 1;

  // loop over all permutations (if there are only 2 indices left)
  if (indexVector.size() == 2) {
    // get number of permutations
    int nperm = NPermutations();

    for (int iperm = nperm-1; iperm >= 0; --iperm) {
      int offset = 0;
      for (KLFitter::Particles::ParticleType itype = KLFitter::Particles::kParton; itype < ptype; ++itype)
        offset += (*fParticles)->NParticles(itype);
      int index1 = indexVector[0] + offset;
      int index2 = indexVector[1] + offset;

      // get permutation
      const std::vector<int>& permutation = fPermutationTable[iperm];

      // check indices
      if (permutation[index1] >= permutation[index2]) {
        fPermutationTable.erase(fPermutationTable.begin() + iperm);

        fParticlesTable.erase(fParticlesTable.begin() + iperm);
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
int KLFitter::Permutations::InvariantParticleGroupPermutations(KLFitter::Particles::ParticleType ptype, std::vector<int> indexVectorPosition1,  std::vector<int> indexVectorPosition2) {
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
    if (index1 < 0 || index1 >= (*fParticles)->NParticles(ptype)) {
      std::cout << "KLFitter::Permutations::InvariantParticleGroupPermutations(). Index out of range." << index1 << " " << (*fParticles)->NParticles(ptype) << std::endl;
      return 0;
    }
    int index2 = indexVectorPosition2[i];
    if (index2 < 0 || index2 >= (*fParticles)->NParticles(ptype)) {
      std::cout << "KLFitter::Permutations::InvariantParticleGroupPermutations(). Index out of range." << index2 << " " << (*fParticles)->NParticles(ptype) << std::endl;
      return 0;
    }
  }

  // swap indices
  indexVectorPosition1.clear();
  std::set<int>::iterator it_indexSetPosition1Begin = indexSetPosition1.begin();
  std::set<int>::iterator it_indexSetPosition1End   = indexSetPosition1.end();
  for (; it_indexSetPosition1Begin != it_indexSetPosition1End; it_indexSetPosition1Begin++)
    indexVectorPosition1.push_back(*it_indexSetPosition1Begin);

  // check particles table
  if (fParticlesTable.empty()) {
    std::cout << "KLFitter::Permutations::InvariantParticleGroupPermutations(). Table does not exist yet." << std::endl;
    return 0;
  }

  // check permutation table
  if (fPermutationTable.empty()) {
    std::cout << "KLFitter::Permutations::InvariantParticleGroupPermutations(). Table of parton permutations doesn ot exist." << std::endl;
    return 0;
  }

  // no error
  int err = 1;

  // loop over all permutations
  // get number of permutations
  int nperm = NPermutations();

  for (int iperm1 = nperm-1; iperm1 >= 1; --iperm1) {
    int offset = 0;
    for (KLFitter::Particles::ParticleType itype = KLFitter::Particles::kParton; itype < ptype; ++itype)
      offset += (*fParticles)->NParticles(itype);

    // get permutation
    const std::vector<int>& permutation1 = fPermutationTable[iperm1];

    for (int iperm2 = iperm1-1; iperm2 >= 0; --iperm2) {
      // get second permutation
      const std::vector<int>& permutation2 = fPermutationTable[iperm2];

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
        fPermutationTable.erase(fPermutationTable.begin() + iperm2);

        fParticlesTable.erase(fParticlesTable.begin() + iperm2);
      }
    }  // second permutation
  }  // first permutation

  // return error code
  return err;
}

// ---------------------------------------------------------
int KLFitter::Permutations::RemoveParticlePermutations(KLFitter::Particles::ParticleType ptype, int index, int position) {
  // check if particles are defined
  if (!CheckParticles())
    return 0;

  // check index
  if (index < 0 || index >= (*fParticles)->NParticles(ptype)) {
    std::cout << "KLFitter::Permutations::RemoveParticlePermutations(). Index out of range." << std::endl;
    return 0;
  }

  // check particles table
  if (fParticlesTable.empty()) {
    std::cout << "KLFitter::Permutations::RemoveParticlePermutations(). Table does not exist yet." << std::endl;
    return 0;
  }

  // check permutation table
  if (fPermutationTable.empty()) {
    std::cout << "KLFitter::Permutations::RemoveParticlePermutations(). Table of parton permutations does not exist." << std::endl;
    return 0;
  }

  // get offset for the particle type
  int offset = 0;
  for (KLFitter::Particles::ParticleType itype = KLFitter::Particles::kParton; itype < ptype; ++itype)
    offset += (*fParticles)->NParticles(itype);
  position += offset;

  // loop over all permutations
  for (int iPerm(NPermutations()-1); iPerm >= 0; --iPerm) {
    const std::vector<int>& permutation = fPermutationTable[iPerm];

    if (permutation[position] == index) {
      fPermutationTable.erase(fPermutationTable.begin() + iPerm);

      fParticlesTable.erase(fParticlesTable.begin() + iPerm);
    }
  }

  // return error code;
  return 1;
}

// ---------------------------------------------------------
int KLFitter::Permutations::CheckParticles() {
  // check if particles are defined
  if (!fParticles) {
    std::cout << "KLFitter::Permutations::CheckParticles(). Pointer to particles not defined yet." << std::endl;
    return 0;
  }

  if (!(*fParticles)) {
    std::cout << "KLFitter::Permutations::CheckParticles(). Particles not defined yet." << std::endl;
    return 0;
  }

  // no error
  return 1;
}

// ---------------------------------------------------------
std::vector<std::vector<int> >* KLFitter::Permutations::PermutationTable() {
  return &fPermutationTable;
}

// ---------------------------------------------------------
std::vector<int> KLFitter::Permutations::Get_int_vector(int i) {
  std::vector<int> vtmp;
  vtmp.push_back(i);
  return vtmp;
}

// ---------------------------------------------------------
std::vector<int> KLFitter::Permutations::Get_int_plus_vector(int i, std::vector<int> v) {
  std::vector<int> vtmp = Get_int_vector(i);
  for (unsigned int j(0), jend(v.size()); j < jend; ++j) {
    vtmp.push_back(v[j]);
  }
  return vtmp;
}

// ---------------------------------------------------------
std::vector<std::vector<int> > KLFitter::Permutations::Get_M_from_N(unsigned int N, unsigned int M, unsigned int start) {
  std::vector<std::vector<int> > v(0);
  for (unsigned int i(start); i < N-(M-1); ++i) {
    if (M == 1) {
      v.push_back(Get_int_vector(i));
    } else {
      std::vector<std::vector<int> > vnext = Get_M_from_N(N, M-1, i+1);
      for (unsigned int j(0), jend(vnext.size()); j < jend; ++j) {
        v.push_back(Get_int_plus_vector(i, vnext[j]));
      }
    }
  }
  return v;
}
