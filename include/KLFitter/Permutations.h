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

#ifndef KLFITTER_PERMUTATIONS_H_
#define KLFITTER_PERMUTATIONS_H_

#include <vector>

#include "KLFitter/Particles.h"

// ---------------------------------------------------------

/**
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter {
/**
  * \class KLFitter::Permutations
  * \brief A class permuting jets, electrons, muons and photons.
  *
  * The class gets a pointer to the orignal set of particles and a
  * pointer to the currently used permutations. It can calculate all
  * permutations and created a table. The pointer of the current
  * permutation is set to the entry in the table.
  */
class Permutations final {
 public:
  /** \name Constructors and destructors */
  /* @{ */

  /**
    * The default constructor.
    * @param p A pointer to the pointer to the original set of particles.
    * @param pp A pointer to the pointer to the permutated set of particles.
    */
  Permutations(KLFitter::Particles** p, KLFitter::Particles** pp);

  /**
    * Copy constructor
    */
  Permutations(const Permutations& obj);

  /**
    * The default destructor.
    */
  ~Permutations();

  /**
    * Assignment operator
    */
  Permutations& operator=(const Permutations& obj);

  /* @} */
  /** \name Member functions (Get)  */
  /* @{ */

  /**
    * Return the original particles.
    * @return A pointer to the particles.
    */
  KLFitter::Particles* Particles() { return *fParticles; }

  /**
    * Return the current permutation of jets and leptons.
    * @return A pointer to the permuted particles.
    */
  KLFitter::Particles* ParticlesPermuted() { return *fParticlesPermuted; }

  /**
    * Return the permutation table.
    * @return A pointer to the permutation table.
    */
  std::vector<std::vector<int>* >* PermutationTable();

  /**
    * Return the number of permutations.
    */
  int NPermutations() { return static_cast<int>(fParticlesTable -> size()); }

  /**
    * Return the current permutation index.
    * @return The current permutation index.
    */
  int PermutationIndex() { return fPermutationIndex; }

  std::vector<std::vector<int>*>* TablePartons() { return fTablePartons; }

  std::vector<std::vector<int>*>* TableElectrons() { return fTableElectrons; }

  std::vector<std::vector<int>*>* TableMuons() { return fTableMuons; }

  std::vector<std::vector<int>*>* TablePhotons() { return fTablePhotons; }

  /* @} */
  /** \name Member functions (Set)  */
  /* @{ */

  /**
    * Set the original particles.
    * @param particles A set of particles.
    * @return An error code.
    */
  int SetParticles(KLFitter::Particles* particles);

  /**
    * Set the permutation.
    * @param index The permutation index.
    * @return An error code.
    */
  int SetPermutation(int index);

  /* @} */
  /** \name Member functions (misc)  */
  /* @{ */

  /**
    * Create all possible permutations of jets and leptons.
    * However, make permutations with exactly nPartonsInPermutations.
    */
  int CreatePermutations(int nPartonsInPermutations = -1);

  /**
    * Remove permutations in which all indices in the vector indexVector are exchanged
    * for the given particle type.
    * This is useful to reduce the number of permutations if
    * interchanging for example jets doesn't have any effect, e.g.,
    * if two jets come from a W (top).
    * @param ptype The type of the particle.
    * @param indexVector Vector of indices.
    * @return An error code.
    */
  int InvariantParticlePermutations(KLFitter::Particles::ParticleType ptype, std::vector<int> indexVector);

  /**
    * Remove permutations in which all indices in the vector indexVectorPosition1 are exchanged with the corresponding indices in indexVectorPosition2
    * for the given particle type.
    * This is useful to reduce the number of permutations if
    * interchanging a whole set of particles doesn't have any effect, e.g.,
    * the particles coming from the two hadronic top quarks in the fully hadronic channel.
    * @param ptype The type of the particle.
    * @param indexVectorPosition1 Vector of indices of first set of particle.
    * @param indexVectorPosition2 Vector of corresponding indices for second set of particle.
    * @return An error code.
    */
  int InvariantParticleGroupPermutations(KLFitter::Particles::ParticleType ptype, std::vector<int> indexVectorPosition1,  std::vector<int> indexVectorPosition2);

  /**
    * Remove permutations in which a certain particles is in a certain position.
    * This is useful to reduce the number of permutations if for example
    * a b-tagged jet is forbidden in the position of a light jet.
    * @param ptype The type of the particle.
    * @param index The index of the particle.
    * @param position The position in which it is forbidden.
    * @return An error code.
    */
  int RemoveParticlePermutations(KLFitter::Particles::ParticleType ptype, int index, int position);

  /**
    * Reset Permutations.
    * @return An error code.
    */
  int Reset();

  /**
    * Creates table of permutations.
    */
  int CreateSubTable(int Nobj,  std::vector<std::vector<int>*>* table, int Nmax = -1);

  /* @} */

 private:
  /**
    * Check if particles are defined.
    */
  int CheckParticles();

 private:
  /**
    * Helper functions to efficienctly create permutations of N particles of only M selected particles.
    */

  std::vector<int> Get_int_vector(int i);
  std::vector<int> Get_int_plus_vector(int i, std::vector<int> v);
  std::vector<std::vector<int> > Get_M_from_N(unsigned int N, unsigned int M, unsigned int start = 0);

  /**
    * A pointer to the pointer of original particles.
    */
  KLFitter::Particles** fParticles;

  /**
    * A pointer to the pointer of permuted particles.
    */
  KLFitter::Particles** fParticlesPermuted;

  /**
    * A table of permuted particles (jets and leptons).
    */
  std::vector <KLFitter::Particles*>* fParticlesTable;

  /**
    * A table of permutations. Needed for the math.
    */
  std::vector < std::vector <int>*>* fPermutationTable;

  /**
    * The permutation index
    */
  int fPermutationIndex;

  std::vector<std::vector<int>*>* fTablePartons;
  std::vector<std::vector<int>*>* fTableElectrons;
  std::vector<std::vector<int>*>* fTableMuons;
  std::vector<std::vector<int>*>* fTablePhotons;
};
}  // namespace KLFitter

#endif  // KLFITTER_PERMUTATIONS_H_
