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

#include "KLFitter/ParticleCollection.h"

// ---------------------------------------------------------

namespace KLFitter {
/**
 * Class to calculate and provide permutations of particles. The
 * class gets a pointer to the orignal set of particles and a
 * pointer to the currently used permutations. It can calculate
 * all permutations and create a table. The pointer of the
 * current permutation is set to the entry in the table.
 */
class Permutations final {
 public:
  /**
   * The default constructor.
   * @param p A pointer to the pointer to the original set of particles.
   * @param pp A pointer to the pointer to the permutated set of particles.
   */
  Permutations(KLFitter::ParticleCollection** p, KLFitter::ParticleCollection** pp);

  /// The (defaulted) copy constructor.
  explicit Permutations(const Permutations& o);

  /// The (defaulted) destructor.
  ~Permutations();

  /// The (defaulted) assignment operator.
  Permutations& operator=(const Permutations& obj);

  /** \name Member functions (Get)  */
  /** @{ */

  /**
   * Return the original particles.
   * @return A pointer to the particles.
   */
  const KLFitter::ParticleCollection* Particles() const { return *m_particles; }

  /**
   * Return the current permutation of jets and leptons.
   * @return A pointer to the permuted particles.
   */
  KLFitter::ParticleCollection* ParticlesPermuted() { return *m_particles_permuted; }

  /**
   * Return the number of permutations.
   */
  int NPermutations() const { return static_cast<int>(m_particles_table.size()); }

  /**
   * Return the current permutation index.
   * @return The current permutation index.
   */
  int PermutationIndex() const { return m_permutation_index; }

  /** @} */
  /** \name Member functions (Set)  */
  /** @{ */

  /**
   * Set the original particles.
   * @param particles A set of particles.
   * @return An error code.
   */
  int SetParticles(KLFitter::ParticleCollection* particles);

  /**
   * Set the permutation.
   * @param index The permutation index.
   * @return An error code.
   */
  int SetPermutation(int index);

  /** @} */
  /** \name Member functions (misc)  */
  /** @{ */

  /**
   * Create all possible permutations of jets and leptons.
   */
  int CreatePermutations();

  /**
   * Remove permutations in which all indices in the vector
   * indexVector are exchanged for the given particle type. This
   * is useful to reduce the number of permutations if
   * interchanging for example jets doesn't have any effect,
   * e.g., if two jets come from a W (top).
   * @param ptype The type of the particle.
   * @param indexVector Vector of indices.
   * @return An error code.
   */
  int InvariantParticlePermutations(Particles::Type ptype, std::vector<int> indexVector);

  /**
   * Remove permutations in which all indices in the vector
   * indexVectorPosition1 are exchanged with the corresponding
   * indices in indexVectorPosition2 for the given particle type.
   * This is useful to reduce the number of permutations if
   * interchanging a whole set of particles doesn't have any
   * effect, e.g., the particles coming from the two hadronic top
   * quarks in the fully hadronic channel.
   * @param ptype The type of the particle.
   * @param indexVectorPosition1 Vector of indices of first set of particle.
   * @param indexVectorPosition2 Vector of corresponding indices for second set of particle.
   * @return An error code.
   */
  int InvariantParticleGroupPermutations(Particles::Type ptype, std::vector<int> indexVectorPosition1,  std::vector<int> indexVectorPosition2);

  /**
   * Remove permutations in which a certain particles is in a certain position.
   * This is useful to reduce the number of permutations if for example
   * a b-tagged jet is forbidden in the position of a light jet.
   * @param ptype The type of the particle.
   * @param index The index of the particle.
   * @param position The position in which it is forbidden.
   * @return An error code.
   */
  int RemoveParticlePermutations(Particles::Type ptype, int index, int position);

  /**
   * Reset Permutations.
   * @return An error code.
   */
  int Reset();

  /** @} */

 private:
  /// Check if particles are defined.
  int CheckParticles();

  /// A pointer to the pointer of original particles.
  KLFitter::ParticleCollection** m_particles;

  /// A pointer to the pointer of permuted particles.
  KLFitter::ParticleCollection** m_particles_permuted;

  /// A table of permuted particles (jets and leptons).
  std::vector<KLFitter::ParticleCollection> m_particles_table;

  /// A table of permutations. Needed for the math.
  std::vector<std::vector<int> > m_permutation_table;

  /// The permutation index
  int m_permutation_index;

  std::vector<std::vector<int> > m_table_partons;
  std::vector<std::vector<int> > m_table_electrons;
  std::vector<std::vector<int> > m_table_muons;
  std::vector<std::vector<int> > m_table_photons;
  std::vector<std::vector<int> > m_table_tracks;
};
}  // namespace KLFitter

#endif  // KLFITTER_PERMUTATIONS_H_
