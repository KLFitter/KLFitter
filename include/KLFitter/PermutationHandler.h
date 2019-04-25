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

#ifndef KLFITTER_PERMUTATIONHANDLER_H_
#define KLFITTER_PERMUTATIONHANDLER_H_

#include <deque>
#include <vector>

#include "KLFitter/ParticleCollection.h"

// ---------------------------------------------------------

namespace KLFitter {
/**
 * Simple structure to hold multiple lists, which in turn hold
 * indices of particles. These lists can then be used to
 * calculate all possible permutations.
 */
struct Permutation {
  std::deque<int> partons{};
  std::deque<int> electrons{};
  std::deque<int> muons{};
  std::deque<int> photons{};
  std::deque<int> tracks{};

  /// Boolean whether permutation is vetoed.
  bool vetoed{false};

  /// Calculate the next possible permutation.
  bool next_permutation();

  /// Print information about the permutation.
  std::string print();
};

/**
 * Class to calculate and provide permutations of particles. The
 * class gets a pointer to the orignal set of particles and a
 * pointer to the currently used permutations. It can calculate
 * all permutations and create a table. The pointer of the
 * current permutation is set to the entry in the table.
 */
class PermutationHandler final {
 public:
  /**
   * The default constructor.
   * @param p A pointer to the pointer to the original set of particles.
   * @param pp A pointer to the pointer to the permutated set of particles.
   */
  PermutationHandler(KLFitter::ParticleCollection** p, KLFitter::ParticleCollection** pp);

  /// The (defaulted) copy constructor.
  explicit PermutationHandler(const PermutationHandler& o);

  /// The (defaulted) destructor.
  ~PermutationHandler();

  /// The (defaulted) assignment operator.
  PermutationHandler& operator=(const PermutationHandler& obj);

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

  const std::deque<Permutation>& GetList() { return m_permutation_list; }

  /**
   * Return the number of permutations.
   */
  int NPermutations() const { return static_cast<int>(m_particles_table.size()); }

  /** @} */
  /** \name Member functions (Set)  */
  /** @{ */

  bool IsVetoed(int index);

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
   * @param indices Vector of indices.
   * @return An error code.
   */
  int InvariantParticlePermutations(Particles::Type ptype, std::vector<int> indices);

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

  /// A list of permutations. Needed for the math.
  std::deque<Permutation> m_permutation_list;

  /// The permutation index
  int m_permutation_index;
};
}  // namespace KLFitter

#endif  // KLFITTER_PERMUTATIONHANDLER_H_
