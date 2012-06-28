/*!
 * \class KLFitter::Permutations
 * \brief A class permuting jets, electrons, muons and photons. 
 * \author Kevin Kr&ouml;ninger
 * \version 1.3
 * \date 03.12.2009
 *
 * The class gets a pointer to the orignal set of particles and a
 * pointer to the currently used permutations. It can calculate all
 * permutations and created a table. The pointer of the current
 * permutation is set to the entry in the table.
 */

// --------------------------------------------------------- 

#ifndef PERMUTATIONS
#define PERMUTATIONS

// --------------------------------------------------------- 

#include "Particles.h"

// --------------------------------------------------------- 

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{

  class Permutations
  {
                
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
     * The default destructor.
     */
    virtual ~Permutations(); 

    /* @} */
    /** \name Member functions (Get)  */
    /* @{ */

    /**
     * Return the original particles. 
     * @return A pointer to the particles. 
     */ 
    KLFitter::Particles* Particles()                    
      { return *fParticles; }; 

    /**
     * Return the current permutation of jets and leptons. 
     * @return A pointer to the permuted particles. 
     */ 
    KLFitter::Particles* ParticlesPermuted()
      { return *fParticlesPermuted; }; 

    /**
     * Return the permutation table.
     * @return A pointer to the permutation table. 
     */
    std::vector<std::vector<int>* >* PermutationTable();

    /**
     * Return the number of permutations. 
     */ 
    int NPermutations()
    { return int(fParticlesTable -> size()); }; 

    /**
     * Return the current permutation index.
     * @return The current permutation index.
     */ 
    int PermutationIndex()
    { return fPermutationIndex; }; 

    std::vector<std::vector<int>*>* TablePartons()
      { return fTablePartons; }; 

    std::vector<std::vector<int>*>* TableElectrons()
      { return fTableElectrons; }; 

    std::vector<std::vector<int>*>* TableMuons()
      { return fTableMuons; }; 

    std::vector<std::vector<int>*>* TablePhotons()
      { return fTablePhotons; }; 

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
     */ 
    int CreatePermutations(); 

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
    int CreateSubTable(int Nobj,  std::vector<std::vector<int>*>* table); 

    /* @} */

  private: 

    /**
     * Check if particles are defined. 
     */ 
    int CheckParticles(); 
                
  protected: 

  private: 

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

} // namespace KLFitter 

// --------------------------------------------------------- 

#endif 

