/*!
 * \class KLFitter::InterfaceBase
 * \brief A base class for interfacing input data. 
 * \author Kevin Kr&ouml;ninger
 * \version 1.3
 * \date 03.12.2009
 *
 * This class is a base class for interfaces to input data, e.g. to
 * ROOT trees.
 */

// --------------------------------------------------------- 

#ifndef __INTERFACEBASE__H
#define __INTERFACEBASE__H

// --------------------------------------------------------- 

#include "Particles.h" 

// --------------------------------------------------------- 

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{

	class InterfaceBase
	{
		
	public: 
		
		/** \name Constructors and destructors */ 
		/* @{ */ 
		
		/** 
		 * The default constructor. 
		 */ 
		InterfaceBase(); 
		
		/**
		 * The default destructor.
		 */
		virtual ~InterfaceBase(); 

		/* @} */
		/** \name Member functions (Get)  */
		/* @{ */

		/**
		 * Return a set of particles. 
		 * @return A pointer to the set of particles. 
		 */ 
		KLFitter::Particles * Particles()
			{ return fParticles; }; 

		KLFitter::Particles ** PParticles()
			{ return &fParticles; }; 

		/**
		 * Return a set of particles. 
		 * @return A pointer to the set of particles. 
		 */ 
		KLFitter::Particles * ParticlesTruth()
			{ return fParticlesTruth; }; 

		/**
		 * Return a pointer to a pointer of a set of particles. 
		 * @return A pointer to the pointer of a set of particles. 
		 */ 
		KLFitter::Particles ** PParticlesTruth()
			{ return &fParticlesTruth; }; 

		/**
		 * Return the number of events. 
		 * @return The number of events
		 */ 
		virtual int NEvents()
		{ return 0; }; 

		/**
		 * Get event weight (MC@NLO)
		 * @return Pointer to event weight.
		 */
		virtual double Weight()
			{ return fWeight; };

		/**
		 * Return the measured missing transverse energy. 
		 * @return The missing ET
		 */ 
		virtual double ET_miss()
		{ return 0; }; 
		
		/* @} */
		/** \name Member functions (Set)  */
		/* @{ */

		/* @} */
		/** \name Member functions (misc)  */
		/* @{ */

		/* @} */

	protected: 

		/**
		 * A set of particles. 
		 */ 
		KLFitter::Particles * fParticles; 

		/**
		 * A set of truth particles. 
		 */ 
		KLFitter::Particles * fParticlesTruth; 

		/**
		 * A container for four-vectors. This is needed to properly remove
		 * the particles.
		 */ 
		std::vector <TLorentzVector *> * fLorentzVectorContainer; 

		/**
		 * The event weight
		 */
		double fWeight;

	private: 

	}; 

} // namespace KLFitter 

// --------------------------------------------------------- 

#endif // __INTERFACEBASE__H

