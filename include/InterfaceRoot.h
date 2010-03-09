/*!
 * \class KLFitter::InterfaceRoot
 * \brief A class for interfacing a Root file. 
 * \author Kevin Kr&ouml;ninger
 * \version 1.3
 * \date 03.12.2009
 *
 * This class interfaces to a Root file.
 */

// --------------------------------------------------------- 

#ifndef __INTERFACEROOT__H
#define __INTERFACEROOT__H

// --------------------------------------------------------- 

#include "PREPROC.h"
#include "InterfaceBase.h" 
#include "Particles.h" 
#include "TFile.h"
#include "TTree.h" 

// --------------------------------------------------------- 

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{

	class InterfaceRoot : public InterfaceBase
	{
		
	public: 
		
		/** \name Constructors and destructors */ 
		/* @{ */ 
		
		/** 
		 * The default constructor. 
		 */ 
		InterfaceRoot(); 
		
		/**
		 * The default destructor.
		 */
		virtual ~InterfaceRoot(); 

		/* @} */
		/** \name Member functions (Get)  */
		/* @{ */

		/**
		 * Fill measured particles with data from tree. 
		 * @param index The event number. 
		 * @return An error code.
		 */
		virtual int Event(int KLFITTER_UNUSED(index))
		{ return 0; };  

		/**
		 * Return the number of events. 
		 * @return The number of events
		 */ 
		virtual int NEvents()
		{ return 0; }; 

		/**
		 * Return the measured missing transverse energy. 
		 * @return The missing ET
		 */ 
		virtual double ET_miss()
		{ return 0; }; 

		/**
		 * Return the measured missing transverse energy (x component). 
		 * @return The missing ET (x component)
		 */ 
		virtual double ET_miss_x()
		{ return 0; }; 

		/**
		 * Return the measured missing transverse energy (y component). 
		 * @return The missing ET (y component)
		 */ 
		virtual double ET_miss_y()
		{ return 0; }; 

		/* @} */
		/** \name Member functions (Set)  */
		/* @{ */
		/**
		 * Set a flag. If flag is true the input is Signal MC.
		 * Truth particle container is filled.
		 */ 
		void SetFlagIsSignalMC(bool flag)
		{ fFlagIsSignalMC = flag; }; 
		/* @} */
		/** \name Member functions (misc)  */
		/* @{ */

		/**
		 * Open Root file containing tree.
		 * @param filename The filename. 
		 * @param opt Options.
		 * @return An error code.
		 */ 
		virtual int OpenRootFile(const char * filename, Option_t * opt = "READ"); 

		/**
		 * Close Root file. 
		 * @return An error code.
		 */ 
		virtual int CloseRootFile(); 

		/* @} */

	protected: 
		/**
		 * A flag for using Signal MC as input.
		 */ 
		bool fFlagIsSignalMC;
		/**
		 * The Root file. 
		 */ 
		TFile * fRootFile; 

		/* @} */
	}; 

} // namespace KLFitter 

// --------------------------------------------------------- 

#endif // __INTERFACEROOT__H

