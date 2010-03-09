/*!
 * \class KLFitter::DetectorDummy
 * \brief A class for describing a dummy detector. 
 * \author Kevin Kr&ouml;ninger
 * \version 1.3
 * \date 03.12.2009
 *
 * This class holds the description of a dummy detector. 
 */

// --------------------------------------------------------- 

#ifndef __DETECTORDUMMY__H
#define __DETECTORDUMMY__H

#include "DetectorBase.h"
#include "ResGauss.h" 

// --------------------------------------------------------- 

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{

	class DetectorDummy : public DetectorBase
	{
		
	public: 
		
		/** \name Constructors and destructors */ 
		/* @{ */ 
		
		/** 
		 * The default constructor. 
		 */ 
		DetectorDummy(); 
		
		/**
		 * The default destructor.
		 */
		virtual ~DetectorDummy(); 

		/* @} */
		/** \name Member functions (Set)  */
		/* @{ */

		/* @} */

 private: 

	}; 

} // namespace KLFitter 

// --------------------------------------------------------- 

#endif // __DETECTORDUMMY__H
