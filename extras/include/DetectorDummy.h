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

#ifndef DETECTORDUMMY
#define DETECTORDUMMY

#include "DetectorBase.h"
#include <string>

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
     * @param folder The folder with transfer function parameters.
     */ 
    DetectorDummy(std::string folder = ""); 
                
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

#endif 
