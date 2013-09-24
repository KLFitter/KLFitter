/*!
 * \class KLFitter::PhysicsConstants
 * \brief A class containing physics constants. 
 * \author Kevin Kr&ouml;ninger
 * \version 1.3
 * \date 03.12.2009
 * \added ttH constants: Leonid Serkin
 * \date 01.12.2011
 *
 * This class contains physics constants. 
 */

// --------------------------------------------------------- 

#ifndef PHYSICSCONSTANTS
#define PHYSICSCONSTANTS

// --------------------------------------------------------- 

// --------------------------------------------------------- 

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{

  class PhysicsConstants
  {
                
  public: 
                
    /** \name Constructors and destructors */ 
    /* @{ */ 
                
    /** 
     * The default constructor. 
     */ 
    PhysicsConstants(); 
                
    /**
     * The default destructor.
     */
    virtual ~PhysicsConstants(); 

    /* @} */
    /** \name Member functions (Get)  */
    /* @{ */

    /**
     * Return the mass of the bottom quark in GeV/c2. 
     * @return The mass of the particle in GeV/c2. 
     */ 
    double MassBottom()
    { return fMassBottom; }; 

    /**
     * Return the msas of the W boson in GeV/c2.
     * @return The mass of the particle in GeV/c2. 
     */ 
    double MassW()
    { return fMassW; }; 

    /**
     * Return the msas of the top quark in GeV/c2
     * @return The mass of the particle in GeV/c2. 
     */ 
    double MassTop()
    { return fMassTop; }; 

    /**
     * Return the msas of the Higgs boson in GeV/c2
     * @return The mass of the particle in GeV/c2. 
     */ 
    double MassHiggs()
    { return fMassHiggs; }; 

    /**
     * Return the width of the W boson in GeV/c2. 
     * @return The width of the particle in GeV/c2. 
     */ 
    double GammaW()
    { return fGammaW; }; 

    /**
     * Return the width of the top quark in GeV/c2
     * @return The width of the particle in GeV/c2. 
     */ 
    double GammaTop()
    { return fGammaTop; }; 


    /**
     * Return the width of the Higgs boson in GeV/c2
     * @return The width of the particle in GeV/c2. 
     */ 
    double GammaHiggs()
    { return fGammaHiggs; }; 


    double MassTopUnc()
    { return fMassTopUnc; };


    /* @} */
    /** \name Member functions (Set)  */
    /* @{ */

    /**
     * Set the mass of the bottom quark in GeV/c2. 
     * @param mass The mass of the particle in GeV/c2. 
     * @return An error code. 
     */ 
    int SetMassBottom(double mass); 

    /**
     * Set the mass of the top quark in GeV/c2. 
     * @param mass The mass of the particle in GeV/c2. 
     * @return An error code. 
     */ 
    int SetMassTop(double mass); 

    /**
     * Set the mass of the Higgs boson in GeV/c2. 
     * @param mass The mass of the particle in GeV/c2. 
     * @return An error code. 
     */ 
    int SetMassHiggs(double mass); 

    /**
     * Set the mass of the W boson in GeV/c2. 
     * @param mass The mass of the particle in GeV/c2. 
     * @return An error code. 
     */ 
    int SetMassW(double mass); 

    /**
     * Set the width of the W boson in GeV/c2. 
     * @param gamma The width of the particle in GeV/c2. 
     * @return An error code. 
     */ 
    int SetGammaW(double gamma); 

    /**
     * Set the width of the top quark in GeV/c2.
     * @param gamma The width of the particle in GeV/c2. 
     * @return An error code. 
     */ 
    int SetGammaTop(double gamma); 


    /**
     * Set the width of the Higgs boson in GeV/c2.
     * @param gamma The width of the particle in GeV/c2. 
     * @return An error code. 
     */ 
    int SetGammaHiggs(double gamma); 

    /* @} */
    /** \name Member functions (misc)  */
    /* @{ */

    /**
     * Calculates the top width at NLO. 
     */ 
    void CalculateGammaTop(); 


    /**
     * Calculates the Higgs width using HDECAY. 
     */ 
    void CalculateGammaHiggs(); 

    /* @} */

  protected: 

  private: 

    /**
     * The bottom quark pole mass in GeV/c2.
     */ 
    double fMassBottom;

    /**
     * The W boson pole mass in GeV/c2.
     */ 
    double fMassW;

    /**
     * The top quark pole mass in GeV/c2.
     */ 
    double fMassTop;
 
    /**
     * The Higgs boson mass in GeV/c2.
     */ 
    double fMassHiggs;
    
    /**
     * The W boson width in GeV/c2.
     */ 
    double fGammaW;

    /**
     * The top quark width in GeV/c2.
     */ 
    double fGammaTop;

    /**
     * The Higgs boson width in GeV/c2.
     */ 
    double fGammaHiggs;

    /**
     * The Fermi constant. 
     */ 
    double fGF; 

    /**
     * alpha_S at m_Z
     */ 
    double fAlphaS; 

    /**
     * The top quark mass LHC uncertainty in GeV/c2.
     */ 
    double fMassTopUnc;

  }; 

} // namespace KLFitter 

// --------------------------------------------------------- 

#endif 

