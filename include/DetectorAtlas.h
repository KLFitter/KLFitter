/*!
 * \class KLFitter::DetectorAtlas
 * \brief A class for describing of the ATLAS detector. 
 * \author Kevin Kr&ouml;ninger
 * \version 1.3
 * \date 03.12.2009
 *
 * This class holds the description of the ATLAS detector. 
 */

// --------------------------------------------------------- 

#ifndef DETECTORATLAS
#define DETECTORATLAS

// --------------------------------------------------------- 

#include "DetectorBase.h"

// --------------------------------------------------------- 

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{

  class DetectorAtlas : public DetectorBase
  {
                
  public: 
                
    /** \name Constructors and destructors */ 
    /* @{ */ 
                
    /** 
     * The default constructor. 
     */ 
    DetectorAtlas(); 
                
    /**
     * The default destructor.
     */
    virtual ~DetectorAtlas(); 

    /* @} */
    /** \name Member functions (Get)  */
    /* @{ */

    /**
     * Return the energy resolution of light jets.
     * @param eta The eta of the particle.
     * @return A pointer to the energy resolution object. 
     */ 
    KLFitter::ResolutionBase * ResEnergyLightJet(double eta = 0.); 

    /**
     * Return the energy resolution of b jets.
     * @param eta The eta of the particle.
     * @return A pointer to the energy resolution object. 
     */ 
    KLFitter::ResolutionBase * ResEnergyBJet(double eta = 0.); 

    /**
     * Return the energy resolution of gluon jets.
     * @param eta The eta of the particle.
     * @return A pointer to the energy resolution object. 
     */ 
    KLFitter::ResolutionBase * ResEnergyGluonJet(double eta = 0.); 

    /**
     * Return the energy resolution of electrons.
     * @param eta The eta of the particle.
     * @return A pointer to the energy resolution object. 
     */ 
    KLFitter::ResolutionBase * ResEnergyElectron(double eta = 0.); 

    /**
     * Return the energy resolution of muons.
     * @param eta The eta of the particle.
     * @return A pointer to the energy resolution object. 
     */ 
    KLFitter::ResolutionBase * ResEnergyMuon(double eta = 0.); 

    /**
     * Return the energy resolution of photons.
     * @param eta The eta of the particle.
     * @return A pointer to the energy resolution object. 
     */ 
    KLFitter::ResolutionBase * ResEnergyPhoton(double eta = 0.); 

    /**
     * Return the missing ET resolution. 
     * @return A pointer to the missing ET resolution.
     */ 
    KLFitter::ResolutionBase * ResMissingET(); 

    /**
     * Return the eta resolution of light jets.
     * @param eta The eta of the particle.
     * @return A pointer to the eta resolution object. 
     */ 
    KLFitter::ResolutionBase * ResEtaLightJet(double eta = 0.); 

    /**
     * Return the eta resolution of b jets.
     * @param eta The eta of the particle.
     * @return A pointer to the eta resolution object. 
     */ 
    KLFitter::ResolutionBase * ResEtaBJet(double eta = 0.);

    /**
     * Return the phi resolution of light jets.
     * @param eta The phi of the particle.
     * @return A pointer to the phi resolution object. 
     */ 
    KLFitter::ResolutionBase * ResPhiLightJet(double eta = 0.);

    /**
     * Return the phi resolution of b jets.
     * @param eta The phi of the particle.
     * @return A pointer to the phi resolution object. 
     */ 
    KLFitter::ResolutionBase * ResPhiBJet(double eta = 0.); 

    /* @} */

  private: 
    /**
     * The energy resolution of light jets for different eta regions. 
     */ 
    KLFitter::ResolutionBase * fResEnergyLightJet_eta1; 
    KLFitter::ResolutionBase * fResEnergyLightJet_eta2; 
    KLFitter::ResolutionBase * fResEnergyLightJet_eta3; 

    /**
     * The energy resolution of b jets for different eta regions. 
     */ 
    KLFitter::ResolutionBase * fResEnergyBJet_eta1; 
    KLFitter::ResolutionBase * fResEnergyBJet_eta2; 
    KLFitter::ResolutionBase * fResEnergyBJet_eta3; 

    /**
     * The energy resolution of gluon jets for different eta regions.
     */ 
    KLFitter::ResolutionBase * fResEnergyGluonJet_eta1; 
    KLFitter::ResolutionBase * fResEnergyGluonJet_eta2; 
    KLFitter::ResolutionBase * fResEnergyGluonJet_eta3; 

    /**
     * The energy resolution of electrons for different eta regions. 
     */ 
    KLFitter::ResolutionBase * fResEnergyElectron_eta1; 
    KLFitter::ResolutionBase * fResEnergyElectron_eta2; 
    KLFitter::ResolutionBase * fResEnergyElectron_eta3; 

    /**
     * The energy resolution of muons for different eta regions.
     */ 
    KLFitter::ResolutionBase * fResEnergyMuon_eta1; 
    KLFitter::ResolutionBase * fResEnergyMuon_eta2; 
    KLFitter::ResolutionBase * fResEnergyMuon_eta3; 

    /**
     * The energy resolution of photons for different eta regions.
     */ 
    KLFitter::ResolutionBase * fResEnergyPhoton_eta1; 
    KLFitter::ResolutionBase * fResEnergyPhoton_eta2; 
    KLFitter::ResolutionBase * fResEnergyPhoton_eta3; 

    /**
     * The eta resolution of light jets for different eta regions. 
     */ 
    KLFitter::ResolutionBase * fResEtaLightJet_eta1;
    KLFitter::ResolutionBase * fResEtaLightJet_eta2;
    KLFitter::ResolutionBase * fResEtaLightJet_eta3;

    /**
     * The eta resolution of b jets for different eta regions. 
     */ 
    KLFitter::ResolutionBase * fResEtaBJet_eta1;
    KLFitter::ResolutionBase * fResEtaBJet_eta2;
    KLFitter::ResolutionBase * fResEtaBJet_eta3;

    /**
     * The phi resolution of light jets for different eta regions. 
     */ 
    KLFitter::ResolutionBase * fResPhiLightJet_eta1;
    KLFitter::ResolutionBase * fResPhiLightJet_eta2;
    KLFitter::ResolutionBase * fResPhiLightJet_eta3;

    /**
     * The phi resolution of b jets for different eta regions. 
     */ 
    KLFitter::ResolutionBase * fResPhiBJet_eta1;
    KLFitter::ResolutionBase * fResPhiBJet_eta2;
    KLFitter::ResolutionBase * fResPhiBJet_eta3;

  }; 

} // namespace KLFitter 

// --------------------------------------------------------- 

#endif 
