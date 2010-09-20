/*!
 * \class KLFitter::Detector
 * \brief A class for describing a detector devided in different eta regions. 
 * \author Johannes Erdmann
 * \version 1.3
 * \date 03.12.2009
 *
 * This class holds the description of a detector with different eta regions. 
 */

// --------------------------------------------------------- 

#ifndef DETECTOR
#define DETECTOR

// --------------------------------------------------------- 

#include "DetectorBase.h"

#include <string>
#include <vector>

// --------------------------------------------------------- 

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{

  class Detector : public DetectorBase
  {
                
  public: 

    /** \name Constructors and destructors */ 
    /* @{ */ 
                
    /** 
     * The default constructor. 
     */ 
    Detector(std::vector<double> MaxValuesForEtaRegions,
             std::string eResTypeLightJet,
             std::string eResTypeBJet,
             std::string eResTypeGluonJet,
             std::string eResTypeElectron,
             std::string eResTypeMuon,
             std::string eResTypePhoton,
             std::vector<double> eResParLightJet,
             std::vector<double> eResParBJet,
             std::vector<double> eResParGluonJet,
             std::vector<double> eResParElectron,
             std::vector<double> eResParMuon,
             std::vector<double> eResParPhoton);

    /** 
     * Constructor with eta and phi resolutions for light and b-jets.
     */ 
    Detector(std::vector<double> MaxValuesForEtaRegions,
             std::string eResTypeLightJet,
             std::string eResTypeBJet,
             std::string eResTypeGluonJet,
             std::string eResTypeElectron,
             std::string eResTypeMuon,
             std::string eResTypePhoton,
             std::string etaResTypeLightJet,
             std::string etaResTypeBJet,
             std::string phiResTypeLightJet,
             std::string phiResTypeBJet,
             std::string metResType,
             std::vector<double> eResParLightJet,
             std::vector<double> eResParBJet,
             std::vector<double> eResParGluonJet,
             std::vector<double> eResParElectron,
             std::vector<double> eResParMuon,
             std::vector<double> eResParPhoton,
             std::vector<double> etaResParLightJet,
             std::vector<double> etaResParBJet,
             std::vector<double> phiResParLightJet,
             std::vector<double> phiResParBJet,
             std::vector<double> metResPar);
                
    /**
     * The default destructor.
     */
    virtual ~Detector(); 

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
     * vector that holds the maximum values for the several eta regions
     */
    std::vector<double> fMaxValuesForEtaRegions;


    /**
     * number of eta regions
     */
    unsigned int fNEtaRegions;

    /**
     * The energy resolution of light jets for different eta regions. 
     */
    std::vector<KLFitter::ResolutionBase*> * fResEnergyLightJet_etaRegions; 

    /**
     * The energy resolution of b jets for different eta regions. 
     */ 
    std::vector<KLFitter::ResolutionBase*> * fResEnergyBJet_etaRegions; 

    /**
     * The energy resolution of gluon jets for different eta regions.
     */ 
    std::vector<KLFitter::ResolutionBase*> * fResEnergyGluonJet_etaRegions; 

    /**
     * The energy resolution of electrons for different eta regions. 
     */ 
    std::vector<KLFitter::ResolutionBase*> * fResEnergyElectron_etaRegions; 

    /**
     * The energy resolution of muons for different eta regions.
     */ 
    std::vector<KLFitter::ResolutionBase*> * fResEnergyMuon_etaRegions; 

    /**
     * The energy resolution of photons for different eta regions.
     */ 
    std::vector<KLFitter::ResolutionBase*> * fResEnergyPhoton_etaRegions; 

    /**
     * The missing ET resolution. 
     */ 
    std::vector<KLFitter::ResolutionBase*> * fResMissingET_etaRegions; 

    /**
     * The eta resolution of light jets for different eta regions.
     */ 
    std::vector<KLFitter::ResolutionBase*> * fResEtaLightJet_etaRegions; 

    /**
     * The eta resolution of b jets for different eta regions.
     */ 
    std::vector<KLFitter::ResolutionBase*> * fResEtaBJet_etaRegions; 

    /**
     * The phi resolution of light jets for different eta regions.
     */ 
    std::vector<KLFitter::ResolutionBase*> * fResPhiLightJet_etaRegions; 

    /**
     * The phi resolution of b jets for different eta regions.
     */ 
    std::vector<KLFitter::ResolutionBase*> * fResPhiBJet_etaRegions; 

    /**
     * Get the vector of parameters that correspond to the eta region
     */
    std::vector<double> getDoubleGaussParameters(std::vector<double> parameters, int iEtaRegion);
  }; 

} // namespace KLFitter 

// --------------------------------------------------------- 

#endif 
