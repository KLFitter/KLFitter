/*!
 * \class KLFitter::ReadConfigFile
 * \brief A class for reading a config file for runKLFitter files
 * \author Sven Ebert
 * \author Johannes Erdmann
 * \version 1.3
 * \date 23.04.2010
 *
 * This class reads a config file with running settings.
 */

// --------------------------------------------------------- 

#ifndef READCONFIGFILE
#define READCONFIGFILE

// --------------------------------------------------------- 

#include <string>
#include "DetectorBase.h" 
#include "LikelihoodTopLeptonJets.h" 
#include "LikelihoodTopLeptonJetsUDSep.h" 

// --------------------------------------------------------- 

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{

  class ReadConfigFile
  {
                
  public: 
                
    /* @} */
    /** \name Constructors and destructors */ 
    /* @{ */ 
                
    /** 
     * The default constructor. 
     */ 
    ReadConfigFile(std::string filename); 
    //Second possible constructor.
    ReadConfigFile(std::string filename, bool * validconfig); 
                
    /**
     * The default destructor.
     */
    virtual ~ReadConfigFile(); 

    /* @} */
    /** \name Member functions (Get)  */
    /* @{ */

    /**
     * Return the LeptonType flag
     * @return The flag.
     */ 
    KLFitter::LikelihoodTopLeptonJets::LeptonType GetLeptonType() { return LeptonType; }

    /**
     * Return the DO_BATCH flag
     * @return The flag.
     */ 
    bool GetDO_BATCH() { return DO_BATCH; }

    /**
     * Return BTaggingMethod
     * @return The flag.
     */ 
    KLFitter::LikelihoodBase::BtaggingMethod GetBTaggingMethod() { return BTaggingMethod; }

    /**
     * Return FlagIntegrate
     * @return The flag.
     */ 
    bool GetFlagIntegrate() { return FlagIntegrate; }

    /**
     * Return FlagTopMassFixed
     * @return The flag.
     */ 
    bool GetFlagTopMassFixed() { return FlagTopMassFixed; }

    /**
     * Return FlagHiggsMassFixed
     * @return The flag.
     */ 
    bool GetFlagHiggsMassFixed() { return FlagHiggsMassFixed; }

    /**
     * Return FlagUseJetMass
     * @return The flag.
     */ 
    bool GetFlagUseJetMass() { return FlagUseJetMass; }

    /**
     * Return FlagWriteSignalMCTruth
     * @return The flag.
     */ 
    bool GetFlagWriteSignalMCTruth() { return FlagWriteSignalMCTruth; }

    /**
     * Return FlagIsSignalMC
     * @return The flag.
     */ 
    double GetTopMass() { return TopPoleMass; }
    
    /**
     * Return FlagIsSignalMC
     * @return The flag.
     */ 
    double GetHiggsMass() { return HiggsMass; }

    /**
     * Return Fthe BTagging cut value
     * @return The value.
     */ 
    double GetCutBTagging() { return CutBTagging; }

    /**
     * Return input_path
     * @return The path.
     */
    std::string GetInputPath() { return input_path; }

    /**
     * Return output_path
     * @return The path.
     */
    std::string GetOutputPath() { return output_path; }

    /**
     * Return IsBkg
     * @return The flag.
     */
    bool GetIsBkg() { return IsBkg; }

    /**
     * Return BeamCMEnergy
     * @return The flag.
     */
    KLFitter::DetectorBase::BeamCMEnergy GetBeamCMEnergy() { return BeamEnergy; }	


    /**
     * Flag for using the truth selection.
     */
    bool GetFlagTruthSel() { return FlagTruthSel; }

    /**
     * Flag for using the light quark permutation.
     */
    KLFitter::LikelihoodTopLeptonJetsUDSep::LJetSeparationMethod GetLJetSeparationMethod() { return LJetSeparationMethod; }

    /**
     * Flag for doing the comparison to TopKLFitter.
     */
    bool GetFlagAthenaComp() { return FlagAthenaComp; }

  private: 

    int ReadConfig(std::string filename);
    int GetTrueOrFalse(std::string line, size_t found);
    int GetValue(double * ret, std::string line, size_t found);
    int GetPath(std::string * ret, std::string line, size_t found);
    bool IsNumber(std::string::iterator a, int * number);
    bool CheckIOPath();

    /**
     * The LeptonType flag.
     */ 
    KLFitter::LikelihoodTopLeptonJets::LeptonType LeptonType;

    /**
     * The DO_BATCH flag.
     */ 
    bool DO_BATCH;

    /**
     * The BTaggingMethod.
     */ 
    KLFitter::LikelihoodBase::BtaggingMethod BTaggingMethod;

    /**
     * The FlagIntegrate.
     */ 
    bool FlagIntegrate;

    /**
     * The FlagTopMassFixed.
     */ 
    bool FlagTopMassFixed;

    /**
     * The FlagHiggsMassFixed.
     */ 
    bool FlagHiggsMassFixed;
    
    /**
     * The FlagUseJetMass.
     */ 
    bool FlagUseJetMass;

    /**
     * The FlagWriteSignalMCTruth.
     */ 
    bool FlagWriteSignalMCTruth;

    /**
     * The top mass.
     */ 
    double TopPoleMass;

    /**
     * The Higgs mass.
     */ 
    double HiggsMass;
    
    /**
     * The cut value.
     */ 
    double CutBTagging;
    /**
     * Path of the inputfile.
     */ 
    std::string input_path;

    /**
     * Path of the outputfile.
     */ 
    std::string output_path;

    /**
     * For background (w/o truth info!)
     */ 
    bool IsBkg;
    /**
     * Flag for using 7TeV transferfunctions
     */ 
    KLFitter::DetectorBase::BeamCMEnergy BeamEnergy;

    /**
     * Flag for using the truth selection.
     */
    bool FlagTruthSel;

    /**
     * Flag for using the light quark (and b Quark) reweighting.
     */
    KLFitter::LikelihoodTopLeptonJetsUDSep::LJetSeparationMethod LJetSeparationMethod;

    /**
     * Flag for doing the comparison to TopKLFitter.
     */
    bool FlagAthenaComp;

  }; 

} // namespace KLFitter 

// --------------------------------------------------------- 

#endif 
