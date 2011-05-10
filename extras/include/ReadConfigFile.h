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
     * Return the DO_ELEcTRON flag
     * @return The flag.
     */ 
    bool GetDO_ELECTRON() { return DO_ELECTRON; }

    /**
     * Return the DO_MUON flag
     * @return The flag.
     */ 
    bool GetDO_MUON() { return DO_MUON; }

    /**
     * Return the DO_BATCH flag
     * @return The flag.
     */ 
    bool GetDO_BATCH() { return DO_BATCH; }

    /**
     * Return FlagBTagging
     * @return The flag.
     */ 
    bool GetFlagBTagging() { return FlagBTagging; }

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
    double GetTopMass() { return MassTop; }

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
     * Return FlagIs7TeV
     * @return The flag.
     */
    bool GetFlagIs7TeV() { return FlagIs7TeV; }	

    /**
     * Return FlagIs10TeV
     * @return The flag.
     */
    bool GetFlagIs10TeV() { return FlagIs10TeV; }

    /**
     * Flag for using the truth selection.
     */
    bool GetFlagTruthSel() { return FlagTruthSel; }

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
     * The DO_ELECTRON flag.
     */ 
    bool DO_ELECTRON;

    /**
     * The DO_MUON flag.
     */ 
    bool DO_MUON;

    /**
     * The DO_BATCH flag.
     */ 
    bool DO_BATCH;

    /**
     * The FlagBTagging.
     */ 
    bool FlagBTagging;

    /**
     * The FlagIntegrate.
     */ 
    bool FlagIntegrate;

    /**
     * The FlagTopMassFixed.
     */ 
    bool FlagTopMassFixed;

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
    double MassTop;
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
    bool FlagIs7TeV;
    /**
     * Flag for using 10TeV transferfunctions
     */ 
    bool FlagIs10TeV;

    /**
     * Flag for using the truth selection.
     */
    bool FlagTruthSel;

    /**
     * Flag for doing the comparison to TopKLFitter.
     */
    bool FlagAthenaComp;

  }; 

} // namespace KLFitter 

// --------------------------------------------------------- 

#endif 
