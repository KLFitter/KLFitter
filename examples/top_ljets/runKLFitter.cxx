#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include "Fitter.h" 
#include "PhysicsConstants.h" 
#include "DetectorAtlas_7TeV.h" 
#include "DetectorAtlas_10TeV.h" 
#include "DetectorDummy.h" 
#include "InterfaceRoot.h" 
#include "InterfaceDummy.h" 
#include "InterfaceD3PD.h" 
#include "Particles.h" 
#include "Permutations.h"
#include "InterfaceOutput.h" 
#include "LikelihoodBase.h" 
#include "LikelihoodTopLeptonJets.h" 
#include "MatchingTool.h" 
#include "SelectionTool.h" 
#include "TString.h"
#include "TSystem.h"  
#include "ReadConfigFile.h"

#include "TFile.h"

int main(int argc, char **argv) 
{
  bool valid=false;
  // parameters (1st: 'first event', 2nd: 'last event - 1'
  if (argc!=3) { // i.e. 2 parameters ...
    std::cout << "number of parameters is unequal 2" << std::endl;
    return 1;
  }

  // read configurating
  //      KLFitter::ReadConfigFile * configReader = new KLFitter::ReadConfigFile("examples/top_ljets/config.conf");
  KLFitter::ReadConfigFile * configReader = new KLFitter::ReadConfigFile("config.conf", &valid);
  double DO_ELECTRON = configReader->GetDO_ELECTRON();
  double DO_MUON     = configReader->GetDO_MUON();
  double DO_BATCH    = configReader->GetDO_BATCH();
  bool FlagIntegrate    = configReader->GetFlagIntegrate();
  bool FlagTopMassFixed = configReader->GetFlagTopMassFixed();
  bool FlagUseJetMass   = configReader->GetFlagUseJetMass();
  bool FlagIsSignalMC   = configReader->GetFlagIsSignalMC();
  bool FlagIs7TeV   = configReader->GetFlagIs7TeV();
  bool FlagIs10TeV   = configReader->GetFlagIs10TeV();
  double CutBTagging   = configReader->GetCutBTagging();

  double MassTop = configReader->GetTopMass();
  std::string input_file=configReader->GetInputPath();
  std::string output_file=configReader->GetOutputPath();
  
  delete configReader;
  if(!valid){ return 0;}//std::cout<<"Error: InputPath=OutputPath. Will not overwrite InputFile!"<<std::endl;return 0;}

  
  // create new fitter 
  KLFitter::Fitter * myFitter = new KLFitter::Fitter(); 

  // open Root file 
  KLFitter::InterfaceRoot * myInterfaceRoot = new KLFitter::InterfaceD3PD();
  std::vector<std::string> inputfiles = myInterfaceRoot->ReadInputFiles(input_file.c_str());
  std::cout << "Input Files: " << std::endl;
  for(unsigned int i=0; i<inputfiles.size(); i++)
  	std::cout << inputfiles.at(i) << std::endl;
  myInterfaceRoot -> OpenRootFiles(inputfiles);
	
  // create detector
  KLFitter::DetectorBase * myDetector;
  if (FlagIs7TeV && !FlagIs10TeV)
    myDetector = new KLFitter::DetectorAtlas_7TeV("../../transferfunctions/ttbar"); 
  else if (!FlagIs7TeV && FlagIs10TeV)
    myDetector = new KLFitter::DetectorAtlas_10TeV("../../transferfunctions/ttbar");
  else{std::cout<<"Error: Detector could not be created, please check the transferfunction flags"<<std::endl;return 1;}

  if (!myFitter -> SetDetector(myDetector))
    return 0; 
        
  // create likelihood for ttbar->e+jets channel 
  KLFitter::LikelihoodTopLeptonJets * myLikelihood = new KLFitter::LikelihoodTopLeptonJets(); 

  myLikelihood -> PhysicsConstants() -> SetMassTop(MassTop); 
  // b-tagging settings: kNotag / kVeto / kWorkingPoint
  myLikelihood -> SetBTagging(KLFitter::LikelihoodBase::kNotag); 
  // Make sure to set btag and efficiency if btagghing set to a working point!
  myLikelihood -> SetbtagEff(0.6); // between 0 and 1
  myLikelihood -> SetbtagRej(900.); // hopefully greater than 1
  myLikelihood -> SetFlagIntegrate(FlagIntegrate); 
  myLikelihood -> SetFlagTopMassFixed(FlagTopMassFixed);
  myLikelihood -> SetFlagUseJetMass(FlagUseJetMass);
  myLikelihood -> SetCutBTag(CutBTagging);
  myInterfaceRoot -> SetFlagIsSignalMC(FlagIsSignalMC);
        
  if (DO_ELECTRON)
    myLikelihood -> SetLeptonType(1); // set lepton type to electron 
  if (DO_MUON)
    myLikelihood -> SetLeptonType(2); // set lepton type to muon

  if (!myFitter -> SetLikelihood(myLikelihood))
    return 0; 

  if (DO_BATCH)
    gSystem->Exec("cp /work/pcatlas07/erdmann/e+jets/input.root .");

  // create interface for output 
  KLFitter::InterfaceOutput * myInterfaceOutput = new KLFitter::InterfaceOutput(); 

  // open Root file 
  if (!myInterfaceOutput -> OpenRootFile(output_file.c_str()))
    return 0; 

  // create selection tool
  KLFitter::SelectionTool * mySelectionTool = new KLFitter::SelectionTool(); 
  mySelectionTool -> SelectElectronEta(2.5);
  mySelectionTool -> SelectMuonEta(2.5);
  mySelectionTool -> SelectPhotonEta(2.5);
  mySelectionTool -> SelectJetEta(2.5);
  mySelectionTool -> RequireNJetsPt(25.0, 4, -1); 
  //mySelectionTool -> RequireNJetsPt(40.0, 3, -1); 
  mySelectionTool -> SetMaxNJetsForFit(4);
  if (DO_ELECTRON)
    mySelectionTool -> RequireNElectronsPt(20.0, 1); 
  if (DO_MUON)
    mySelectionTool -> RequireNMuonsPt(20.0, 1);
  mySelectionTool -> RequireMET(20.); 

  // set fitter and particles 
  myInterfaceOutput -> SetFitter(myFitter); 
  myInterfaceOutput -> SetParticlesMeasured( myInterfaceRoot -> PParticles() ); 
  myInterfaceOutput -> SetSelectionTool(mySelectionTool); 

  // switch off filling of histograms
  //myFitter -> Likelihood() -> MCMCSetFlagFillHistograms(false);

  int nevents = myInterfaceRoot -> NEvents(); 
  int minEv = TString(argv[1]).Atoi();
  int maxEv = TString(argv[2]).Atoi();
  if (maxEv>nevents) {
    maxEv = nevents;
    std::cout << "parameter warning: last event parameter reset to maximum event number available (" << nevents << ")" << std::endl;
  }
  if (minEv>=maxEv) {
    std::cout << "parameter error: first event larger than last event  -->  exiting" << std::endl;
    return 1;
  }

  // for print-out for first fitted event
  bool firstevent(true);

  // loop over events
  for (int ievent = minEv; ievent < maxEv; ++ievent)
    {
      if ((ievent-minEv+1)%1000 == 0)
        std::cout << " event " << (ievent+1) << std::endl; 


      // get first event
      if (!myInterfaceRoot -> Event(ievent))
        return 0; 

      // get the event weight
      double weight = myInterfaceRoot->Weight();

      // read single event from root file and get particles 
      KLFitter::Particles * measuredparticles = myInterfaceRoot -> Particles(); 
      
      // select event
      if (!mySelectionTool -> SelectEvent(measuredparticles, myInterfaceRoot -> ET_miss())) 
        continue; 

      // get particles from selection tool
      KLFitter::Particles * particles = mySelectionTool -> ParticlesSelected(); 

      // add particles to fitter 
      if (!myFitter -> SetParticles(particles))
        return 0;       

      // add ETmiss to fitter
      if (!myFitter -> SetET_miss_XY( myInterfaceRoot -> ET_miss_x(), myInterfaceRoot -> ET_miss_y() ) )
        return 0;

      // copy information into output class
      myInterfaceOutput -> SetEventWeight(weight); 
      myInterfaceOutput -> FillTreeMeasured(); 
      myInterfaceOutput -> FillTreeSelected(); 
      
      myInterfaceOutput -> FillTreeMap(); 

      if (firstevent) {
        std::cout << "----------------------------------------------------------------------" << std::endl;
        std::cout << "FIT RESULTS FOR THE FIRST EVENT" << std::endl;
      }

      // loop over all permutations 
      for (int iperm = 0; iperm < myFitter -> Permutations() -> NPermutations(); ++iperm)
        {
          // fit the first permutation
          myFitter -> Fit(iperm); 
          // copy information into output class
          myInterfaceOutput -> FillTreeModelPermutation();

          if (firstevent) {
            std::cout << "----------" << std::endl;
            std::cout << "PERMUTATION NO. " << iperm << std::endl;
            std::cout << "log(likelihood) ................................. : " << myFitter->Likelihood()->LogLikelihood(myFitter->Likelihood()->GetBestFitParameters()) << std::endl;
            std::cout << "event probability (not normalized!).............. : " << exp(myFitter->Likelihood()->LogEventProbability()) << std::endl;
            unsigned int ConvergenceStatusBitWord = myFitter->ConvergenceStatus();
            std::cout << "convergence: MinuitDidNotConverge ............... : " << bool((ConvergenceStatusBitWord & KLFitter::Fitter::MinuitDidNotConvergeMask) != 0) << std::endl;
            std::cout << "convergence: FitAbortedDueToNaN ................. : " << bool((ConvergenceStatusBitWord & KLFitter::Fitter::FitAbortedDueToNaNMask) != 0) << std::endl;
            std::cout << "convergence: AtLeastOneFitParameterAtItsLimit ... : " << bool((ConvergenceStatusBitWord & KLFitter::Fitter::AtLeastOneFitParameterAtItsLimitMask) != 0) << std::endl;
            std::cout << "convergence: InvalidTransferFunctionAtConvergence : " << bool((ConvergenceStatusBitWord & KLFitter::Fitter::InvalidTransferFunctionAtConvergenceMask) != 0) << std::endl;
            std::vector<double> Parameters = myFitter->Likelihood()->GetBestFitParameters();
            std::vector<double> ParameterErrors = myFitter->Likelihood()->GetBestFitParameterErrors();
            std::cout << "fit parameter : central top mass value .......... : " << Parameters[KLFitter::LikelihoodTopLeptonJets::parTopM] << " +- " << ParameterErrors[KLFitter::LikelihoodTopLeptonJets::parTopM] << std::endl;
            std::cout << "fit parameter : neutrino pz ..................... : " << Parameters[KLFitter::LikelihoodTopLeptonJets::parNuPz] << " +- " << ParameterErrors[KLFitter::LikelihoodTopLeptonJets::parNuPz]<< std::endl;
          }
           
        }

      // fill tree
      myInterfaceOutput -> FillTrees(); 

      // do not print following events
      if (firstevent) {
        std::cout << "----------------------------------------------------------------------" << std::endl;
        firstevent = false;
      }
  }
      

  // output cut flow 
  std::cout << " N (all)       : " << mySelectionTool -> CounterEvents() << std::endl;
  if (DO_ELECTRON)
    std::cout << " N (electrons) : " << mySelectionTool -> CounterElectrons() << std::endl;
  if (DO_MUON)
    std::cout << " N (muons  )   : " << mySelectionTool -> CounterMuons() << std::endl;
  std::cout << " N (jets)      : " << mySelectionTool -> CounterJets() << std::endl;
  std::cout << " N (MET)       : " << mySelectionTool -> CounterMET() << std::endl;
  
  // close input file 
  if (!myInterfaceRoot -> CloseRootFile()) 
    return 0;

  // close output file 
  if (!myInterfaceOutput -> CloseRootFile()) 
    return 0; 

  // free memory 

  delete myInterfaceOutput; 
  delete myInterfaceRoot; 
  delete myDetector; 
  delete myFitter; 
  delete myLikelihood; 
  delete mySelectionTool;
  
  // no error 
  return 1; 

}
// -----------------------------------------------------------------------------------------------------------
