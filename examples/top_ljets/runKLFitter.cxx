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
#include "InterfaceD3PD_ljets.h" 
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

bool EventTruthSelection(KLFitter::Particles * particles, KLFitter::LikelihoodTopLeptonJets::LeptonType Lepton); 
int GetCombinationNo(TString s);//int matchHadB, int matchLepB, int matchLJ1, int matchLJ2);
TString GetTrueCombinationString(std::vector<int> v0, std::vector<int> v1, std::vector<int> v2, std::vector<int> v3);


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
  KLFitter::LikelihoodTopLeptonJets::LeptonType Lepton = configReader->GetLeptonType();
  bool FlagIntegrate    = configReader->GetFlagIntegrate();
  bool FlagTopMassFixed = configReader->GetFlagTopMassFixed();
  bool FlagWriteSignalMCTruth   = configReader->GetFlagWriteSignalMCTruth();
  KLFitter::DetectorBase::BeamCMEnergy BeamEnergy = configReader->GetBeamCMEnergy();
  KLFitter::LikelihoodBase::BtaggingMethod Btagmethod = configReader->GetBTaggingMethod();
  double CutBTagging   = configReader->GetCutBTagging();
  double MassTop = configReader->GetTopMass();
  std::string input_file=configReader->GetInputPath();
  std::string output_file=configReader->GetOutputPath();
  bool FlagTruthSel = configReader->GetFlagTruthSel();

  delete configReader;
  if(!valid){ return 0;}//std::cout<<"Error: InputPath=OutputPath. Will not overwrite InputFile!"<<std::endl;return 0;}

  //counters
  unsigned int countLMatch = 0;
  unsigned int countJMatch = 0;

  // create new fitter 
  KLFitter::Fitter * myFitter = new KLFitter::Fitter(); 

	// open Root file 
  KLFitter::InterfaceRoot * myInterfaceRoot = new KLFitter::InterfaceD3PD_ljets();
  std::vector<std::string> inputfiles = myInterfaceRoot->ReadInputFiles(input_file.c_str());
  std::cout << "Input Files: " << std::endl;
  for(unsigned int i=0; i<inputfiles.size(); i++)
  	std::cout << inputfiles.at(i) << std::endl;
  myInterfaceRoot -> OpenRootFiles(inputfiles);
  myInterfaceRoot -> SetBtaggingInformation(CutBTagging, 0.696, 134.);

  // create detector
  KLFitter::DetectorBase * myDetector;
  if (BeamEnergy==KLFitter::DetectorBase::k7TeV)
    myDetector = new KLFitter::DetectorAtlas_7TeV("../../transferfunctions/7TeV/ttbar/mc11c"); 
  else if (BeamEnergy==KLFitter::DetectorBase::k10TeV)
    myDetector = new KLFitter::DetectorAtlas_10TeV("../../transferfunctions/10TeV/ttbar");
  else{std::cout<<"Error: Detector could not be created, please check the transferfunction flags"<<std::endl;return 1;}

  if (!myFitter -> SetDetector(myDetector))
    return 0; 

  // create likelihood for ttbar->e+jets channel 
  KLFitter::LikelihoodTopLeptonJets * myLikelihood = new KLFitter::LikelihoodTopLeptonJets(); 
  myLikelihood -> PhysicsConstants() -> SetMassTop(MassTop); 
  // b-tagging settings: (kNotag/kVeto/kWorkingPoint, TaggerCutValue, efficiency[0,1], rejection[>1])
  // Make sure to set btag rejection and efficiency if btagging set to a working
  myLikelihood -> SetBTagging(Btagmethod);
  myLikelihood -> SetFlagIntegrate(FlagIntegrate); 
  myLikelihood -> SetFlagTopMassFixed(FlagTopMassFixed);
  myInterfaceRoot -> WriteSignalMCTruth(FlagWriteSignalMCTruth, KLFitter::InterfaceRoot::kHerwig);

  if (Lepton==KLFitter::LikelihoodTopLeptonJets::kElectron)
    myLikelihood -> SetLeptonType(1); // set lepton type to electron 
  if (Lepton==KLFitter::LikelihoodTopLeptonJets::kMuon)
    myLikelihood -> SetLeptonType(2); // set lepton type to muon

  if (!myFitter -> SetLikelihood(myLikelihood))
    return 0; 

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
//  mySelectionTool -> RequireNBJets(0.601713, 2, 0); 


  mySelectionTool -> SetMaxNJetsForFit(4);
  if (Lepton==KLFitter::LikelihoodTopLeptonJets::kElectron)
    mySelectionTool -> RequireNElectronsPt(20.0, 1); 
  if (Lepton==KLFitter::LikelihoodTopLeptonJets::kMuon)
    mySelectionTool -> RequireNMuonsPt(20.0, 1);
  mySelectionTool -> RequireMET(20.); 

  // create matching tool
  KLFitter::MatchingTool * myMatchingTool = 0x0;
  if (FlagWriteSignalMCTruth)
    if (FlagTruthSel && FlagWriteSignalMCTruth)
      myMatchingTool = new KLFitter::MatchingTool( myFitter -> PParticles(), myInterfaceRoot -> PParticlesTruth() ); 

  // set fitter and truth particles 
  myInterfaceOutput -> SetFitter(myFitter); 
  if (FlagWriteSignalMCTruth)
    myInterfaceOutput -> SetParticlesTruth( myInterfaceRoot -> PParticlesTruth() ); 
  myInterfaceOutput -> SetParticlesMeasured( myInterfaceRoot -> PParticles() ); 
  if (myMatchingTool)
    myInterfaceOutput -> SetMatchingTool(myMatchingTool); 
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
      KLFitter::Particles * truthparticles = 0x0;
      if (FlagWriteSignalMCTruth)
        truthparticles = myInterfaceRoot -> ParticlesTruth(); 

      // truth event selection
      if (FlagWriteSignalMCTruth)
        if (FlagTruthSel)
          if (!EventTruthSelection(truthparticles, Lepton))
            continue; 

      // select event
      if (!mySelectionTool -> SelectEvent(measuredparticles, myInterfaceRoot -> ET_miss())) 
        continue; 

      // get particles from selection tool
      KLFitter::Particles * particles = mySelectionTool -> ParticlesSelected(); 

      // add particles to fitter 
      if (!myFitter -> SetParticles(particles))
        return 0;       
      // add ETmiss and sumET to fitter
      if (!myFitter -> SetET_miss_XY_SumET( myInterfaceRoot -> ET_miss_x(), myInterfaceRoot -> ET_miss_y(), myInterfaceRoot -> Sum_ET() ) )
        return 0;

      if (myMatchingTool)
        {
          // perform matching
          myMatchingTool -> MatchTruthAll(KLFitter::Particles::kParton); 
          if (Lepton==KLFitter::LikelihoodTopLeptonJets::kElectron)
            myMatchingTool -> MatchTruthAll(KLFitter::Particles::kElectron); 
          if (Lepton==KLFitter::LikelihoodTopLeptonJets::kMuon)
            myMatchingTool -> MatchTruthAll(KLFitter::Particles::kMuon);
                                        
          // identify true permutation
          std::vector<int> v0 = myMatchingTool->ListMatchedTruth(0, KLFitter::Particles::kParton);
          std::vector<int> v1 = myMatchingTool->ListMatchedTruth(1, KLFitter::Particles::kParton);
          std::vector<int> v2 = myMatchingTool->ListMatchedTruth(2, KLFitter::Particles::kParton);
          std::vector<int> v3 = myMatchingTool->ListMatchedTruth(3, KLFitter::Particles::kParton);
          //                      int trueCombi = GetCombinationNo(GetTrueCombinationString(v0, v1, v2, v3));
                                        
          //                      // require the event to be matched
          //                      if (trueCombi==-1)
          //                              continue;
        }

      // copy information into output class
      myInterfaceOutput -> SetEventWeight(weight); 
      myInterfaceOutput -> FillTreeMeasured(); 
      myInterfaceOutput -> FillTreeSelected(); 
      if (FlagWriteSignalMCTruth)
        {
          myInterfaceOutput -> FillTreeTruth();
          myInterfaceOutput -> FillTreeMatching();
        }
      myInterfaceOutput -> FillTreeMap();

      if (firstevent) {
        printf("----------------------------------------------------------------------------------------------\n");
				printf("--------------------------------FIT RESULTS FOR THE FIRST EVENT-------------------------------\n"); 
      }

      // loop over all permutations 
      for (int iperm = 0; iperm < myFitter -> Permutations() -> NPermutations(); ++iperm)
        {
          // fit the first permutation
          myFitter -> Fit(iperm); 
          // copy information into output class

          myInterfaceOutput -> FillTreeModelPermutation(); 

          if (firstevent) {
						KLFitter::Particles * OutputParticles = myFitter->Likelihood()->ParticlesModel(); 
  					KLFitter::Particles * InputParticles = *myFitter->Likelihood()->PParticlesPermuted();
						std::vector<double> Par = myFitter->Likelihood()->GetBestFitParameters();
            std::vector<double> ParErrors = myFitter->Likelihood()->GetBestFitParameterErrors();
						unsigned int ConvergenceStatusBitWord = myFitter->ConvergenceStatus();
						printf("----------------------------------------------------------------------------------------------\n");
						printf("----------------------------------------Permutation %2i----------------------------------------\n",iperm);
						printf("----------------------------------------------------------------------------------------------\n");
						printf("                  | hadronic b quark | leptonic b quark  |  light quark 1   |  light quark 2  |\n");
						printf("Input Energies    | %16.2f | %17.2f | %16.2f | %15.2f |\n",
						InputParticles->Parton(0)->E(), InputParticles->Parton(1)->E(),
						InputParticles->Parton(2)->E(), InputParticles->Parton(3)->E() );
						printf("Output Energies   | %16.2f | %17.2f | %16.2f | %15.2f |\n",
						OutputParticles->Parton(0)->E(), OutputParticles->Parton(1)->E(),
						OutputParticles->Parton(2)->E(), OutputParticles->Parton(3)->E() );
						printf("----------------------------------------------------------------------------------------------\n");
						printf("                  |  Log(Likelihood) | Event Probability |   Top Pole Mass  |   Neutrino pz   |\n");
						printf("Fitting Variables | %16.2f | %17.2E | %6.2f +- %6.2f | %6.2f +- %5.2f |\n",
						myFitter->Likelihood()->LogLikelihood(myFitter->Likelihood()->GetBestFitParameters()),
						exp(myFitter->Likelihood()->LogEventProbability()),
						Par[KLFitter::LikelihoodTopLeptonJets::parTopM],ParErrors[KLFitter::LikelihoodTopLeptonJets::parTopM],
						Par[KLFitter::LikelihoodTopLeptonJets::parNuPz],ParErrors[KLFitter::LikelihoodTopLeptonJets::parNuPz]);
						printf("----------------------------------------------------------------------------------------------\n");
						printf("                  | Minuit Not Conv. | Fit Aborted: NaN  | >=1 Par at Limit | Invalid TF@Conv.|\n");
						printf("Status Code       | %16i | %17i | %16i | %15i |\n",
						bool((ConvergenceStatusBitWord & KLFitter::Fitter::MinuitDidNotConvergeMask) != 0), 
						bool((ConvergenceStatusBitWord & KLFitter::Fitter::FitAbortedDueToNaNMask) != 0),
						bool((ConvergenceStatusBitWord & KLFitter::Fitter::AtLeastOneFitParameterAtItsLimitMask) != 0),
						bool((ConvergenceStatusBitWord & KLFitter::Fitter::InvalidTransferFunctionAtConvergenceMask) != 0));
				    }
        }

      // fill tree
      myInterfaceOutput -> FillTrees(); 

      // do not print following events
      if (firstevent) {
        printf("----------------------------------------------------------------------------------------------\n");
        firstevent = false;
      }

      // check if the particles have been matched
      if (myMatchingTool)
        {
          bool leptonsMatched = true;
          if (Lepton==KLFitter::LikelihoodTopLeptonJets::kElectron) {
            for (int iElec = 0; iElec < truthparticles -> NElectrons(); iElec++) {
              int nMatchedElectron = myMatchingTool -> NMatchedTruth(iElec, KLFitter::Particles::kElectron); 
              if (nMatchedElectron == 0)
                leptonsMatched = false;
            }
          }
          if (Lepton==KLFitter::LikelihoodTopLeptonJets::kMuon) {
            for (int iMu = 0; iMu < truthparticles -> NMuons(); iMu++) {
              int nMatchedMuon = myMatchingTool -> NMatchedTruth(iMu, KLFitter::Particles::kMuon); 
              if (nMatchedMuon == 0)
                leptonsMatched = false;
            }
          }

          if (!leptonsMatched)
            continue;

          countLMatch++;

          bool jetsMatched = true;
          std::set<TString> set_listOfMatches;
          set_listOfMatches.clear();

          for (int iQuark = 0; iQuark < truthparticles -> NPartons()-2; iQuark++) {
            int nMatchedJet = myMatchingTool -> NMatchedTruth(iQuark, KLFitter::Particles::kParton); 
            std::vector<int> listOfMatches = myMatchingTool -> ListMatchedTruth(iQuark, KLFitter::Particles::kParton);

            // check if quark is matched is to exactly one jet
            if (nMatchedJet != 1)
              jetsMatched = false;

            // check if another quark has already been matched to that jet
            TString string_listOfMatches = "";
            for (unsigned int j=0; j<listOfMatches.size(); j++)
              string_listOfMatches += listOfMatches.at(j);
                                
            if (set_listOfMatches.find(string_listOfMatches) != set_listOfMatches.end())
              jetsMatched = false;
            else
              set_listOfMatches.insert(string_listOfMatches);

          }

          set_listOfMatches.clear();
                        
          if (!jetsMatched)
            continue;
                        
          countJMatch++;
        }
    }
  
  // output cut flow 
  std::cout << " N (all)       : " << mySelectionTool -> CounterEvents() << std::endl;
  if (Lepton==KLFitter::LikelihoodTopLeptonJets::kElectron)
    std::cout << " N (electrons) : " << mySelectionTool -> CounterElectrons() << std::endl;
  if (Lepton==KLFitter::LikelihoodTopLeptonJets::kMuon)
    std::cout << " N (muons  )   : " << mySelectionTool -> CounterMuons() << std::endl;
  std::cout << " N (jets)      : " << mySelectionTool -> CounterJets() << std::endl;
  std::cout << " N (MET)       : " << mySelectionTool -> CounterMET() << std::endl;
  std::cout << " N (selected)  : " << mySelectionTool -> CounterSelected() << std::endl;
  if (FlagWriteSignalMCTruth)
    {
      if (Lepton==KLFitter::LikelihoodTopLeptonJets::kElectron)
        std::cout << " N (e  matched) : " << countLMatch << std::endl;
      if (Lepton==KLFitter::LikelihoodTopLeptonJets::kMuon)
        std::cout << " N (mu matched) : " << countLMatch << std::endl;
      std::cout << " N (j matched) : " << countJMatch << std::endl;
    }

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
  if (myMatchingTool)
    delete myMatchingTool;

  // no error 
  return 1; 

}

// -----------------------------------------------------------------------------------------------------------

bool EventTruthSelection(KLFitter::Particles * particles, KLFitter::LikelihoodTopLeptonJets::LeptonType Lepton)
{
  // ---------------------------------------------------------
  // truth selection
  // ---------------------------------------------------------

  int nTruthElectrons = particles -> NElectrons();
  int nTruthMuons = particles -> NMuons();
  int nTruthTaus = particles -> NTaus();
  int nTruthPartons = particles -> NPartons();

  if (Lepton==KLFitter::LikelihoodTopLeptonJets::kElectron) {
    // require exactly 1 truth electron
    if (nTruthElectrons != 1)
      return false;
                
    // require no truth muon or tau
    if (nTruthMuons != 0) return false;
    if (nTruthTaus != 0) return false;
  }

  if (Lepton==KLFitter::LikelihoodTopLeptonJets::kMuon) {
    // require exactly 1 truth muon
    if (nTruthMuons != 1)
      return false;
                
    // require no truth electron or tau
    if (nTruthElectrons != 0) return false;
    if (nTruthTaus != 0) return false;
  }

  // require exactly 6 truth quarks (including the 2 tops)
  if (nTruthPartons !=6) return false;

  return true;

}

///////////////////////////////////////////////////////////
int GetCombinationNo(TString s) { // int matchHadB, int matchLepB, int matchLJ1, int matchLJ2){

  // Match String to KLFitter Combinatorics
  if (s == "0 1 2 3") return  0;
  if (s == "0 2 1 3") return  1;
  if (s == "0 3 1 2") return  2;
  if (s == "3 0 1 2") return  3;
  if (s == "3 2 0 1") return  4;
  if (s == "3 1 0 2") return  5;
  if (s == "2 1 0 3") return  6;
  if (s == "2 3 0 1") return  7;
  if (s == "2 0 1 3") return  8;
  if (s == "1 0 2 3") return  9;
  if (s == "1 3 0 2") return 10;
  if (s == "1 2 0 3") return 11;
  // return permutation n+1000 if the two light quarks are reversed
  if (s == "0 1 3 2") return 1000;
  if (s == "0 2 3 1") return 1001;
  if (s == "0 3 2 1") return 1002;
  if (s == "3 0 2 1") return 1003;
  if (s == "3 2 1 0") return 1004;
  if (s == "3 1 2 0") return 1005;
  if (s == "2 1 3 0") return 1006;
  if (s == "2 3 1 0") return 1007;
  if (s == "2 0 3 1") return 1008;
  if (s == "1 0 3 2") return 1009;
  if (s == "1 3 2 0") return 1010;
  if (s == "1 2 3 0") return 1011;

  return -1;

}

///////////////////////////////////////////////////////////
TString GetTrueCombinationString(std::vector<int> v0, std::vector<int> v1, std::vector<int> v2, std::vector<int> v3) {
  // combination string
  TString s = "";
  for (unsigned int i = 0; i < v0.size(); i++)
    if (v0.at(i) == 1)
      s += i;
  s += " ";
  for (unsigned int i = 0; i < v1.size(); i++)
    if (v1.at(i) == 1)
      s += i;
  s += " ";
  for (unsigned int i = 0; i < v2.size(); i++)
    if (v2.at(i) == 1)
      s += i;
  s += " ";
  for (unsigned int i = 0; i < v3.size(); i++)
    if (v3.at(i) == 1)
      s += i;
  return s;
}
