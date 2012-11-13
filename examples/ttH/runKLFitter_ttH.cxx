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
#include "InterfaceD3PD_ttH.h" 
#include "Particles.h" 
#include "Permutations.h"
#include "InterfaceOutput.h" 
#include "LikelihoodBase.h" 
#include "LikelihoodTTHLeptonJets.h" 
#include "MatchingTool.h" 
#include "SelectionTool.h" 
#include "TString.h"
#include "TSystem.h"  
#include "ReadConfigFile.h"

#include "TFile.h"

bool EventTruthSelection(KLFitter::Particles * particles, KLFitter::LikelihoodTTHLeptonJets::LeptonType Lepton); 
int GetCombinationNo(TString s);
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
  KLFitter::ReadConfigFile * configReader = new KLFitter::ReadConfigFile("config_ttH.conf", &valid);
  KLFitter::LikelihoodTTHLeptonJets::LeptonType Lepton = (KLFitter::LikelihoodTTHLeptonJets::LeptonType) configReader->GetLeptonType();
 
  bool FlagIntegrate    = configReader->GetFlagIntegrate();
  bool FlagTopMassFixed = configReader->GetFlagTopMassFixed();
  bool FlagHiggsMassFixed = configReader->GetFlagHiggsMassFixed();
  bool FlagWriteSignalMCTruth   = configReader->GetFlagWriteSignalMCTruth();
  KLFitter::DetectorBase::BeamCMEnergy BeamEnergy = configReader->GetBeamCMEnergy();
  KLFitter::LikelihoodBase::BtaggingMethod Btagmethod = configReader->GetBTaggingMethod();
  double CutBTagging   = configReader->GetCutBTagging();
  double MassTop = configReader->GetTopMass();
  double MassHiggs = configReader->GetHiggsMass();
  std::string input_file=configReader->GetInputPath();
  std::string output_file=configReader->GetOutputPath();
  bool FlagTruthSel = configReader->GetFlagTruthSel();

  delete configReader;
  if(!valid){ return 0;}//std::cout<<"Error: InputPath=OutputPath. Will not overwrite InputFile!"<<std::endl;return 0;}

  //counters
  unsigned int countLMatch = 0;
  unsigned int countJMatch = 0;
  unsigned int count_bfromH = 0;
  unsigned int count_bbarfromH = 0;
  unsigned int count_0bsfromH = 0;
  unsigned int count_1bsfromH = 0;
  unsigned int count_2bsfromH = 0;
  unsigned int count_b_overlap = 0;
  unsigned int count_bbar_overlap = 0;
  unsigned int count_both_overlap = 0;
 
  // create new fitter 
  KLFitter::Fitter * myFitter = new KLFitter::Fitter(); 

  // open Root file 
  KLFitter::InterfaceRoot * myInterfaceRoot = new KLFitter::InterfaceD3PD_ttH();

  std::vector<std::string> inputfiles = myInterfaceRoot->ReadInputFiles(input_file.c_str());
  std::cout << "Input Files: " << std::endl;
  for(unsigned int i=0; i<inputfiles.size(); i++)
    std::cout << inputfiles.at(i) << std::endl;
  myInterfaceRoot -> OpenRootFiles(inputfiles);
  myInterfaceRoot -> SetBtaggingInformation(CutBTagging, 0.696, 134.);
  
  // create detector
  KLFitter::DetectorBase * myDetector;
  if (BeamEnergy==KLFitter::DetectorBase::k7TeV)
    myDetector = new KLFitter::DetectorAtlas_7TeV("../../transferfunctions/7TeV/ttbar/mc11c"); // mc11c transfer functions
  else if (BeamEnergy==KLFitter::DetectorBase::k10TeV)
    myDetector = new KLFitter::DetectorAtlas_10TeV("../../transferfunctions/10TeV/ttbar");
  else{std::cout<<"Error: Detector could not be created, please check the transferfunction flags"<<std::endl;return 1;}
  
  if (!myFitter -> SetDetector(myDetector))
    return 0; 

  // create likelihood for ttbarH->e+jets channel 
  KLFitter::LikelihoodTTHLeptonJets * myLikelihood = new KLFitter::LikelihoodTTHLeptonJets(); 
  myLikelihood -> PhysicsConstants() -> SetMassTop(MassTop); 
  myLikelihood -> PhysicsConstants() -> SetMassHiggs(MassHiggs); 
 
  // Make sure to set btag rejection and efficiency if btagging set to a working
  myLikelihood -> SetBTagging(Btagmethod);  // b-tagging settings: (kNotag/kVeto/kWorkingPoint, TaggerCutValue, efficiency[0,1], rejection[>1])
  
  myLikelihood -> SetFlagIntegrate(FlagIntegrate); 
  myLikelihood -> SetFlagTopMassFixed(FlagTopMassFixed);
  myLikelihood -> SetFlagHiggsMassFixed(FlagHiggsMassFixed);

  myInterfaceRoot -> WriteSignalMCTruth(FlagWriteSignalMCTruth, KLFitter::InterfaceRoot::kAcer); // hack to pick Pythia mathing

  if (Lepton==KLFitter::LikelihoodTTHLeptonJets::kElectron)
    myLikelihood -> SetLeptonType(1); // set lepton type to electron 
  if (Lepton==KLFitter::LikelihoodTTHLeptonJets::kMuon)
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
  mySelectionTool -> RequireNJetsPt(25.0, 6, -1); 
 
  // Require NBjets
  mySelectionTool -> RequireNBJets(0.601713, 4, 0);
	
  // lepton selection
  if (Lepton==KLFitter::LikelihoodTTHLeptonJets::kElectron) {
    mySelectionTool -> RequireNElectronsPt(25.0, 1); 
    mySelectionTool -> RequireMET(30.); 
    mySelectionTool -> RequireMWT(30.); 
  }
  if (Lepton==KLFitter::LikelihoodTTHLeptonJets::kMuon) {
    mySelectionTool -> RequireNMuonsPt(20.0, 1);
    mySelectionTool -> RequireTriangular(60.); 
  }
  
  // NJets to fit
  mySelectionTool -> SetMaxNJetsForFit(6);

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
      if ((ievent-minEv+1)%100 == 0)
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
      if (!mySelectionTool -> SelectEvent(measuredparticles, myInterfaceRoot -> ET_miss(), myInterfaceRoot->MWT()))
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
          if (Lepton==KLFitter::LikelihoodTTHLeptonJets::kElectron)
            myMatchingTool -> MatchTruthAll(KLFitter::Particles::kElectron); 
          if (Lepton==KLFitter::LikelihoodTTHLeptonJets::kMuon)
            myMatchingTool -> MatchTruthAll(KLFitter::Particles::kMuon);
                                        
          // identify true permutation
          std::vector<int> v0 = myMatchingTool->ListMatchedTruth(0, KLFitter::Particles::kParton); // b from H
          std::vector<int> v1 = myMatchingTool->ListMatchedTruth(1, KLFitter::Particles::kParton); // bbar from H
          std::vector<int> v2 = myMatchingTool->ListMatchedTruth(2, KLFitter::Particles::kParton); // hadronic b quark
          std::vector<int> v3 = myMatchingTool->ListMatchedTruth(3, KLFitter::Particles::kParton); // leptonic b quark
          std::vector<int> v4 = myMatchingTool->ListMatchedTruth(4, KLFitter::Particles::kParton); // light quark 1
          std::vector<int> v5 = myMatchingTool->ListMatchedTruth(5, KLFitter::Particles::kParton); // light quark 2
	  
	  // std::cout << "True combination: " << GetTrueCombinationString(v0, v1, v2, v3, v4, v5) << std::endl;
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
	    printf("                  | hadronic b quark | leptonic b quark  |  light quark 1   |  light quark 2  |    Higgs b 1   |  Higgs b 2  |    \n");
	    printf("Input Energies    | %16.2f | %17.2f | %16.2f | %15.2f | %15.2f | %15.2f |   \n",
		   InputParticles->Parton(0)->E(), InputParticles->Parton(1)->E(), InputParticles->Parton(2)->E(), 
		   InputParticles->Parton(3)->E(), InputParticles->Parton(4)->E(), InputParticles->Parton(5)->E()  );

	    printf("Output Energies   | %16.2f | %17.2f | %16.2f | %15.2f | %15.2f | %15.2f |    \n",
		   OutputParticles->Parton(0)->E(), OutputParticles->Parton(1)->E(), OutputParticles->Parton(2)->E(), 
		   OutputParticles->Parton(3)->E(), OutputParticles->Parton(4)->E(), OutputParticles->Parton(5)->E() );

	    printf("----------------------------------------------------------------------------------------------\n");
	    printf("                  |  Log(Likelihood) | Event Probability |   Top Pole Mass  |   Neutrino pz   |   Higgs Mass  |  \n");
	    printf("Fitting Variables | %16.2f | %17.2E | %6.2f +- %6.2f | %6.2f +- %5.2f | %6.2f +- %5.2f | \n",
		   myFitter->Likelihood()->LogLikelihood(myFitter->Likelihood()->GetBestFitParameters()),
		   exp(myFitter->Likelihood()->LogEventProbability()),
		   Par[KLFitter::LikelihoodTTHLeptonJets::parTopM],ParErrors[KLFitter::LikelihoodTTHLeptonJets::parTopM],
		   Par[KLFitter::LikelihoodTTHLeptonJets::parNuPz],ParErrors[KLFitter::LikelihoodTTHLeptonJets::parNuPz],
		   Par[KLFitter::LikelihoodTTHLeptonJets::parHiggsM],ParErrors[KLFitter::LikelihoodTTHLeptonJets::parHiggsM] );


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
	  // check if b and bbar from H overlap
	  std::vector<int> BfromH_list = myMatchingTool -> ListMatchedTruth(0, KLFitter::Particles::kParton);
	  std::vector<int> BbarfromH_list = myMatchingTool -> ListMatchedTruth(1, KLFitter::Particles::kParton);
	  
	  int ntruth = truthparticles->NPartons();
	  int nmeasured = particles->NPartons();
	  int nmatches_jets_bfromH = 0;
	  int nmatches_partons_bfromH = 0;
	  int nmatches_jets_bbarfromH = 0;					
	  int nmatches_partons_bbarfromH = 0;					
	  int index_bfromH = -1;
	  int index_bbarfromH = -1;
	  
	  for (int i = 0; i < nmeasured; ++i) {
	    int temp_b = myMatchingTool->MatchVectors(i, 0, KLFitter::Particles::kParton);
	    int temp_bbar = myMatchingTool->MatchVectors(i, 1, KLFitter::Particles::kParton);
	    nmatches_jets_bfromH += temp_b;
	    nmatches_jets_bbarfromH += temp_bbar;
	    if (temp_b > 0)
	      index_bfromH = i;
	    if (temp_bbar > 0)
	      index_bbarfromH = i;
	  }
	  
	  if (nmatches_jets_bfromH > 0) {
	    for (int j = 0; j < ntruth; ++j) {
	      int temp = myMatchingTool->MatchVectors(index_bfromH, j, KLFitter::Particles::kParton);
	      nmatches_partons_bfromH+=temp;
	    }
	  }
	  
	  if (nmatches_jets_bbarfromH > 0) {
	    for (int j = 0; j < ntruth; ++j) {
	      int temp = myMatchingTool->MatchVectors(index_bbarfromH, j, KLFitter::Particles::kParton);
	      nmatches_partons_bbarfromH+=temp;
	    }
	  }
	  
	  int nMatchedBfromH = myMatchingTool->NMatchedTruth(0, KLFitter::Particles::kParton);
	  int nMatchedBbarfromH = myMatchingTool->NMatchedTruth(1, KLFitter::Particles::kParton);
	  
	  // veto if overlap
	  if (nmatches_jets_bfromH == 1 && nmatches_partons_bfromH == 1)
	    nMatchedBfromH = 1;
	  else
	    nMatchedBfromH = 0;
	  
	  if (nmatches_jets_bbarfromH == 1 && nmatches_partons_bbarfromH == 1)
	    nMatchedBbarfromH = 1;
	  else
	    nMatchedBbarfromH = 0;
	  
	  if ( (nmatches_jets_bfromH > 1) || (nmatches_partons_bfromH > 1))
	    count_b_overlap++;
	  if ( (nmatches_jets_bbarfromH > 1) || (nmatches_partons_bbarfromH > 1))
	    count_bbar_overlap++;
	  if ( ((nmatches_jets_bfromH > 1) || (nmatches_partons_bfromH > 1))
	       && ((nmatches_jets_bbarfromH > 1) || (nmatches_partons_bbarfromH > 1)))
	    count_both_overlap++;
	  
	  
	  // increase counter
	  if (nMatchedBfromH > 0)
	    count_bfromH++;
	  if (nMatchedBbarfromH > 0)
	    count_bbarfromH++;
	  
	  if (nMatchedBfromH == 0 && nMatchedBbarfromH == 0)
	    count_0bsfromH++;
	  else if ((nMatchedBfromH > 0 && nMatchedBbarfromH == 0) || (nMatchedBfromH == 0 && nMatchedBbarfromH > 0))
	    count_1bsfromH++;
	  else if (nMatchedBfromH > 0 && nMatchedBbarfromH > 0)
	    count_2bsfromH++;
	  
          bool leptonsMatched = false;
          if (Lepton==KLFitter::LikelihoodTTHLeptonJets::kElectron) {
            for (int iElec = 0; iElec < truthparticles -> NElectrons(); iElec++) {
              int nMatchedElectron = myMatchingTool -> NMatchedTruth(iElec, KLFitter::Particles::kElectron); 
              if (nMatchedElectron > 0)
                leptonsMatched = true;
            }
          }
	  
          if (Lepton==KLFitter::LikelihoodTTHLeptonJets::kMuon) {
            for (int iMu = 0; iMu < truthparticles -> NMuons(); iMu++) {
              int nMatchedMuon = myMatchingTool -> NMatchedTruth(iMu, KLFitter::Particles::kMuon); 
              if (nMatchedMuon > 0)
                leptonsMatched = true;
            }
          }
	  
	  //          if (!leptonsMatched)
	  //            continue;
	  
          countLMatch++;
	  
          bool jetsMatched = true;
          std::set<TString> set_listOfMatches;
          set_listOfMatches.clear();

          for (int iQuark = 2; iQuark < truthparticles -> NPartons()-3; iQuark++) {
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
  std::cout << " Selection: " << std::endl;
  std::cout << " ---------- " << std::endl;
  std::cout << " N (after preselection)       : " << mySelectionTool -> CounterEvents() << std::endl;
  if (Lepton==KLFitter::LikelihoodTTHLeptonJets::kElectron)
    std::cout << " N (after electron cuts)      : " << mySelectionTool -> CounterElectrons() << std::endl;
  if (Lepton==KLFitter::LikelihoodTTHLeptonJets::kMuon)
    std::cout << " N (after muon cuts )         : " << mySelectionTool -> CounterMuons() << std::endl;
  std::cout << " N (after jet cuts)           : " << mySelectionTool -> CounterJets() << std::endl;
  std::cout << " N (after MET cuts)           : " << mySelectionTool -> CounterMET() << std::endl;
  std::cout << " N (after MWT cuts)           : " << mySelectionTool -> CounterMWT() << std::endl;
  std::cout << " N (after triangular cuts)    : " << mySelectionTool -> CounterTriangular() << std::endl;
  std::cout << " N (after b-jet cuts)         : " << mySelectionTool -> CounterBJets() << std::endl;
  std::cout << " N (after all cuts)           : " << mySelectionTool -> CounterSelected() << std::endl;
  
  if (FlagWriteSignalMCTruth)
    {
      std::cout << std::endl;
      std::cout << " Matching: " << std::endl;
      std::cout << " --------- " << std::endl;
      std::cout << " N (selected)                  : " << mySelectionTool -> CounterSelected() << std::endl;
      if (Lepton==KLFitter::LikelihoodTTHLeptonJets::kElectron)
        std::cout << " N (with matched e )           : " << countLMatch << std::endl;
      if (Lepton==KLFitter::LikelihoodTTHLeptonJets::kMuon)
        std::cout << " N (with matched muon)         : " << countLMatch << std::endl;
      std::cout << " N (with b from H matched)     : " << count_bfromH << std::endl;
      std::cout << " N (with bbar from H matched)  : " << count_bbarfromH << std::endl;
      std::cout << " N (with 0 b's from H matched) : " << count_0bsfromH << std::endl;
      std::cout << " N (with 1 b's from H matched) : " << count_1bsfromH << std::endl;
      std::cout << " N (with 2 b's from H matched) : " << count_2bsfromH << std::endl;
      std::cout << " N (with b overlap)            : " << count_b_overlap << std::endl;
      std::cout << " N (with bbar overlap)         : " << count_bbar_overlap << std::endl;
      std::cout << " N (with both overlap)         : " << count_both_overlap << std::endl;
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

bool EventTruthSelection(KLFitter::Particles * particles, KLFitter::LikelihoodTTHLeptonJets::LeptonType Lepton)
{
  // ---------------------------------------------------------
  // truth selection
  // ---------------------------------------------------------

  int nTruthElectrons = particles -> NElectrons();
  int nTruthMuons = particles -> NMuons();
  int nTruthTaus = particles -> NTaus();
  //int nTruthPartons = particles -> NPartons();

  if (Lepton==KLFitter::LikelihoodTTHLeptonJets::kElectron) {
    // require exactly 1 truth electron
    if (nTruthElectrons != 1)
      return false;
                
    // require no truth muon or tau
    if (nTruthMuons != 0) return false;
    if (nTruthTaus != 0) return false;
  }

  if (Lepton==KLFitter::LikelihoodTTHLeptonJets::kMuon) {
    // require exactly 1 truth muon
    if (nTruthMuons != 1)
      return false;
                
    // require no truth electron or tau
    if (nTruthElectrons != 0) return false;
    if (nTruthTaus != 0) return false;
  }

  // require exactly 8 truth quarks (including the 2 tops)
  //  if (nTruthPartons !=6) return false;

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
TString GetTrueCombinationString(std::vector<int> v0, std::vector<int> v1, std::vector<int> v2, std::vector<int> v3, std::vector<int> v4, std::vector<int> v5) {
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
  s += " ";
  for (unsigned int i = 0; i < v4.size(); i++)
    if (v4.at(i) == 1)
      s += i;
  s += " ";
  for (unsigned int i = 0; i < v5.size(); i++)
    if (v5.at(i) == 1)
      s += i;
  return s;
}

