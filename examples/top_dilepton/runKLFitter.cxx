#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <exception>
#include "Fitter.h" 
#include "PhysicsConstants.h" 
#include "DetectorAtlas_7TeV.h" 
#include "DetectorAtlas_10TeV.h" 
#include "DetectorDummy.h" 
#include "InterfaceRoot.h" 
#include "InterfaceDummy.h" 
#include "InterfaceD3PD_dilepton.h" 
#include "Particles.h" 
#include "Permutations.h"
#include "InterfaceOutput_dilepton.h" 
#include "LikelihoodBase.h" 
#include "LikelihoodTopDilepton.h" 
#include "MatchingTool.h" 
#include "SelectionTool.h" 
#include "TString.h"
#include "TSystem.h"  
#include "readparameters.h"

#include "BAT/BCModel.h"
#include "BAT/BCH1D.h"
#include "BAT/BCParameter.h"

#include "TFile.h"

bool EventTruthSelection(KLFitter::Particles * particles, KLFitter::LikelihoodTopDilepton::LeptonType Lepton1, KLFitter::LikelihoodTopDilepton::LeptonType Lepton2); 
int GetCombinationNo(TString s);//int matchHadB, int matchLepB, int matchLJ1, int matchLJ2);
TString GetTrueCombinationString(std::vector<int> v0, std::vector<int> v1, std::vector<int> v2, std::vector<int> v3);


int main(int argc, char **argv) 
{
  char* paramfile = argv[1];
  // parameters (1st: config_file, 2nd: 'first event', 3rd: 'last event - 1'
  if (argc!=4) { // i.e. 3 parameters ...
    std::cout << "number of parameters is unequal 3. Usage: ./runKLFitter.exe [config] [first-event] [last-event]" << std::endl;
    return 1;
  }

  ///// read configurating

  //auxiliar parameters
  std::string tmp_Lepton1, tmp_Lepton2, tmp_Btagmethod;
  int tmp_BeamEnergy;

  // parameters to read
  KLFitter::LikelihoodTopDilepton::LeptonType Lepton1(KLFitter::LikelihoodTopDilepton::kElectron);
  KLFitter::LikelihoodTopDilepton::LeptonType Lepton2(KLFitter::LikelihoodTopDilepton::kElectron);
  bool FlagIntegrate;
  bool FlagTopMassFixed;
  bool FlagWriteSignalMCTruth;
  KLFitter::DetectorBase::BeamCMEnergy BeamEnergy;
  KLFitter::LikelihoodBase::BtaggingMethod Btagmethod;
  double CutBTagging;
  double MassTop;
  std::string input_file;
  std::string output_file;
  bool FlagTruthSel;
  bool FlagSumloglik;

  try {
    readparameters rp(paramfile);
  
    tmp_Lepton1 = rp.get<std::string>("LeptonType1");
    tmp_Lepton2 = rp.get<std::string>("LeptonType2");
    FlagIntegrate    = rp.get<bool>("FlagIntegrate");
    FlagTopMassFixed = rp.get<bool>("FlagTopMassFixed");
    FlagWriteSignalMCTruth   = rp.get<bool>("FlagWriteSignalMCTruth");
    tmp_BeamEnergy =  rp.get<int>("BeamCMEnergy");
    tmp_Btagmethod = rp.get<std::string>("BTaggingMethod");
    CutBTagging =  rp.get<double>("CutBTagging");
    MassTop = rp.get<double>("TopPoleMass");
    input_file = rp.get<std::string>("PathToInputFile");
    output_file = rp.get<std::string>("PathToOutputFile");
    FlagTruthSel = rp.get<bool>("FlagTruthSel");
    FlagSumloglik = rp.get<bool>("FlagSumloglik");
  }
  catch (std::exception& e) {	    
    std::cerr<<e.what()<<std::endl; 
    return 0;
  } 

  // LeptonTypes 1 & 2
  if (tmp_Lepton1 == "electron")
    Lepton1 = KLFitter::LikelihoodTopDilepton::kElectron;
  else if (tmp_Lepton1 == "muon")
    Lepton1 = KLFitter::LikelihoodTopDilepton::kMuon;

  if (tmp_Lepton2 == "electron")
    Lepton2 = KLFitter::LikelihoodTopDilepton::kElectron;
  else if (tmp_Lepton2 == "muon")
    Lepton2 = KLFitter::LikelihoodTopDilepton::kMuon;

  // BeamEnergy
  if(tmp_BeamEnergy==7)
    BeamEnergy=KLFitter::DetectorBase::k7TeV;
  if(tmp_BeamEnergy==10)
    BeamEnergy=KLFitter::DetectorBase::k10TeV;

  // Btagmethod
  if(tmp_Btagmethod == "Notag")
    Btagmethod=KLFitter::LikelihoodBase::kNotag;
  if(tmp_Btagmethod == "Veto")
    Btagmethod=KLFitter::LikelihoodBase::kVeto;
  if(tmp_Btagmethod == "VetoNoFit")
    Btagmethod=KLFitter::LikelihoodBase::kVetoNoFit;		  
  if(tmp_Btagmethod == "WorkingPoint")
    Btagmethod=KLFitter::LikelihoodBase::kWorkingPoint;


  //counters
  unsigned int countLMatch = 0;
  unsigned int countJMatch = 0;

  // create new fitter 
  KLFitter::Fitter * myFitter = new KLFitter::Fitter(); 

  // set minimization method to MCMC
  myFitter->SetMinimizationMethod(KLFitter::Fitter::kSimulatedAnnealing);

  // open Root file 
  KLFitter::InterfaceRoot * myInterfaceRoot = new KLFitter::InterfaceD3PD_dilepton();
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


  ///////////////////////////////////////////////////////////
  // sigma_nueta dependency on mtop (a + b*mtop) => push_back (a,b)

  std::vector<double> m_nu_params_ee;
  std::vector<double> m_nu_params_emu;
  std::vector<double> m_nu_params_mumu;

  m_nu_params_ee.push_back(1.207);
  m_nu_params_ee.push_back(-2.4e-04);
  
  m_nu_params_emu.push_back(1.438);
  m_nu_params_emu.push_back(-8.55e-04);
  
  m_nu_params_mumu.push_back(1.295);
  m_nu_params_mumu.push_back(-4.00e-04);
  ///////////////////////////////////////////////////////////

  // create likelihood for ttbar->e+jets channel 
  KLFitter::LikelihoodTopDilepton * myLikelihood = new KLFitter::LikelihoodTopDilepton(); 
  myLikelihood -> PhysicsConstants() -> SetMassTop(MassTop); 
  // b-tagging settings: (kNotag/kVeto/kWorkingPoint, TaggerCutValue, efficiency[0,1], rejection[>1])
  // Make sure to set btag rejection and efficiency if btagging set to a working
  myLikelihood -> SetBTagging(Btagmethod);
  myLikelihood -> SetFlagIntegrate(FlagIntegrate); 
  myLikelihood -> SetFlagTopMassFixed(FlagTopMassFixed);
  myInterfaceRoot -> WriteSignalMCTruth(FlagWriteSignalMCTruth, KLFitter::InterfaceRoot::kHerwig);

  if (Lepton1==KLFitter::LikelihoodTopDilepton::kElectron && Lepton2==KLFitter::LikelihoodTopDilepton::kElectron){
    std::cout << "in the ee channel" << std::endl;
    myLikelihood -> SetLeptonType(1,1); // ee channel
    myLikelihood -> SetEtaNuParams(m_nu_params_ee);
  }
  else if (Lepton1==KLFitter::LikelihoodTopDilepton::kElectron && Lepton2==KLFitter::LikelihoodTopDilepton::kMuon){
    std::cout << "in the emu channel" << std::endl;
    myLikelihood -> SetLeptonType(1,2); // emu channel
    myLikelihood -> SetEtaNuParams(m_nu_params_emu);
  }
  else if (Lepton1==KLFitter::LikelihoodTopDilepton::kMuon && Lepton2==KLFitter::LikelihoodTopDilepton::kMuon){
    std::cout << "in the mumu channel" << std::endl;
    myLikelihood -> SetLeptonType(2,2); // mumu channel
    myLikelihood -> SetEtaNuParams(m_nu_params_mumu);
  }
  else
    std::cout << "Error: SetLeptonType1/2 wrong! check lepton types in config..." << std::endl;

  //if flag true, set "sumloglikelihood" option, instead of the default "bestpermutation" one
  myLikelihood -> SetDoSumLogLik(FlagSumloglik);

  if (!myFitter -> SetLikelihood(myLikelihood))
    return 0; 

  // create interface for output 
  KLFitter::InterfaceOutput_dilepton * myInterfaceOutput = new KLFitter::InterfaceOutput_dilepton(); 

  // open Root file 
  if (!myInterfaceOutput -> OpenRootFile(output_file.c_str()))
    return 0; 

  // create selection tool
  KLFitter::SelectionTool * mySelectionTool = new KLFitter::SelectionTool(); 
  mySelectionTool -> SelectElectronEta(2.5);
  mySelectionTool -> SelectMuonEta(2.5);
  mySelectionTool -> SelectPhotonEta(2.5);
  mySelectionTool -> SelectJetEta(2.5);
  mySelectionTool -> RequireNJetsPt(25.0, 2, -1);   // dilepton case
//  mySelectionTool -> RequireNBJets(0.601713, 2, 0);  // optional


  mySelectionTool -> SetMaxNJetsForFit(2);                            // dilepton case (can be set to >2)
  if (Lepton1==KLFitter::LikelihoodTopDilepton::kElectron && Lepton2==KLFitter::LikelihoodTopDilepton::kElectron){
    mySelectionTool -> RequireNElectronsPt(20.0, 2);
  }
  else if (Lepton1==KLFitter::LikelihoodTopDilepton::kElectron && Lepton2==KLFitter::LikelihoodTopDilepton::kMuon){
    mySelectionTool -> RequireNElectronsPt(20.0, 1); 
    mySelectionTool -> RequireNMuonsPt(20.0, 1);
  }
  else if (Lepton1==KLFitter::LikelihoodTopDilepton::kMuon && Lepton2==KLFitter::LikelihoodTopDilepton::kMuon){
    mySelectionTool -> RequireNMuonsPt(20.0, 2);
  }

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
  int minEv = TString(argv[2]).Atoi();
  int maxEv = TString(argv[3]).Atoi();
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
          if (!EventTruthSelection(truthparticles, Lepton1, Lepton2))
            continue; 

      // select event
      if (!mySelectionTool -> SelectEvent(measuredparticles, myInterfaceRoot -> ET_miss())){
	//std::cout << "*** event did not pass selection -> check! ***" << std::endl;
        continue; 
      }

      // get particles from selection tool
      KLFitter::Particles * particles = mySelectionTool -> ParticlesSelected(); 

      // add particles to fitter 
      if (!myFitter -> SetParticles(particles))
        return 0;  
      // add truth particles to fitter 
      if (FlagWriteSignalMCTruth)
	if (!myFitter -> SetMyParticlesTruth(truthparticles))
        return 0;  
      // add ETmiss and sumET to fitter
      if (!myFitter -> SetET_miss_XY_SumET( myInterfaceRoot -> ET_miss_x(), myInterfaceRoot -> ET_miss_y(), myInterfaceRoot -> Sum_ET() ) )
        return 0;

      if (myMatchingTool)
        {
          // perform matching
          myMatchingTool -> MatchTruthAll(KLFitter::Particles::kParton); 
	  if (Lepton1==KLFitter::LikelihoodTopDilepton::kElectron && Lepton2==KLFitter::LikelihoodTopDilepton::kElectron){
	    myMatchingTool -> MatchTruthAll(KLFitter::Particles::kElectron); 
	  }
	  else if (Lepton1==KLFitter::LikelihoodTopDilepton::kElectron && Lepton2==KLFitter::LikelihoodTopDilepton::kMuon){
	    myMatchingTool -> MatchTruthAll(KLFitter::Particles::kElectron); 
	    myMatchingTool -> MatchTruthAll(KLFitter::Particles::kMuon);
	  }
	  else if (Lepton1==KLFitter::LikelihoodTopDilepton::kMuon && Lepton2==KLFitter::LikelihoodTopDilepton::kMuon){
	   myMatchingTool -> MatchTruthAll(KLFitter::Particles::kMuon);
	  }
	  
                                  
          // identify true permutation
          std::vector<int> v0 = myMatchingTool->ListMatchedTruth(0, KLFitter::Particles::kParton);
          std::vector<int> v1 = myMatchingTool->ListMatchedTruth(1, KLFitter::Particles::kParton);
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
						double inputE_lep1(0.);
						double inputE_lep2(0.);
						double outputE_lep1(0.);
						double outputE_lep2(0.);
						if (Lepton1==KLFitter::LikelihoodTopDilepton::kElectron && Lepton2==KLFitter::LikelihoodTopDilepton::kElectron){
						  inputE_lep1 = InputParticles->Electron(0)->E();
						  inputE_lep2 = InputParticles->Electron(1)->E();
						  outputE_lep1 = OutputParticles->Electron(0)->E();
						  outputE_lep2 = OutputParticles->Electron(1)->E();
						}
						else if (Lepton1==KLFitter::LikelihoodTopDilepton::kElectron && Lepton2==KLFitter::LikelihoodTopDilepton::kMuon){
						  inputE_lep1 = InputParticles->Electron(0)->E();
						  inputE_lep2 = InputParticles->Muon(0)->E();
						  outputE_lep1 = OutputParticles->Electron(0)->E();
						  outputE_lep2 = OutputParticles->Muon(0)->E();
						}
						else if (Lepton1==KLFitter::LikelihoodTopDilepton::kMuon && Lepton2==KLFitter::LikelihoodTopDilepton::kMuon){
						  inputE_lep1 = InputParticles->Muon(0)->E();
						  inputE_lep2 = InputParticles->Muon(1)->E();
						  outputE_lep1 = OutputParticles->Muon(0)->E();
						  outputE_lep2 = OutputParticles->Muon(1)->E();
						}

						std::vector<double> Par = myFitter->Likelihood()->GetBestFitParameters();
						printf("----------------------------------------------------------------------------------------------\n");
						printf("----------------------------------------Permutation %2i----------------------------------------\n",iperm);
						printf("----------------------------------------------------------------------------------------------\n");
						printf("                  | b1 quark | b2 quark  | lepton 1  | lepton 2  |\n");
						printf("Input Energies    | %16.2f   | %17.2f    | %16.2f    | %15.2f |\n",
						       InputParticles->Parton(0)->E(), InputParticles->Parton(1)->E(),
						       inputE_lep1, inputE_lep2);
						printf("Output Energies   | %16.2f   | %17.2f    | %16.2f    | %15.2f |\n",
						       OutputParticles->Parton(0)->E(), OutputParticles->Parton(1)->E(),
						       outputE_lep1, outputE_lep2);
						printf("----------------------------------------------------------------------------------------------\n");
						printf("                  |  Log(Likelihood) | Event Probability  |\n");
						printf("Fitting Variables | %16.2f | %17.2E |\n",
						       myFitter->Likelihood()->LogLikelihood(myFitter->Likelihood()->GetBestFitParameters()),
						       exp(myFitter->Likelihood()->LogEventProbability()) );
						printf("----------------------------------------------------------------------------------------------\n");
						printf("                  |  Top Pole Mass  |  Nu Eta  |  AntiNu Eta  |\n");
						printf("Fitting Variables | %16.2f | %16.2f | %16.2f |\n",
						       Par[KLFitter::LikelihoodTopDilepton::parTopM],Par[KLFitter::LikelihoodTopDilepton::parNuEta],Par[KLFitter::LikelihoodTopDilepton::parAntiNuEta]);
						printf("----------------------------------------------------------------------------------------------\n");
						 
						std::vector<double> LHCompVec = myFitter->Likelihood()->LogLikelihoodComponents( myFitter->Likelihood()->GetBestFitParameters() );
						
						std::cout << "LHComp(0)= " << LHCompVec.at(0) << " LHComp(1)= " << LHCompVec.at(1)  << " LHComp(2)= " << LHCompVec.at(2)  <<
						  " LHComp(3)= " <<  LHCompVec.at(3) << " LHComp(4)= " <<  LHCompVec.at(4) << std::endl;
						std::cout << "LHComp(5)= " << LHCompVec.at(5) << " LHComp(6)= " << LHCompVec.at(6) <<" LHComp(7)= " << LHCompVec.at(7) << std::endl;
						
	  }

	  // // example: get marginalized histogram wrt Par(0)==mtop if using kMarkovChainMC!!
// 	  BCParameter * a = myFitter->Likelihood()->GetParameter(0);
	
// 	  if(myFitter->Likelihood()->GetMarginalized(a)){
// 	    if(myFitter->Likelihood()->GetMarginalized(a)->GetHistogram()->Integral() > 0){
// 	      TH1D* h_mtop_marginalized = (TH1D *)myFitter->Likelihood()->GetMarginalized(a)->GetHistogram();
// 	    }
// 	  }//get marginalized
	  
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

	  if (Lepton1==KLFitter::LikelihoodTopDilepton::kElectron && Lepton2==KLFitter::LikelihoodTopDilepton::kElectron){
	    for (int iElec = 0; iElec < truthparticles -> NElectrons(); iElec++) {
              int nMatchedElectron = myMatchingTool -> NMatchedTruth(iElec, KLFitter::Particles::kElectron); 
              if (nMatchedElectron == 0)
                leptonsMatched = false;
            }
	  }
	  else if (Lepton1==KLFitter::LikelihoodTopDilepton::kElectron && Lepton2==KLFitter::LikelihoodTopDilepton::kMuon){
	    for (int iElec = 0; iElec < truthparticles -> NElectrons(); iElec++) {
              int nMatchedElectron = myMatchingTool -> NMatchedTruth(iElec, KLFitter::Particles::kElectron); 
              if (nMatchedElectron == 0)
                leptonsMatched = false;
            }

	    for (int iMu = 0; iMu < truthparticles -> NMuons(); iMu++) {
              int nMatchedMuon = myMatchingTool -> NMatchedTruth(iMu, KLFitter::Particles::kMuon); 
              if (nMatchedMuon == 0)
                leptonsMatched = false;
            }
	  }
	  else if (Lepton1==KLFitter::LikelihoodTopDilepton::kMuon && Lepton2==KLFitter::LikelihoodTopDilepton::kMuon){
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

          for (int iQuark = 0; iQuark < truthparticles -> NPartons()-2; iQuark++) {  // considering tops in Partons
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
  if (Lepton1==KLFitter::LikelihoodTopDilepton::kElectron && Lepton2==KLFitter::LikelihoodTopDilepton::kElectron){
    std::cout << " N (electrons) : " << mySelectionTool -> CounterElectrons() << std::endl;
  }
  else if (Lepton1==KLFitter::LikelihoodTopDilepton::kElectron && Lepton2==KLFitter::LikelihoodTopDilepton::kMuon){
    std::cout << " N (electrons) : " << mySelectionTool -> CounterElectrons() << std::endl;
    std::cout << " N (muons)   : " << mySelectionTool -> CounterMuons() << std::endl;
  }
  else if (Lepton1==KLFitter::LikelihoodTopDilepton::kMuon && Lepton2==KLFitter::LikelihoodTopDilepton::kMuon){
    std::cout << " N (muons)   : " << mySelectionTool -> CounterMuons() << std::endl;
  }
  std::cout << " N (jets)      : " << mySelectionTool -> CounterJets() << std::endl;
  std::cout << " N (MET)       : " << mySelectionTool -> CounterMET() << std::endl;
  std::cout << " N (selected)  : " << mySelectionTool -> CounterSelected() << std::endl;                  
  if (FlagWriteSignalMCTruth)
    {
      if (Lepton1==KLFitter::LikelihoodTopDilepton::kElectron && Lepton2==KLFitter::LikelihoodTopDilepton::kElectron){
	std::cout << " N (e  matched) : " << countLMatch << std::endl;
      }
      else if (Lepton1==KLFitter::LikelihoodTopDilepton::kElectron && Lepton2==KLFitter::LikelihoodTopDilepton::kMuon){
	std::cout << " N (e&mu  matched) : " << countLMatch << std::endl;
      }
      else if (Lepton1==KLFitter::LikelihoodTopDilepton::kMuon && Lepton2==KLFitter::LikelihoodTopDilepton::kMuon){
	std::cout << " N (mu matched) : " << countLMatch << std::endl;
      }
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

bool EventTruthSelection(KLFitter::Particles * particles, KLFitter::LikelihoodTopDilepton::LeptonType Lepton1, KLFitter::LikelihoodTopDilepton::LeptonType Lepton2)
{
  // ---------------------------------------------------------
  // truth selection
  // ---------------------------------------------------------

  int nTruthElectrons = particles -> NElectrons();
  int nTruthMuons = particles -> NMuons();
  int nTruthTaus = particles -> NTaus();
  int nTruthPartons = particles -> NPartons();

  if (Lepton1==KLFitter::LikelihoodTopDilepton::kElectron && Lepton2==KLFitter::LikelihoodTopDilepton::kElectron){
    
    // require exactly 2 truth electron
    if (nTruthElectrons != 2) return false;
                
    // require no truth muon or tau
    if (nTruthMuons != 0) return false;
    if (nTruthTaus != 0) return false;
  }
  else if (Lepton1==KLFitter::LikelihoodTopDilepton::kElectron && Lepton2==KLFitter::LikelihoodTopDilepton::kMuon){
    // require exactly 1 truth electron
    if (nTruthElectrons != 1) return false;
    // require exactly 1 truth muon
    if (nTruthMuons != 1) return false;

    // require no truth tau
    if (nTruthTaus != 0) return false;
  }
  else if (Lepton1==KLFitter::LikelihoodTopDilepton::kMuon && Lepton2==KLFitter::LikelihoodTopDilepton::kMuon){
    // require exactly 2 truth muon
    if (nTruthMuons != 2) return false;
                
    // require no truth electron or tau
    if (nTruthElectrons != 0) return false;
    if (nTruthTaus != 0) return false;
  }

  // require exactly 4 truth quarks (including the 2 tops)
  if (nTruthPartons !=4) return false;

  return true;

}

///////////////////////////////////////////////////////////
int GetCombinationNo(TString s) { // int matchB1, int matchB2

  // Match String to KLFitter Combinatorics
  if (s == "0 1") return  0;
  if (s == "1 0") return  1;
 
  return -1;

}
///////////////////////////////////////////////////////////
TString GetTrueCombinationString(std::vector<int> v0, std::vector<int> v1) {
  // combination string
  TString s = "";
  for (unsigned int i = 0; i < v0.size(); i++)
    if (v0.at(i) == 1)
      s += i;
  s += " ";
  for (unsigned int i = 0; i < v1.size(); i++)
    if (v1.at(i) == 1)
      s += i;
 
  return s;
}
///////////////////////////////////////////////////////////
