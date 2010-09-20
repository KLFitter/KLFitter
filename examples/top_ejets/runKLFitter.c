#include <iostream> 
#include <vector>
#include <set>
#include <Fitter.h> 
#include <PhysicsConstants.h> 
#include <DetectorDummy.h> 
#include <InterfaceRoot.h> 
#include <InterfaceDummy.h> 
#include <Particles.h> 
#include <InterfaceOutput.h> 
#include <LikelihoodTopLeptonJets.h> 
#include <SelectionTool.h> 
#include <TString.h>
#include <TSystem.h>  

int main(int argc, char **argv) 
{
  // parameters (1st: 'first event', 2nd: 'last event - 1'
  if (argc!=3) { // i.e. 2 parameters ...
    std::cout << "Usage: ./runKLFitter.exe nmin nmax" << std::endl << std::endl;
    std::cout << "Run KLFitter with events with event number between" << std::endl;
    std::cout << "nmin and nmax." << std::endl << std::endl; 
    return 1;
  }

  // create new fitter 
  KLFitter::Fitter * myFitter = new KLFitter::Fitter(); 

  // open Root file 
  KLFitter::InterfaceRoot * myInterfaceRoot = new KLFitter::InterfaceDummy(); 
  myInterfaceRoot -> OpenRootFile("input.root");

  // create detector 
  KLFitter::DetectorBase * myDetector = new KLFitter::DetectorDummy(); 

  // set detector 
  if (!myFitter -> SetDetector(myDetector))
    return 0; 
        
  // create likelihood for ttbar->e+jets channel 
  KLFitter::LikelihoodTopLeptonJets * myLikelihood = new KLFitter::LikelihoodTopLeptonJets(); 
  myLikelihood -> PhysicsConstants() -> SetMassTop(172.5);   // set top mass to 172.5 GeV/c2
  myLikelihood -> SetFlagBTagging(false);                    // do not use b-tagging
  myLikelihood -> SetFlagIntegrate(false);                   // do not integrate Likelihood
  myLikelihood -> SetFlagTopMassFixed(false);                // top mass is a free parameter
  myLikelihood -> SetFlagUseJetMass(false);                  // do not use jet mass but PDG values 
  myLikelihood -> SetLeptonType(1);                          // set lepton type to electron

  // set likelihood
  if (!myFitter -> SetLikelihood(myLikelihood))
    return 0; 

  // create interface for output 
  KLFitter::InterfaceOutput * myInterfaceOutput = new KLFitter::InterfaceOutput(); 

  // open Root file 
  if (!myInterfaceOutput -> OpenRootFile("output.root"))
    return 0; 

  // create selection tool
  KLFitter::SelectionTool * mySelectionTool = new KLFitter::SelectionTool(); 
  mySelectionTool -> SelectElectronEta(2.5);
  mySelectionTool -> SelectMuonEta(2.5);
  mySelectionTool -> SelectPhotonEta(2.5);
  mySelectionTool -> SelectJetEta(2.5);
  mySelectionTool -> RequireNJetsPt(20.0, 4, -1); 
  mySelectionTool -> RequireNJetsPt(40.0, 3, -1); 
  mySelectionTool -> SetMaxNJetsForFit(4);
  mySelectionTool -> RequireNElectronsPt(20.0, 1);
  mySelectionTool -> RequireMET(20.); 

  // set fitter and truth particles 
  myInterfaceOutput -> SetFitter(myFitter); 
  myInterfaceOutput -> SetParticlesTruth( myInterfaceRoot -> PParticlesTruth() ); 
  myInterfaceOutput -> SetParticlesMeasured( myInterfaceRoot -> PParticles() ); 
  myInterfaceOutput -> SetSelectionTool(mySelectionTool); 

  int nevents = myInterfaceRoot -> NEvents(); 
  int minEv = TString(argv[1]).Atoi();
  int maxEv = TString(argv[2]).Atoi();
  if (maxEv>=nevents) {
    maxEv = nevents;
    std::cout << "Warning: nmax is larger than number of available events." << std::endl;
    std::cout << "Reset nmax to last event in the tree." << std::endl;
  }
  if (minEv > maxEv) {
    std::cout << "Error: nmin greater than nmax. Exit." << std::endl; 
  }
  if (minEv>=maxEv) {
    std::cout << "Error: nmax is greater than number of available events. Exit." << std::endl;
    return 1;
  }

  // loop over events
  for (int ievent = minEv; ievent < maxEv; ++ievent)
    {
      if ((ievent-minEv+1)%1000 == 0)
        std::cout << " event " << (ievent+1) << std::endl; 

      // get first event
      if (!myInterfaceRoot -> Event(ievent))
        return 0; 

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
      myInterfaceOutput -> FillTreeMeasured(); 
      myInterfaceOutput -> FillTreeSelected(); 
      myInterfaceOutput -> FillTreeTruth();
      myInterfaceOutput -> FillTreeMap(); 

      // loop over all permutations 
      for (int iperm = 0; iperm < myFitter -> Permutations() -> NPermutations(); ++iperm)
        {
          // fit the first permutation
          myFitter -> Fit(iperm); 
                                        
          // copy information into output class
          myInterfaceOutput -> FillTreeModelPermutation(); 
        }
      // fill tree
      myInterfaceOutput -> FillTrees(); 
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

  // no error 
  return 1; 

}

// -----------------------------------------------------------------------------------------------------------
