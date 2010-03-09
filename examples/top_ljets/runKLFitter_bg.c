#include <iostream> 
#include <vector>
#include <set>
#include <Fitter.h> 
#include <PhysicsConstants.h> 
#include <DetectorAtlas.h> 
#include <DetectorDummy.h> 
#include <InterfaceRoot.h> 
#include <InterfaceDummy.h> 
#include <InterfaceGoTopTree.h> 
#include <Particles.h> 
#include <InterfaceOutput.h> 
#include <LikelihoodTopLeptonJets.h> 
#include <MatchingTool.h> 
#include <SelectionTool.h> 
#include <TString.h>
#include <TSystem.h>  


bool DO_ELECTRON = true;
bool DO_MUON     = false;
bool IS_BATCH    = false;

int GetCombinationNo(TString s);//int matchHadB, int matchLepB, int matchLJ1, int matchLJ2);

int main(int argc, char **argv) 
{
	
	// parameters (1st: 'first event', 2nd: 'last event - 1'
	if (argc!=3) { // i.e. 2 parameters ...
		std::cout << "number of parameters is unequal 2" << std::endl;
		return 1;
	}

	// create new fitter 
	KLFitter::Fitter * myFitter = new KLFitter::Fitter(); 

	// open Root file 
	KLFitter::InterfaceRoot * myInterfaceRoot = new KLFitter::InterfaceGoTopTree(); 
	myInterfaceRoot -> OpenRootFile("/work1/onacken/data/KLFitter/input/W+njets/WenuP2345_r808_r838/input.root");
	//myInterfaceRoot -> OpenRootFile("/work1/onacken/data/KLFitter/input/W+njets/4jets/e/input.root");
	
	// create detector 
	KLFitter::DetectorBase * myDetector = new KLFitter::DetectorAtlas(); 
	if (!myFitter -> SetDetector(myDetector))
		return 0; 
	
	// create likelihood for ttbar->e+jets channel 
	KLFitter::LikelihoodTopLeptonJets * myLikelihood = new KLFitter::LikelihoodTopLeptonJets(); 

	myLikelihood -> PhysicsConstants() -> SetMassTop(172.5); 
	myLikelihood -> SetFlagBTagging(false); 
	myLikelihood -> SetFlagIntegrate(false); 
	myLikelihood -> SetFlagTopMassFixed(false);
	myLikelihood -> SetFlagUseJetMass(false);
	myLikelihood -> SetLeptonType(1); // set lepton type to electron 
	myInterfaceRoot -> SetFlagIsSignalMC(false);

	if (!myFitter -> SetLikelihood(myLikelihood))
		return 0; 

	if (IS_BATCH)
  	gSystem->Exec("cp /work/pcatlas07/erdmann/e+jets/input.root .");

	// create interface for output 
	KLFitter::InterfaceOutput * myInterfaceOutput = new KLFitter::InterfaceOutput(); 

	// open Root file 
	if (!myInterfaceOutput -> OpenRootFile("output.root"))
		return 0; 

	// create selection tool
	KLFitter::SelectionTool * mySelectionTool = new KLFitter::SelectionTool(); 
	mySelectionTool -> RequireNJetsPt(20.0, 4, -1); 
	mySelectionTool -> RequireNJetsPt(40.0, 3); 
	mySelectionTool -> RequireNElectronsPt(20.0, 1, -1); 
	mySelectionTool -> RequireMET(20.); 

	// set fitter and truth particles 
	myInterfaceOutput -> SetFitter(myFitter); 
	myInterfaceOutput -> SetParticlesMeasured( myInterfaceRoot -> PParticles() ); 
	myInterfaceOutput -> SetSelectionTool(mySelectionTool); 

	// switch off filling of histograms
	//myFitter -> Likelihood() -> MCMCSetFlagFillHistograms(false);

	int nevents = myInterfaceRoot -> NEvents(); 
	int minEv = TString(argv[1]).Atoi();
	int maxEv = TString(argv[2]).Atoi();
	if (maxEv>=nevents) {
		maxEv = nevents;
		std::cout << "parameter warning: last event parameter reset to maximum event number available" << std::endl;
	}
	if (minEv>=maxEv) {
		std::cout << "parameter error: first event larger than last event  -->  exiting" << std::endl;
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

	// output cut flow 
	std::cout << " N (all)       : " << mySelectionTool -> CounterEvents() << std::endl;
	std::cout << " N (electrons) : " << mySelectionTool -> CounterElectrons() << std::endl;
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

	// no error 
	return 1; 

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
