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
#include <LikelihoodTTGamma_RadTopProd.h> 
#include <LikelihoodTTGamma_HadWRadDecay.h>
#include <LikelihoodTTGamma_LepWRadDecay.h>
#include <LikelihoodTTGamma_HadTopRadDecay.h>
#include <LikelihoodTTGamma_LepTopRadDecay.h>
#include <MatchingTool.h> 
#include <SelectionTool.h> 
#include <TString.h>
#include <TSystem.h>  
#include <PhotonType.h>

bool DO_ELECTRON = true;
bool DO_MUON     = false;
bool IS_BATCH    = false;

int GetCombinationNo(TString s);//int matchHadB, int matchLepB, int matchLJ1, int matchLJ2);
TString GetTrueCombinationString(std::vector<int> v0, std::vector<int> v1, std::vector<int> v2, std::vector<int> v3);


int main(int argc, char **argv) 
{
	// counters for matching
	unsigned int countEMatch = 0;
	unsigned int countJMatch = 0;
	unsigned int countGMatch = 0;

	// parameters (1st: 'first event (starts with 0)', 2nd: 'last event - 1'
	if (argc!=3) { // i.e. 2 parameters ...
		std::cout << "number of parameters is unequal 2" << std::endl;
		return 1;
	}

	// create new fitters
	KLFitter::Fitter * myFitter_RadTopProd     = new KLFitter::Fitter(); 
	KLFitter::Fitter * myFitter_HadTopRadDecay = new KLFitter::Fitter(); 
	KLFitter::Fitter * myFitter_LepTopRadDecay = new KLFitter::Fitter(); 
	KLFitter::Fitter * myFitter_HadWRadDecay   = new KLFitter::Fitter(); 
	KLFitter::Fitter * myFitter_LepWRadDecay   = new KLFitter::Fitter(); 

	// open Root file 
	//KLFitter::InterfaceRoot * myInterfaceRoot = new KLFitter::InterfaceGoTopTree(); 
	KLFitter::InterfaceRoot * myInterfaceRoot = new KLFitter::InterfaceDummy(); 
	myInterfaceRoot -> OpenRootFile("input.root");

	// create detector
	//      KLFitter::DetectorBase * myDetector = new KLFitter::DetectorAtlas(); 
	KLFitter::DetectorBase * myDetector = new KLFitter::DetectorDummy(); 
	if (!myFitter_RadTopProd     -> SetDetector(myDetector)) return 0; 
	if (!myFitter_HadTopRadDecay -> SetDetector(myDetector)) return 0; 
	if (!myFitter_LepTopRadDecay -> SetDetector(myDetector)) return 0; 
	if (!myFitter_HadWRadDecay   -> SetDetector(myDetector)) return 0; 
	if (!myFitter_LepWRadDecay   -> SetDetector(myDetector)) return 0; 
	
	// create likelihoods for different ttbargamma hypothesis
	KLFitter::LikelihoodTTGamma_RadTopProd     * myLikelihood_RadTopProd     = new KLFitter::LikelihoodTTGamma_RadTopProd();
	KLFitter::LikelihoodTTGamma_HadTopRadDecay * myLikelihood_HadTopRadDecay = new KLFitter::LikelihoodTTGamma_HadTopRadDecay();
	KLFitter::LikelihoodTTGamma_LepTopRadDecay * myLikelihood_LepTopRadDecay = new KLFitter::LikelihoodTTGamma_LepTopRadDecay();
	KLFitter::LikelihoodTTGamma_HadWRadDecay   * myLikelihood_HadWRadDecay   = new KLFitter::LikelihoodTTGamma_HadWRadDecay();
	KLFitter::LikelihoodTTGamma_LepWRadDecay   * myLikelihood_LepWRadDecay   = new KLFitter::LikelihoodTTGamma_LepWRadDecay();

	myLikelihood_RadTopProd     -> PhysicsConstants() -> SetMassTop(172.5); 
	myLikelihood_RadTopProd     -> SetFlagBTagging(false); 
	myLikelihood_RadTopProd     -> SetFlagIntegrate(false); 
	myLikelihood_RadTopProd     -> SetFlagTopMassFixed(true);
	myLikelihood_RadTopProd     -> SetFlagUseJetMass(false);
	myLikelihood_HadTopRadDecay -> PhysicsConstants() -> SetMassTop(172.5); 
	myLikelihood_HadTopRadDecay -> SetFlagBTagging(false); 
	myLikelihood_HadTopRadDecay -> SetFlagIntegrate(false); 
	myLikelihood_HadTopRadDecay -> SetFlagTopMassFixed(true);
	myLikelihood_HadTopRadDecay -> SetFlagUseJetMass(false);
	myLikelihood_LepTopRadDecay -> PhysicsConstants() -> SetMassTop(172.5); 
	myLikelihood_LepTopRadDecay -> SetFlagBTagging(false); 
	myLikelihood_LepTopRadDecay -> SetFlagIntegrate(false); 
	myLikelihood_LepTopRadDecay -> SetFlagTopMassFixed(true);
	myLikelihood_LepTopRadDecay -> SetFlagUseJetMass(false);
	myLikelihood_HadWRadDecay   -> PhysicsConstants() -> SetMassTop(172.5); 
	myLikelihood_HadWRadDecay   -> SetFlagBTagging(false); 
	myLikelihood_HadWRadDecay   -> SetFlagIntegrate(false); 
	myLikelihood_HadWRadDecay   -> SetFlagTopMassFixed(true);
	myLikelihood_HadWRadDecay   -> SetFlagUseJetMass(false);
	myLikelihood_LepWRadDecay   -> PhysicsConstants() -> SetMassTop(172.5); 
	myLikelihood_LepWRadDecay   -> SetFlagBTagging(false); 
	myLikelihood_LepWRadDecay   -> SetFlagIntegrate(false); 
	myLikelihood_LepWRadDecay   -> SetFlagTopMassFixed(true);
	myLikelihood_LepWRadDecay   -> SetFlagUseJetMass(false);

	myInterfaceRoot -> SetFlagIsSignalMC(true);

	if (DO_ELECTRON) {
		myLikelihood_RadTopProd     -> SetLeptonType(1); // set lepton type to electron 
		myLikelihood_HadTopRadDecay -> SetLeptonType(1);
		myLikelihood_LepTopRadDecay -> SetLeptonType(1);
		myLikelihood_HadWRadDecay   -> SetLeptonType(1);
		myLikelihood_LepWRadDecay   -> SetLeptonType(1);
	}
	if (DO_MUON) {
		myLikelihood_RadTopProd     -> SetLeptonType(2); // set lepton type to muon
		myLikelihood_HadTopRadDecay -> SetLeptonType(2);
		myLikelihood_LepTopRadDecay -> SetLeptonType(2);
		myLikelihood_HadWRadDecay   -> SetLeptonType(2);
		myLikelihood_LepWRadDecay   -> SetLeptonType(2);
	}

	if (!myFitter_RadTopProd     -> SetLikelihood(myLikelihood_RadTopProd    )) return 0; 
	if (!myFitter_HadTopRadDecay -> SetLikelihood(myLikelihood_HadTopRadDecay)) return 0; 
	if (!myFitter_LepTopRadDecay -> SetLikelihood(myLikelihood_LepTopRadDecay)) return 0; 
	if (!myFitter_HadWRadDecay   -> SetLikelihood(myLikelihood_HadWRadDecay  )) return 0; 
	if (!myFitter_LepWRadDecay   -> SetLikelihood(myLikelihood_LepWRadDecay  )) return 0; 

	if (IS_BATCH)
  	gSystem->Exec("cp /work/pcatlas07/erdmann/e+jets/input.root .");

	// create interfaces for output 
	KLFitter::InterfaceOutput * myInterfaceOutput_RadTopProd     = new KLFitter::InterfaceOutput(); 
	KLFitter::InterfaceOutput * myInterfaceOutput_HadTopRadDecay = new KLFitter::InterfaceOutput(); 
	KLFitter::InterfaceOutput * myInterfaceOutput_LepTopRadDecay = new KLFitter::InterfaceOutput(); 
	KLFitter::InterfaceOutput * myInterfaceOutput_HadWRadDecay   = new KLFitter::InterfaceOutput(); 
	KLFitter::InterfaceOutput * myInterfaceOutput_LepWRadDecay   = new KLFitter::InterfaceOutput(); 

	// open Root files
	if (!myInterfaceOutput_RadTopProd     -> OpenRootFile("output_RadTopProd.root"))     return 0; 
	if (!myInterfaceOutput_HadTopRadDecay -> OpenRootFile("output_HadTopRadDecay.root")) return 0; 
	if (!myInterfaceOutput_LepTopRadDecay -> OpenRootFile("output_LepTopRadDecay.root")) return 0; 
	if (!myInterfaceOutput_HadWRadDecay   -> OpenRootFile("output_HadWRadDecay.root"))   return 0; 
	if (!myInterfaceOutput_LepWRadDecay   -> OpenRootFile("output_LepWRadDecay.root"))   return 0; 

	// create selection tool
	KLFitter::SelectionTool * mySelectionTool = new KLFitter::SelectionTool(); 
	mySelectionTool -> RequireNJetsPt(20.0, 4, -1); 
	//	mySelectionTool -> RequireNJetsPt(40.0, 3); 
	if (DO_ELECTRON)
		mySelectionTool -> RequireNElectronsPt(20.0, 1, -1); 
	if (DO_MUON)
		mySelectionTool -> RequireNMuonsPt(20.0, 1, -1);
	mySelectionTool -> RequireNPhotonsPt(10.0, 1, -1);
	mySelectionTool -> RequireMET(20.); 

	// create matching tool
	KLFitter::MatchingTool * myMatchingTool = new KLFitter::MatchingTool( myFitter_RadTopProd -> PParticles(), myInterfaceRoot -> PParticlesTruth() ); 

	// set fitter and truth particles 
	myInterfaceOutput_RadTopProd     -> SetFitter(myFitter_RadTopProd); 
	myInterfaceOutput_RadTopProd     -> SetParticlesTruth( myInterfaceRoot -> PParticlesTruth() ); 
	myInterfaceOutput_RadTopProd     -> SetParticlesMeasured( myInterfaceRoot -> PParticles() ); 
	myInterfaceOutput_RadTopProd     -> SetMatchingTool(myMatchingTool); 
	myInterfaceOutput_RadTopProd     -> SetSelectionTool(mySelectionTool); 
	myInterfaceOutput_HadTopRadDecay -> SetFitter(myFitter_HadTopRadDecay); 
	myInterfaceOutput_HadTopRadDecay -> SetParticlesTruth( myInterfaceRoot -> PParticlesTruth() ); 
	myInterfaceOutput_HadTopRadDecay -> SetParticlesMeasured( myInterfaceRoot -> PParticles() ); 
	myInterfaceOutput_HadTopRadDecay -> SetMatchingTool(myMatchingTool); 
	myInterfaceOutput_HadTopRadDecay -> SetSelectionTool(mySelectionTool); 
	myInterfaceOutput_LepTopRadDecay -> SetFitter(myFitter_LepTopRadDecay); 
	myInterfaceOutput_LepTopRadDecay -> SetParticlesTruth( myInterfaceRoot -> PParticlesTruth() ); 
	myInterfaceOutput_LepTopRadDecay -> SetParticlesMeasured( myInterfaceRoot -> PParticles() ); 
	myInterfaceOutput_LepTopRadDecay -> SetMatchingTool(myMatchingTool); 
	myInterfaceOutput_LepTopRadDecay -> SetSelectionTool(mySelectionTool); 
	myInterfaceOutput_HadWRadDecay   -> SetFitter(myFitter_HadWRadDecay); 
	myInterfaceOutput_HadWRadDecay   -> SetParticlesTruth( myInterfaceRoot -> PParticlesTruth() ); 
	myInterfaceOutput_HadWRadDecay   -> SetParticlesMeasured( myInterfaceRoot -> PParticles() ); 
	myInterfaceOutput_HadWRadDecay   -> SetMatchingTool(myMatchingTool); 
	myInterfaceOutput_HadWRadDecay   -> SetSelectionTool(mySelectionTool); 
	myInterfaceOutput_LepWRadDecay   -> SetFitter(myFitter_LepWRadDecay); 
	myInterfaceOutput_LepWRadDecay   -> SetParticlesTruth( myInterfaceRoot -> PParticlesTruth() ); 
	myInterfaceOutput_LepWRadDecay   -> SetParticlesMeasured( myInterfaceRoot -> PParticles() ); 
	myInterfaceOutput_LepWRadDecay   -> SetMatchingTool(myMatchingTool); 
	myInterfaceOutput_LepWRadDecay   -> SetSelectionTool(mySelectionTool); 

	// make truth photon classifyer
	KLFitter::PhotonType * photonType = new KLFitter::PhotonType();
	photonType->SetDeltaTopMass(8.0);
	photonType->SetDeltaWMass(8.0);
	photonType->SetTruthParticles(myInterfaceRoot->PParticlesTruth());
	photonType->SetPhysicsConstants(myLikelihood_RadTopProd->PhysicsConstants());

	// switch off filling of histograms
	//myFitter -> Likelihood() -> MCMCSetFlagFillHistograms(false);

	int nevents = myInterfaceRoot -> NEvents(); 
	int minEv = TString(argv[1]).Atoi();
	int maxEv = TString(argv[2]).Atoi();
	if (maxEv>=nevents) {
		maxEv = nevents;
		std::cout << "parameter warning: last event parameter reset to maximum event number available (" << nevents << ")" << std::endl;
	}
	if (minEv>=maxEv) {
		std::cout << "parameter error: first event larger than last event  -->  exiting" << std::endl;
		return 1;
	}

	// loop over events
	for (int ievent = minEv; ievent < maxEv; ++ievent)
		{
			if ((ievent-minEv+1)%50 == 0)
				std::cout << " event " << (ievent+1) << std::endl; 

			// get first event
			if (!myInterfaceRoot -> Event(ievent))
				return 0; 

			// get the event weight
			double weight = myInterfaceRoot->Weight();

			// read single event from root file and get particles 
			KLFitter::Particles * measuredparticles = myInterfaceRoot -> Particles(); 
			KLFitter::Particles * truthparticles = myInterfaceRoot -> ParticlesTruth(); 

			// select event
			if (!mySelectionTool -> SelectEvent(measuredparticles, myInterfaceRoot -> ET_miss())) 
				continue; 

			// get particles from selection tool
			KLFitter::Particles * particles = mySelectionTool -> ParticlesSelected(); 

			// add particles to fitters
			if (!myFitter_RadTopProd     -> SetParticles(particles)) return 0; 	
			if (!myFitter_HadTopRadDecay -> SetParticles(particles)) return 0; 	
			if (!myFitter_LepTopRadDecay -> SetParticles(particles)) return 0; 	
			if (!myFitter_HadWRadDecay   -> SetParticles(particles)) return 0; 	
			if (!myFitter_LepWRadDecay   -> SetParticles(particles)) return 0; 	

			// add ETmiss to fitters
			if (!myFitter_RadTopProd     -> SetET_miss_XY( myInterfaceRoot -> ET_miss_x(), myInterfaceRoot -> ET_miss_y() ) ) return 0;
			if (!myFitter_HadTopRadDecay -> SetET_miss_XY( myInterfaceRoot -> ET_miss_x(), myInterfaceRoot -> ET_miss_y() ) ) return 0;
			if (!myFitter_LepTopRadDecay -> SetET_miss_XY( myInterfaceRoot -> ET_miss_x(), myInterfaceRoot -> ET_miss_y() ) ) return 0;
			if (!myFitter_HadWRadDecay   -> SetET_miss_XY( myInterfaceRoot -> ET_miss_x(), myInterfaceRoot -> ET_miss_y() ) ) return 0;
			if (!myFitter_LepWRadDecay   -> SetET_miss_XY( myInterfaceRoot -> ET_miss_x(), myInterfaceRoot -> ET_miss_y() ) ) return 0;

			// perform matching
			myMatchingTool -> MatchTruthAll(KLFitter::Particles::kParton); 
			if (DO_ELECTRON)
				myMatchingTool -> MatchTruthAll(KLFitter::Particles::kElectron); 
			if (DO_MUON)
				myMatchingTool -> MatchTruthAll(KLFitter::Particles::kMuon);
			myMatchingTool -> MatchTruthAll(KLFitter::Particles::kPhoton);

			// identify true permutation
			std::vector<int> v0 = myMatchingTool->ListMatchedTruth(0, KLFitter::Particles::kParton);
			std::vector<int> v1 = myMatchingTool->ListMatchedTruth(1, KLFitter::Particles::kParton);
			std::vector<int> v2 = myMatchingTool->ListMatchedTruth(2, KLFitter::Particles::kParton);
			std::vector<int> v3 = myMatchingTool->ListMatchedTruth(3, KLFitter::Particles::kParton);
//			int trueCombi = GetCombinationNo(GetTrueCombinationString(v0, v1, v2, v3));
//
//			// require the event to be matched
//			if (trueCombi==-1)
//				continue;

			// classify truth photon
			photonType->Classify();

			// copy information into output classes
			myInterfaceOutput_RadTopProd     -> SetEventWeight(weight); 
			myInterfaceOutput_RadTopProd     -> SetPhotonType(photonType->IsNotClassified(),
																												photonType->IsRadTopProd(),
																												photonType->IsHadTopRadDecay(),
																												photonType->IsLepTopRadDecay(),
																												photonType->IsHadWRadDecay(),
																												photonType->IsLepWRadDecay());
			myInterfaceOutput_RadTopProd     -> FillTreeMeasured(); 
			myInterfaceOutput_RadTopProd     -> FillTreeSelected(); 
			myInterfaceOutput_RadTopProd     -> FillTreeTruth();
			myInterfaceOutput_RadTopProd     -> FillTreeMatching();
			myInterfaceOutput_RadTopProd     -> FillTreeMap(); 
			myInterfaceOutput_HadTopRadDecay -> SetPhotonType(photonType->IsNotClassified(),
																												photonType->IsRadTopProd(),
																												photonType->IsHadTopRadDecay(),
																												photonType->IsLepTopRadDecay(),
																												photonType->IsHadWRadDecay(),
																												photonType->IsLepWRadDecay());
			myInterfaceOutput_HadTopRadDecay -> SetEventWeight(weight); 
			myInterfaceOutput_HadTopRadDecay -> FillTreeMeasured(); 
			myInterfaceOutput_HadTopRadDecay -> FillTreeSelected(); 
			myInterfaceOutput_HadTopRadDecay -> FillTreeTruth();
			myInterfaceOutput_HadTopRadDecay -> FillTreeMatching();
			myInterfaceOutput_HadTopRadDecay -> FillTreeMap(); 
			myInterfaceOutput_LepTopRadDecay -> SetPhotonType(photonType->IsNotClassified(),
																												photonType->IsRadTopProd(),
																												photonType->IsHadTopRadDecay(),
																												photonType->IsLepTopRadDecay(),
																												photonType->IsHadWRadDecay(),
																												photonType->IsLepWRadDecay());
			myInterfaceOutput_LepTopRadDecay -> SetEventWeight(weight); 
			myInterfaceOutput_LepTopRadDecay -> FillTreeMeasured(); 
			myInterfaceOutput_LepTopRadDecay -> FillTreeSelected(); 
			myInterfaceOutput_LepTopRadDecay -> FillTreeTruth();
			myInterfaceOutput_LepTopRadDecay -> FillTreeMatching();
			myInterfaceOutput_LepTopRadDecay -> FillTreeMap(); 
			myInterfaceOutput_HadWRadDecay -> SetPhotonType(photonType->IsNotClassified(),
																											photonType->IsRadTopProd(),
																											photonType->IsHadTopRadDecay(),
																											photonType->IsLepTopRadDecay(),
																											photonType->IsHadWRadDecay(),
																											photonType->IsLepWRadDecay());
			myInterfaceOutput_HadWRadDecay -> SetEventWeight(weight); 
			myInterfaceOutput_HadWRadDecay -> FillTreeMeasured(); 
			myInterfaceOutput_HadWRadDecay -> FillTreeSelected(); 
			myInterfaceOutput_HadWRadDecay -> FillTreeTruth();
			myInterfaceOutput_HadWRadDecay -> FillTreeMatching();
			myInterfaceOutput_HadWRadDecay -> FillTreeMap(); 
			myInterfaceOutput_LepWRadDecay -> SetPhotonType(photonType->IsNotClassified(),
																											photonType->IsRadTopProd(),
																											photonType->IsHadTopRadDecay(),
																											photonType->IsLepTopRadDecay(),
																											photonType->IsHadWRadDecay(),
																											photonType->IsLepWRadDecay());
			myInterfaceOutput_LepWRadDecay -> SetEventWeight(weight); 
			myInterfaceOutput_LepWRadDecay -> FillTreeMeasured(); 
			myInterfaceOutput_LepWRadDecay -> FillTreeSelected(); 
			myInterfaceOutput_LepWRadDecay -> FillTreeTruth();
			myInterfaceOutput_LepWRadDecay -> FillTreeMatching();
			myInterfaceOutput_LepWRadDecay -> FillTreeMap(); 

			// loop over all permutations 
			for (int iperm = 0; iperm < myFitter_RadTopProd -> Permutations() -> NPermutations(); ++iperm)
				{
					// fit the first permutation
					myFitter_RadTopProd     -> Fit(iperm); 
					myFitter_HadTopRadDecay -> Fit(iperm); 
					myFitter_LepTopRadDecay -> Fit(iperm); 
					myFitter_HadWRadDecay   -> Fit(iperm);
					myFitter_LepWRadDecay   -> Fit(iperm); 

					// copy information into output class
					myInterfaceOutput_RadTopProd     -> FillTreeModelPermutation(); 
					myInterfaceOutput_HadTopRadDecay -> FillTreeModelPermutation(); 
					myInterfaceOutput_LepTopRadDecay -> FillTreeModelPermutation(); 
					myInterfaceOutput_HadWRadDecay   -> FillTreeModelPermutation(); 
					myInterfaceOutput_LepWRadDecay   -> FillTreeModelPermutation(); 
				}

			// fill tree
			myInterfaceOutput_RadTopProd     -> FillTrees(); 
			myInterfaceOutput_HadTopRadDecay -> FillTrees(); 
			myInterfaceOutput_LepTopRadDecay -> FillTrees(); 
			myInterfaceOutput_HadWRadDecay   -> FillTrees(); 
			myInterfaceOutput_LepWRadDecay   -> FillTrees(); 

			// check if the particles have been matched

			bool electronsMatched = true;
			for (int iElec = 0; iElec < truthparticles -> NElectrons(); iElec++) {
				int nMatchedElec = myMatchingTool -> NMatchedTruth(iElec, KLFitter::Particles::kElectron); 
				if (nMatchedElec == 0)
					electronsMatched = false;
			}

			if (!electronsMatched)
				continue;

			countEMatch++;

			bool jetsMatched = true;
			std::set<TString> set_listOfMatches;
			set_listOfMatches.clear();

			for (int iQuark = 0; iQuark < truthparticles -> NPartons(); iQuark++) {
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

			bool photonsMatched = true;
			for (int iPhoton = 0; iPhoton < truthparticles -> NPhotons(); iPhoton++) {
				int nMatchedPhoton = myMatchingTool -> NMatchedTruth(iPhoton, KLFitter::Particles::kPhoton); 
				if (nMatchedPhoton == 0)
					photonsMatched = false;
			}

			if (!photonsMatched)
				continue;

			countGMatch++;
		}

	// output cut flow 
	std::cout << " N (all)       : " << mySelectionTool -> CounterEvents() << std::endl;
	std::cout << " N (electrons) : " << mySelectionTool -> CounterElectrons() << std::endl;
	std::cout << " N (jets)      : " << mySelectionTool -> CounterJets() << std::endl;
	std::cout << " N (photons)   : " << mySelectionTool -> CounterPhotons() << std::endl;
	std::cout << " N (MET)       : " << mySelectionTool -> CounterMET() << std::endl;
	std::cout << " N (e matched) : " << countEMatch << std::endl;
	std::cout << " N (j matched) : " << countJMatch << std::endl;
	std::cout << " N (g matched) : " << countGMatch << std::endl;

	// close input file 
	if (!myInterfaceRoot -> CloseRootFile()) 
		return 0;

	// close output files
	if (!myInterfaceOutput_RadTopProd     -> CloseRootFile()) return 0; 
	if (!myInterfaceOutput_HadTopRadDecay -> CloseRootFile()) return 0; 
	if (!myInterfaceOutput_LepTopRadDecay -> CloseRootFile()) return 0; 
	if (!myInterfaceOutput_HadWRadDecay   -> CloseRootFile()) return 0; 
	if (!myInterfaceOutput_LepWRadDecay   -> CloseRootFile()) return 0; 

	// free memory 
	delete myInterfaceOutput_RadTopProd; 
	delete myInterfaceOutput_HadTopRadDecay; 
	delete myInterfaceOutput_LepTopRadDecay; 
	delete myInterfaceOutput_HadWRadDecay; 
	delete myInterfaceOutput_LepWRadDecay; 
	delete myInterfaceRoot; 
	delete myDetector; 
	delete myFitter_RadTopProd; 
	delete myFitter_HadTopRadDecay; 
	delete myFitter_LepTopRadDecay; 
	delete myFitter_HadWRadDecay; 
	delete myFitter_LepWRadDecay; 
	delete myLikelihood_RadTopProd; 
	delete myLikelihood_HadTopRadDecay; 
	delete myLikelihood_LepTopRadDecay; 
	delete myLikelihood_HadWRadDecay; 
	delete myLikelihood_LepWRadDecay; 

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
