#include <iostream> 
#include <vector>
#include <set>
#include <Fitter.h> 
#include <PhysicsConstants.h> 
#include <DetectorAtlas.h> 
#include <DetectorDummy.h> 
#include <InterfaceRoot.h> 
#include <InterfaceDummy.h> 
#include <InterfaceDPD.h> 
#include <Particles.h> 
#include <InterfaceOutput.h> 
#include <LikelihoodTTHElectron.h> 
#include <MatchingTool.h> 
#include <TString.h>

bool EventSelection(KLFitter::Particles * particles, double ET_miss, std::vector<int> * cutflow);
bool EventTruthSelection(KLFitter::Particles * particles, std::vector<int> * cutflow); 
bool EventSummary(KLFitter::Particles * particles); 

int main()
{
	// create new fitter 
	KLFitter::Fitter * myFitter = new KLFitter::Fitter(); 

	// create detector 
	KLFitter::DetectorBase * myDetector = new KLFitter::DetectorAtlas(); 
	if (!myFitter -> SetDetector(myDetector))
		return 0; 
	
	// create likelihood for tt~H with H -> bb~ and tt~ -> e+jets channel 
	KLFitter::LikelihoodTTHElectron * myLikelihood = new KLFitter::LikelihoodTTHElectron(); 

	myLikelihood -> PhysicsConstants() -> SetMassTop(172.5); 
	myLikelihood -> SetFlagBTagging(false); 
	myLikelihood -> SetFlagIntegrate(false); 

	if (!myFitter -> SetLikelihood(myLikelihood))
		return 0; 

	// create interface for reading in data
	KLFitter::InterfaceRoot * myInterfaceRoot = new KLFitter::InterfaceDPD(); 
	
	// open Root file 
	//	if (!myInterfaceRoot -> OpenRootFile("/work/pcatlas19/dschiep/samples/ttbar/ttbar_e357_s462_r635.root"))
	if (!myInterfaceRoot -> OpenRootFile("input.root"))
		return 0;

	// create interface for output 
	KLFitter::InterfaceOutput * myInterfaceOutput = new KLFitter::InterfaceOutput(); 

	// open Root file 
	if (!myInterfaceOutput -> OpenRootFile("output.root"))
		return 0; 

	// create matching tool
	KLFitter::MatchingTool * myMatchingTool = new KLFitter::MatchingTool( myFitter -> PParticles(), myInterfaceRoot -> PParticlesTruth() ); 

	// set fitter and truth particles 
	myInterfaceOutput -> SetFitter(myFitter); 
	myInterfaceOutput -> SetParticlesTruth( myInterfaceRoot -> PParticlesTruth() ); 
	myInterfaceOutput -> SetMatchingTool(myMatchingTool); 

	// vector for keeping track of the cut flow
	std::vector<int> * cutflow = new std::vector<int>(7, 0);

	int nevents = myInterfaceRoot -> NEvents(); 
	//	nevents = 500;

	// loop over events
	for (int ievent = 0; ievent < nevents; ++ievent)
		{
			// get first event
			if (!myInterfaceRoot -> Event(ievent))
				return 0; 

			// get the event weight
			double weight = myInterfaceRoot->Weight();

			// read single event from root file and get particles 
			KLFitter::Particles * particles = myInterfaceRoot -> Particles(); 

			// read single event from root file and get particles 
			KLFitter::Particles * truthparticles = myInterfaceRoot -> ParticlesTruth(); 

			// cut flow variable
			cutflow->at(0)++;

			// truth event selection
			if (!EventTruthSelection(truthparticles, cutflow))
				continue; 

			// event selection
			if (!EventSelection(particles, myInterfaceRoot -> ET_miss() , cutflow))
				continue; 

			// print event summary 
			//			EventSummary(particles); 
			std::cout << " processing event " << ievent << " / " << nevents << std::endl; 

			// add particles to fitter 
			if (!myFitter -> SetParticles(particles))
				return 0; 	

			// perform matching
			myMatchingTool -> MatchTruthAll(KLFitter::Particles::kParton); 
			myMatchingTool -> MatchTruthAll(KLFitter::Particles::kElectron); 

			// copy information into output class
			myInterfaceOutput -> SetEventWeight(weight); 
			myInterfaceOutput -> FillTreeMeasured(); 
			myInterfaceOutput -> FillTreeTruth();
			myInterfaceOutput -> FillTreeMatching();

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

			// check if the particles have been matched

			bool electronsMatched = true;
			for (int iElec = 0; iElec < particles -> NElectrons(); iElec++) {
				int nMatchedElec = myMatchingTool -> NMatchedTruth(iElec, KLFitter::Particles::kElectron); 
				if (nMatchedElec == 0)
					electronsMatched = false;
			}

			if (!electronsMatched)
				continue;

			cutflow->at(5)++;

			bool jetsMatched = true;
			std::set<TString> set_listOfMatches;
			set_listOfMatches.clear();
			for (int iQuark = 0; iQuark < particles -> NPartons(); iQuark++) {
				int nMatchedJet = myMatchingTool -> NMatchedTruth(iQuark, KLFitter::Particles::kParton); 
				std::vector<int> listOfMatches = myMatchingTool -> ListMatchedTruth(iQuark, KLFitter::Particles::kParton);

				// check if matched is quark to exactly one jet
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

			cutflow->at(6)++;

		}

	// cut flow output
	printf("events before selection:      %5d\n", cutflow->at(0));
	printf("events after truth cuts:      %5d\n", cutflow->at(1));
	printf("events with 1 good electron:  %5d\n", cutflow->at(2));
	printf("events with 6 good jets:      %5d\n", cutflow->at(3));
	printf("events with MET > 20 GeV:     %5d\n", cutflow->at(4));
	printf("events with matched electron: %5d\n", cutflow->at(5));
	printf("events with 4 matched jets:   %5d\n", cutflow->at(6));

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
	delete cutflow;

	// no error 
	return 1; 

}

// -----------------------------------------------------------------------------------------------------------

bool EventSelection(KLFitter::Particles * particles, double ET_miss, std::vector<int> * cutflow)
{

	// ---------------------------------------------------------
	// electron selection 
	// ---------------------------------------------------------
	
	int nelectrons = particles -> NElectrons();

	// remove electrons outside of valied eta-region 
	for (int i = 0; i < nelectrons; ++i)
		{
			if (TMath::Abs(particles -> Electron(i) -> Eta()) > 2.5)
				{
					particles -> RemoveParticle(i, KLFitter::Particles::kElectron); 
					nelectrons--; 
					i--; 
					continue; 
				}
		}

	// remove electrons with pt < 20 GeV/c
	for (int i = 0; i < nelectrons; ++i)
		{
			if (TMath::Abs(particles -> Electron(i) -> Pt()) < 20.0)
				{
					particles -> RemoveParticle(i, KLFitter::Particles::kElectron); 
					nelectrons--; 
					i--; 
					continue; 
				}
		}
	
	if (nelectrons != 1)
		return 0; 

	cutflow->at(2)++;

	// ---------------------------------------------------------
	// jet selection 
	// ---------------------------------------------------------

	// get number of jets 
	int njets = particles -> NPartons(); 
	
	// remove jets outside of valied eta-region 
	for (int i = 0; i < njets; ++i)
		{
			if (TMath::Abs(particles -> Parton(i) -> Eta()) > 2.5)
				{
					particles -> RemoveParticle(i, KLFitter::Particles::kParton); 
					njets--; 
					i--; 
					continue; 
				}
		}

	// remove jets with pt < 20 GeV/c
	for (int i = 0; i < njets; ++i)
		{
			if (TMath::Abs(particles -> Parton(i) -> Pt()) < 20.0)
				{
					particles -> RemoveParticle(i, KLFitter::Particles::kParton); 
					njets--; 
					i--; 
					continue; 
				}
		}

	// count number of jets
	int njets_pt40 = 0; 

	for (int i = 0; i < njets; ++i)
		{
			if (TMath::Abs(particles -> Parton(i) -> Pt()) > 40.0)
				njets_pt40++; 
		}

	if (njets != 6)
		return 0; 

	// remove all jets but the first four
	// i.e. take only the four highest in pT into account
	//	while (njets>4)
	//		{
	//			particles -> RemoveParticle(njets-1, KLFitter::Particles::kParton); 
	//			njets--; 
	//			continue; 
	//		}

	cutflow->at(3)++;

	// ---------------------------------------------------------
	// ET miss selection 
	// ---------------------------------------------------------

	if (ET_miss < 20.)
		return 0;

	cutflow->at(4)++;

	// event pass
	return 1; 

}

// -----------------------------------------------------------------------------------------------------------

bool EventSummary(KLFitter::Particles * particles)
{
	// print summary 
	std::cout << "========================================================== " << std::endl; 
	std::cout << "Event summary " << std::endl; 
	std::cout << "========================================================== " << std::endl; 

	// print jets 
	if (particles -> NParticles(KLFitter::Particles::kParton) > 0) 
		{
			std::cout << "Jets: " << std::endl; 
			for (int i = 0; i < particles -> NParticles(KLFitter::Particles::kParton); ++i)
				std::cout << i << " " << particles -> NameParticle(i, KLFitter::Particles::kParton) << std::endl; 
		}

	// print electrons 
	if (particles -> NParticles(KLFitter::Particles::kElectron) > 0) 
		{
			std::cout << std::endl; 
			std::cout << "Electrons: " << std::endl; 
			for (int i = 0; i < particles -> NParticles(KLFitter::Particles::kElectron); ++i)
				std::cout << i << " " << particles -> NameParticle(i, KLFitter::Particles::kElectron) << std::endl; 
		}

	// print muons 
	if (particles -> NParticles(KLFitter::Particles::kMuon) > 0) 
		{
			std::cout << std::endl; 
			std::cout << "Muons: " << std::endl; 
			for (int i = 0; i < particles -> NParticles(KLFitter::Particles::kMuon); ++i)
				std::cout << i << " " << particles -> NameParticle(i, KLFitter::Particles::kMuon) << std::endl; 
		}

	// print neutrinos 
	if (particles -> NParticles(KLFitter::Particles::kNeutrino) > 0) 
		{
			std::cout << std::endl; 
			std::cout << "Neutrinos: " << std::endl; 
			for (int i = 0; i < particles -> NParticles(KLFitter::Particles::kNeutrino); ++i)
				std::cout << i << " " << particles -> NameParticle(i, KLFitter::Particles::kNeutrino) << std::endl; 
		}
	// print bosons 
	if (particles -> NParticles(KLFitter::Particles::kBoson) > 0) 
		{
			std::cout << std::endl; 
			std::cout << "Bosons: " << std::endl; 
			for (int i = 0; i < particles -> NParticles(KLFitter::Particles::kBoson); ++i)
				std::cout << i << " " << particles -> NameParticle(i, KLFitter::Particles::kBoson) << std::endl; 
		}
	std::cout << "========================================================== " << std::endl; 

	// no error 
	return 1; 
}

// -----------------------------------------------------------------------------------------------------------

bool EventTruthSelection(KLFitter::Particles * particles, std::vector<int> * cutflow)
{
	// ---------------------------------------------------------
	// truth selection
	// ---------------------------------------------------------

	int nTruthElectrons = particles -> NElectrons();
	int nTruthMuons = particles -> NMuons();
	int nTruthTaus = particles -> NTaus();
	int nTruthPartons = particles -> NPartons();

	// require exactly 1 truth electron
	if (nTruthElectrons != 1)
		return false;

	// require no truth muon or tau
	if (nTruthMuons != 0) return false;
	if (nTruthTaus != 0) return false;

	// require exactly 6 truth quarks (including the 2 tops)
	if (nTruthPartons !=6) return false;

	cutflow->at(1)++;

	return true;

}
