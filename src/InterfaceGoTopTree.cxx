#include "InterfaceGoTopTree.h" 
#include <iostream> 
#include <cmath>

// --------------------------------------------------------- 
KLFitter::InterfaceGoTopTree::InterfaceGoTopTree()
{
	fTree = 0; 

	fWeight = 0;

	Muon_N = 0;  
	Muon_E = 0;  
	Muon_Px = 0;  
	Muon_Py = 0;  
	Muon_Pz = 0;  
	Muon_Pt = 0;  
	Muon_Eta = 0;  
	Muon_Phi = 0;  

	Electron_N = 0;  
	Electron_E = 0;  
	Electron_Px = 0;  
	Electron_Py = 0;  
	Electron_Pz = 0;  
	Electron_Pt = 0;  
	Electron_Eta = 0;  
	Electron_Phi = 0;  

	Jet_N = 0;  
	Jet_E = 0;  
	Jet_Px = 0;  
	Jet_Py = 0;  
	Jet_Pz = 0;  
	Jet_Pt = 0;  
	Jet_Eta = 0;  
	Jet_Phi = 0;  

	Photon_N = 0;  
	Photon_E = 0;  
	Photon_Px = 0;  
	Photon_Py = 0;  
	Photon_Pz = 0;  
	Photon_Pt = 0;  
	Photon_Eta = 0;  
	Photon_Phi = 0;  

	MET_Et = 0; 
	MET_Phi = 0; 
	MET_Etx = 0; 
	MET_Ety = 0; 

	TruthPart_N = 0; 

	/*
	fTruth_Botlep_N = 0;  
	fTruth_Botlep_E = 0;  
	fTruth_Botlep_Px = 0;  
	fTruth_Botlep_Py = 0;  
	fTruth_Botlep_Pz = 0;  
	fTruth_Botlep_Eta = 0;  
	fTruth_Botlep_Pt = 0;  
	fTruth_Botlep_Phi = 0;  

	fTruth_Bothad_N = 0;  
	fTruth_Bothad_E = 0;  
	fTruth_Bothad_Px = 0;  
	fTruth_Bothad_Py = 0;  
	fTruth_Bothad_Pz = 0;  
	fTruth_Bothad_Pt = 0;  
	fTruth_Bothad_Eta = 0;  
	fTruth_Bothad_Phi = 0;  

	fTruth_Tophad_N = 0;  
	fTruth_Tophad_E = 0;  
	fTruth_Tophad_Px = 0;  
	fTruth_Tophad_Py = 0;  
	fTruth_Tophad_Pz = 0;  
	fTruth_Tophad_Pt = 0;  
	fTruth_Tophad_Eta = 0;  
	fTruth_Tophad_Phi = 0;  

	fTruth_Toplep_N = 0;  
	fTruth_Toplep_E = 0;  
	fTruth_Toplep_Px = 0;  
	fTruth_Toplep_Py = 0;  
	fTruth_Toplep_Pz = 0;  
	fTruth_Toplep_Pt = 0;  
	fTruth_Toplep_Eta = 0;  
	fTruth_Toplep_Phi = 0;  

	fTruth_QfromW_N = 0;  
	fTruth_QfromW_E = 0;  
	fTruth_QfromW_Px = 0;  
	fTruth_QfromW_Py = 0;  
	fTruth_QfromW_Pz = 0;  
	fTruth_QfromW_Pt = 0;  
	fTruth_QfromW_Eta = 0;  
	fTruth_QfromW_Phi = 0;  

	fTruth_Muon_N = 0;  
	fTruth_Muon_E = 0;  
	fTruth_Muon_Px = 0;  
	fTruth_Muon_Py = 0;  
	fTruth_Muon_Pz = 0;  
	fTruth_Muon_Pt = 0;  
	fTruth_Muon_Eta = 0;  
	fTruth_Muon_Phi = 0;  

	fTruth_Tau_N = 0;  
	fTruth_Tau_E = 0;  
	fTruth_Tau_Px = 0;  
	fTruth_Tau_Py = 0;  
	fTruth_Tau_Pz = 0;  
	fTruth_Tau_Pt = 0;  
	fTruth_Tau_Eta = 0;  
	fTruth_Tau_Phi = 0;  

	fTruth_Electron_N = 0;  
	fTruth_Electron_E = 0;  
	fTruth_Electron_Px = 0;  
	fTruth_Electron_Py = 0;  
	fTruth_Electron_Pz = 0;  
	fTruth_Electron_Pt = 0;  
	fTruth_Electron_Eta = 0;  
	fTruth_Electron_Phi = 0;  

	fTruth_Neutrino_PDG = 0;  
	fTruth_Neutrino_N = 0;  
	fTruth_Neutrino_E = 0;  
	fTruth_Neutrino_Px = 0;  
	fTruth_Neutrino_Py = 0;  
	fTruth_Neutrino_Pz = 0;  
	fTruth_Neutrino_Pt = 0;  
	fTruth_Neutrino_Eta = 0;  
	fTruth_Neutrino_Phi = 0;  
	*/ 

	fLorentzVectorContainer = new std::vector <TLorentzVector *>(0); 
}

// --------------------------------------------------------- 
KLFitter::InterfaceGoTopTree::~InterfaceGoTopTree()
{
	while (!fLorentzVectorContainer -> empty())
		{
			TLorentzVector * lv = fLorentzVectorContainer -> front(); 
			fLorentzVectorContainer -> erase(fLorentzVectorContainer  -> begin() ); 
			delete lv; 
		}
	delete fLorentzVectorContainer; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceGoTopTree::NEvents()
{
	if (!fTree)
		return 0; 
	
	else
		return fTree -> GetEntries(); 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceGoTopTree::OpenRootFile(const char * filename, Option_t * opt)
{
	// define error code 
	int err = 1; 

	// open file 
	err *= KLFitter::InterfaceRoot::OpenRootFile(filename, opt); 

	// connect Root tree 
	err *= this  -> ConnectTree("GoTopTree"); 

	// return error code 
	return err; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceGoTopTree::ConnectTree(const char * treename)
{
	// check if file exists 
	if (!fRootFile)
		{
			std::cout << "KLFitter::InterfaceGoTopTree::ConnectTree(). No Root file defined." << std::endl; 
			return 0; 
		} 

	// check if file is open 
	if (!fRootFile -> IsOpen())
		{ 
			std::cout << "KLFitter::InterfaceGoTopTree::ConnectTree(). Root file not open."<< std::endl; 
			return 0; 
		}

	// get tree from file 
	fTree = (TTree *) fRootFile -> Get(treename); 

	if (!fTree)
		{
			std::cout << "KLFitter::InterfaceGoTopTree::ConnectTree(). Tree not found." << std::endl; 
			return 0; 
		}

	// set branch addresses
	fTree -> SetBranchAddress("Event_EventNumber",  &Event_EventNumber); 
	fTree -> SetBranchAddress("Event_Weight", &fWeight); 

	fTree -> SetBranchAddress("Muon_N",  &Muon_N); 
	fTree -> SetBranchAddress("Muon_E",  &Muon_E); 
	fTree -> SetBranchAddress("Muon_Px", &Muon_Px); 
	fTree -> SetBranchAddress("Muon_Py", &Muon_Py); 
	fTree -> SetBranchAddress("Muon_Pz", &Muon_Pz); 
	fTree -> SetBranchAddress("Muon_Pt", &Muon_Pt); 
	fTree -> SetBranchAddress("Muon_Eta", &Muon_Eta); 
	fTree -> SetBranchAddress("Muon_Phi", &Muon_Phi); 

	fTree -> SetBranchAddress("Electron_N",  &Electron_N); 
	fTree -> SetBranchAddress("Electron_E",  &Electron_E); 
	fTree -> SetBranchAddress("Electron_Px", &Electron_Px); 
	fTree -> SetBranchAddress("Electron_Py", &Electron_Py); 
	fTree -> SetBranchAddress("Electron_Pz", &Electron_Pz); 
	fTree -> SetBranchAddress("Electron_Pt", &Electron_Pt); 
	fTree -> SetBranchAddress("Electron_Eta", &Electron_Eta); 
	fTree -> SetBranchAddress("Electron_Phi", &Electron_Phi); 

	fTree -> SetBranchAddress("Jet_N",   &Jet_N); 
	fTree -> SetBranchAddress("Jet_E",   &Jet_E); 
	fTree -> SetBranchAddress("Jet_Px",  &Jet_Px); 
	fTree -> SetBranchAddress("Jet_Py",  &Jet_Py); 
	fTree -> SetBranchAddress("Jet_Pz",  &Jet_Pz); 
	fTree -> SetBranchAddress("Jet_Pt",  &Jet_Pt); 
	fTree -> SetBranchAddress("Jet_Eta", &Jet_Eta); 
	fTree -> SetBranchAddress("Jet_Phi", &Jet_Phi); 

	fTree -> SetBranchAddress("Photon_N",  &Photon_N); 
	fTree -> SetBranchAddress("Photon_E",  &Photon_E); 
	fTree -> SetBranchAddress("Photon_Px", &Photon_Px); 
	fTree -> SetBranchAddress("Photon_Py", &Photon_Py); 
	fTree -> SetBranchAddress("Photon_Pz", &Photon_Pz); 
	fTree -> SetBranchAddress("Photon_Pt", &Photon_Pt); 
	fTree -> SetBranchAddress("Photon_Eta", &Photon_Eta); 
	fTree -> SetBranchAddress("Photon_Phi", &Photon_Phi); 

	fTree -> SetBranchAddress("MET_Et",  &MET_Et); 
	fTree -> SetBranchAddress("MET_Phi", &MET_Phi); 
	fTree -> SetBranchAddress("MET_Etx", &MET_Etx); 
	fTree -> SetBranchAddress("MET_Ety", &MET_Ety); 

	fTree -> SetBranchAddress("Truth_IsProperMCEvent", &Truth_IsProperMCEvent); 
	fTree -> SetBranchAddress("TruthPart_N",           &TruthPart_N); 
	fTree -> SetBranchAddress("TruthPart_PDG",         &TruthPart_PDG); 


	fTree -> SetBranchAddress("TruthPart_NParents", &TruthPart_NParents); 
	fTree -> SetBranchAddress("TruthPart_ParentIdx", &TruthPart_ParentIdx); 
	fTree -> SetBranchAddress("TruthPart_NChildren", &TruthPart_NChildren); 
	fTree -> SetBranchAddress("TruthPart_ChildIdx", &TruthPart_ChildIdx); 
	fTree -> SetBranchAddress("TruthPart_E",   &TruthPart_E); 
	fTree -> SetBranchAddress("TruthPart_Px",  &TruthPart_Px); 
	fTree -> SetBranchAddress("TruthPart_Py",  &TruthPart_Py); 
	fTree -> SetBranchAddress("TruthPart_Pz",  &TruthPart_Pz); 
	fTree -> SetBranchAddress("TruthPart_Eta", &TruthPart_Eta); 
	fTree -> SetBranchAddress("TruthPart_Phi", &TruthPart_Phi); 
	fTree -> SetBranchAddress("TruthPart_Pt",  &TruthPart_Pt); 

	fTree -> SetBranchAddress("Truth_WplusHad",  &Truth_WplusHad); 
	fTree -> SetBranchAddress("Truth_WminusHad",  &Truth_WminusHad); 

	fTree -> SetBranchAddress("TruthIdx_Wplus",  &TruthIdx_Wplus); 
	fTree -> SetBranchAddress("TruthIdx_b",      &TruthIdx_b); 
	fTree -> SetBranchAddress("TruthIdx_bbar",   &TruthIdx_bbar); 
	fTree -> SetBranchAddress("TruthIdx_lminus", &TruthIdx_lminus); 
	fTree -> SetBranchAddress("TruthIdx_lplus",  &TruthIdx_lplus); 
	fTree -> SetBranchAddress("TruthIdx_n",      &TruthIdx_n); 
	fTree -> SetBranchAddress("TruthIdx_nbar",   &TruthIdx_nbar); 
	fTree -> SetBranchAddress("TruthIdx_t",      &TruthIdx_t); 
	fTree -> SetBranchAddress("TruthIdx_Wminus", &TruthIdx_Wminus); 
	fTree -> SetBranchAddress("TruthIdx_tbar",   &TruthIdx_tbar); 
	fTree -> SetBranchAddress("TruthIdx_photon", &TruthIdx_photon);
	fTree -> SetBranchAddress("TruthIdx_QfromWminus", &TruthIdx_QfromWminus); 
	fTree -> SetBranchAddress("TruthIdx_QfromWplus", &TruthIdx_QfromWplus); 
	fTree -> SetBranchAddress("TruthIdx_QbarfromWminus", &TruthIdx_QbarfromWminus); 
	fTree -> SetBranchAddress("TruthIdx_QbarfromWplus", &TruthIdx_QbarfromWplus); 

	// no error 
	return 1;
}

// --------------------------------------------------------- 
int KLFitter::InterfaceGoTopTree::Event(int index)
{

	// check tree 
	if (!fTree)
		{
			std::cout << "KLFitter::InterfaceGoTopTree::GetEvent(). Tree not defined." << std::endl; 
			return 0; 
		} 

	// check event number 
	if (index < 0 || index >= fTree -> GetEntries())
		{
			std::cout << "KLFitter::InterfaceGoTopTree::GetEvent(). Event number negative or too large." << std::endl; 
			return 0; 
		} 

	// get event 
	fTree -> GetEntry(index); 

	// fill particles 
	if (!this -> FillParticles())
		return 0; 

	// no error 
	return 1; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceGoTopTree::FillParticles()
{
	// delete four vectors 
	while (!fLorentzVectorContainer -> empty())
		{
			TLorentzVector * lv = fLorentzVectorContainer -> front(); 
			fLorentzVectorContainer -> erase(fLorentzVectorContainer  -> begin() ); 
			delete lv; 
		}

 	// delete old particle container
	if (fParticles)
		delete fParticles; 

	// delete old truth particles container
	if (fParticlesTruth)
		delete fParticlesTruth; 

	// create new particle container
	fParticles = new KLFitter::Particles(); 

	TLorentzVector * lv = 0; 

	// fill jets
	for (int i = 0; i < Jet_N; ++i)
		{
			fParticles -> AddParticle(lv = new TLorentzVector(Jet_Px -> at(i), Jet_Py -> at(i), Jet_Pz -> at(i), Jet_E -> at(i)), KLFitter::Particles::kParton); 

			// debugKK

			fLorentzVectorContainer -> push_back(lv); 
		}

	// fill electrons
	for (int i = 0; i < Electron_N; ++i)
		{
			fParticles -> AddParticle(lv = new TLorentzVector(Electron_Px -> at(i), Electron_Py -> at(i), Electron_Pz -> at(i), Electron_E -> at(i)), KLFitter::Particles::kElectron); 
			fLorentzVectorContainer -> push_back(lv); 
		}

	// fill muons
	for (int i = 0; i < Muon_N; ++i)
		{
			fParticles -> AddParticle(lv = new TLorentzVector(Muon_Px -> at(i), Muon_Py -> at(i), Muon_Pz -> at(i), Muon_E -> at(i)), KLFitter::Particles::kMuon); 
			fLorentzVectorContainer -> push_back(lv); 
		}

	// fill photons
	for (int i = 0; i < Photon_N; ++i)
		{
			fParticles -> AddParticle(lv = new TLorentzVector(Photon_Px -> at(i), Photon_Py -> at(i), Photon_Pz -> at(i), Photon_E -> at(i)), KLFitter::Particles::kPhoton); 
			fLorentzVectorContainer -> push_back(lv); 
		}

	// check if input is Signal MC
	if (!fFlagIsSignalMC)
		return 1; 
	
	// create truth particle container
	fParticlesTruth = new KLFitter::Particles(); 
	
	// check if event is proper ttbar event
	if (!Truth_IsProperMCEvent)
		return 1; 

	// do not fill mc information if di-leptonic
 	if ( (Truth_WplusHad == false) && (Truth_WminusHad == false) )
 		return 1; 

	int index_Whad = -1;
	int index_Wlep = -1;
	int index_tophad = -1; 
	int index_toplep = -1; 
	int index_bhad = -1; 
	int index_blep = -1; 
	int index_q1 = -1; 
	int index_q2 = -1; 
	int index_l = -1; 
	int index_nu = -1;
	int index_photon = -1;  

	if (Truth_WplusHad) 
		{
			index_Whad = TruthIdx_Wplus; 
			index_tophad = TruthIdx_t; 
			index_bhad = TruthIdx_b; 

			index_Wlep = TruthIdx_Wminus; 
			index_toplep = TruthIdx_tbar; 
			index_blep = TruthIdx_bbar; 
			index_l = TruthIdx_lminus;
			index_nu = TruthIdx_nbar; 
			index_q1 = TruthIdx_QfromWplus; 
			index_q2 = TruthIdx_QbarfromWplus; 

			index_photon = TruthIdx_photon;
		}
	else
		{
			index_Wlep = TruthIdx_Wplus; 
			index_toplep = TruthIdx_t; 
			index_blep = TruthIdx_b; 
			index_l = TruthIdx_lplus;
			index_nu = TruthIdx_n; 
			index_q1 = TruthIdx_QfromWminus; 
			index_q2 = TruthIdx_QbarfromWminus; 

			index_Whad = TruthIdx_Wminus; 
			index_tophad = TruthIdx_tbar; 
			index_bhad = TruthIdx_bbar; 

			index_photon = TruthIdx_photon;
		}

	// create Lorentz-vectors and add to list of particles 
	fParticlesTruth -> AddParticle(lv = new TLorentzVector(TruthPart_Px -> at(index_bhad), 
																												 TruthPart_Py -> at(index_bhad), 
																												 TruthPart_Pz -> at(index_bhad), 
																												 TruthPart_E -> at(index_bhad)), 
																 KLFitter::Particles::kParton, 
																 "hadronic b quark"); 
	fLorentzVectorContainer -> push_back(lv); 

	fParticlesTruth -> AddParticle(lv = new TLorentzVector(TruthPart_Px -> at(index_blep), 
																												 TruthPart_Py -> at(index_blep), 
																												 TruthPart_Pz -> at(index_blep), 
																												 TruthPart_E -> at(index_blep)), 
																 KLFitter::Particles::kParton, 
																 "leptonic b quark"); 
	fLorentzVectorContainer -> push_back(lv); 
	
	fParticlesTruth -> AddParticle(lv = new TLorentzVector(TruthPart_Px -> at(index_q1), 
																												 TruthPart_Py -> at(index_q1), 
																												 TruthPart_Pz -> at(index_q1), 
																												 TruthPart_E -> at(index_q1)), 
																 KLFitter::Particles::kParton, 
																 "light quark 1"); 
	fLorentzVectorContainer -> push_back(lv); 

	fParticlesTruth -> AddParticle(lv = new TLorentzVector(TruthPart_Px -> at(index_q2), 
																												 TruthPart_Py -> at(index_q2), 
																												 TruthPart_Pz -> at(index_q2), 
																												 TruthPart_E -> at(index_q2)), 
																 KLFitter::Particles::kParton, 
																 "light quark 2"); 
	fLorentzVectorContainer -> push_back(lv); 
 

	if (abs(TruthPart_PDG -> at(index_l)) == 11)
		fParticlesTruth -> AddParticle(lv = new TLorentzVector(TruthPart_Px -> at(index_l), 
																													 TruthPart_Py -> at(index_l), 
																													 TruthPart_Pz -> at(index_l), 
																													 TruthPart_E -> at(index_l)),
																	 KLFitter::Particles::kElectron, 
																	 "electron"); 
	else if (	abs(TruthPart_PDG -> at(index_l)) == 13)
		fParticlesTruth -> AddParticle(lv = new TLorentzVector(TruthPart_Px -> at(index_l), 
																													 TruthPart_Py -> at(index_l), 
																													 TruthPart_Pz -> at(index_l), 
																													 TruthPart_E -> at(index_l)),
																	 KLFitter::Particles::kMuon, 
																	 "muon"); 
	else if (	abs(TruthPart_PDG -> at(index_l)) == 15)
		fParticlesTruth -> AddParticle(lv = new TLorentzVector(TruthPart_Px -> at(index_l), 
																													 TruthPart_Py -> at(index_l), 
																													 TruthPart_Pz -> at(index_l), 
																													 TruthPart_E -> at(index_l)),
																	 KLFitter::Particles::kTau, 
																	 "tau"); 
	fLorentzVectorContainer -> push_back(lv); 

	fParticlesTruth -> AddParticle(lv = new TLorentzVector(TruthPart_Px -> at(index_photon), 
																												 TruthPart_Py -> at(index_photon), 
																												 TruthPart_Pz -> at(index_photon), 
																												 TruthPart_E -> at(index_photon)),
																 KLFitter::Particles::kPhoton, 
																 "photon"); 
	fLorentzVectorContainer -> push_back(lv); 

	fParticlesTruth -> AddParticle(lv = new TLorentzVector(TruthPart_Px -> at(index_nu), 
																												 TruthPart_Py -> at(index_nu), 
																												 TruthPart_Pz -> at(index_nu), 
																												 TruthPart_E -> at(index_nu)), 
																 KLFitter::Particles::kNeutrino, 
																 "neutrino"); 
	fLorentzVectorContainer -> push_back(lv); 

	if (index_tophad >= 0) {
		fParticlesTruth -> AddParticle(lv = new TLorentzVector(TruthPart_Px -> at(index_tophad), 
																													 TruthPart_Py -> at(index_tophad), 
																													 TruthPart_Pz -> at(index_tophad), 
																													 TruthPart_E -> at(index_tophad)), 
																	 KLFitter::Particles::kParton, 
																	 "hadronic top quark"); 
		fLorentzVectorContainer -> push_back(lv); 
	}

	if (index_toplep >= 0) {
		fParticlesTruth -> AddParticle(lv = new TLorentzVector(TruthPart_Px -> at(index_toplep), 
																													 TruthPart_Py -> at(index_toplep), 
																													 TruthPart_Pz -> at(index_toplep), 
																													 TruthPart_E -> at(index_toplep)), 
																	 KLFitter::Particles::kParton, 
																	 "leptonic top quark"); 
		fLorentzVectorContainer -> push_back(lv); 	
	}

	// no error 
	return 1;
}

// --------------------------------------------------------- 
