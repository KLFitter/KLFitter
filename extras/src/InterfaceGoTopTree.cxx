#include "InterfaceGoTopTree.h" 
#include "Particles.h"

#include <TLorentzVector.h>

#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>

#include <iostream> 
#include <cmath>

// --------------------------------------------------------- 
KLFitter::InterfaceGoTopTree::InterfaceGoTopTree()
{
  fTree = 0; 

  Event_EventNumber = 0;

  Muon_N = 0;  
  Muon_E = 0;  
  Muon_Px = 0;  
  Muon_Py = 0;  
  Muon_Pz = 0;  
  Muon_Pt = 0;  
  Muon_Eta = 0;  
  Muon_Phi = 0;  
  Muon_IsTopInputs = 0;

  Electron_N = 0;  
  Electron_E = 0;  
  Electron_Px = 0;  
  Electron_Py = 0;  
  Electron_Pz = 0;  
  Electron_Pt = 0;  
  Electron_Eta = 0;  
  Electron_DetEta = 0;  
  Electron_Phi = 0;  
  Electron_IsTopInputs = 0;

  Jet_N = 0;  
  Jet_E = 0;  
  Jet_Px = 0;  
  Jet_Py = 0;  
  Jet_Pz = 0;  
  Jet_Pt = 0;  
  Jet_Eta = 0;
  Jet_DetEta = 0;  
  Jet_Phi = 0;  
  Jet_IsTopInputs = 0;

  Jet_SV0_Weight = 0;  

  Photon_N = 0;  
  Photon_E = 0;  
  Photon_Px = 0;  
  Photon_Py = 0;  
  Photon_Pz = 0;  
  Photon_Pt = 0;  
  Photon_Eta = 0;  
  Photon_Phi = 0;  
  Photon_IsTopInputs = 0;

  MET_Et = 0; 
  MET_Phi = 0; 
  MET_Etx = 0; 
  MET_Ety = 0; 

  Truth_IsProperMCEvent = false;

  TruthPart_N = 0; 
  TruthPart_PDG = 0x0; 
  TruthPart_NParents = 0x0; 
  TruthPart_ParentIdx = 0x0;
  TruthPart_NChildren = 0x0; 
  TruthPart_ChildIdx = 0x0;
  TruthPart_E = 0x0; 
  TruthPart_Px = 0x0; 
  TruthPart_Py = 0x0; 
  TruthPart_Pz = 0x0; 
  TruthPart_Eta = 0x0; 
  TruthPart_Phi = 0x0; 
  TruthPart_Pt = 0x0;

  Truth_WplusHad = false; 
  Truth_WminusHad = false; 

  TruthIdx_b = -1; 
  TruthIdx_bbar = -1; 
  TruthIdx_Wminus = -1; 
  TruthIdx_Wplus = -1; 
  TruthIdx_lminus = -1; 
  TruthIdx_lplus = -1; 
  TruthIdx_n = -1; 
  TruthIdx_nbar = -1; 
  TruthIdx_t = -1; 
  TruthIdx_tbar = -1; 
  TruthIdx_QfromWminus = -1; 
  TruthIdx_QfromWplus = -1; 
  TruthIdx_QbarfromWminus = -1; 
  TruthIdx_QbarfromWplus = -1; 
  TruthIdx_photon = -1; 

  // needed due to incomplete dictionaries in ROOT (reading vector<bool> from TTree won't work without)
  gROOT->ProcessLine("#include <vector>");
  gROOT->ProcessLine(".L loader.C+");
}

// --------------------------------------------------------- 
KLFitter::InterfaceGoTopTree::~InterfaceGoTopTree()
{
}

// --------------------------------------------------------- 
int KLFitter::InterfaceGoTopTree::NEvents()
{
  if(fTree)
  	return fTree->GetEntries();
  if(fChain)
  	return fChain->GetEntries();	        
  else
    return 0; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceGoTopTree::OpenRootFile(const char* filename, Option_t* opt)
{
  // define error code 
  int err = 1; 

  // open file 
  err *= KLFitter::InterfaceRoot::OpenRootFile(filename, opt); 

  // connect Root tree 
  err *= this ->ConnectTree("GoTopTree"); 

  // return error code 
  return err; 
}

//--------------------------------------------------------- 
int KLFitter::InterfaceGoTopTree::OpenRootFiles(std::vector<std::string> filenames, Option_t* opt)
{
  // define error code 
  int err = 1; 
	
	fChain =  new TChain("GoTopTree");
  // open files
  for(unsigned int i=0; i<filenames.size(); i++){ 
  	err *= KLFitter::InterfaceRoot::OpenRootFile(filenames.at(i).c_str(), opt);
  	fChain->Add(filenames.at(i).c_str());
  }	 
		
  // connect Root tree 
  err *= this ->ConnectChain(fChain); 

  // return error code 
  return err; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceGoTopTree::ConnectTree(const char* treename)
{
  // check if file exists 
  if (!fRootFile)
    {
      std::cout << "KLFitter::InterfaceGoTopTree::ConnectTree(). No Root file defined." << std::endl; 
      return 0; 
    } 

  // check if file is open 
  if (!fRootFile->IsOpen())
    { 
      std::cout << "KLFitter::InterfaceGoTopTree::ConnectTree(). Root file not open."<< std::endl; 
      return 0; 
    }

  // get tree from file 
  fTree = (TTree *) fRootFile->Get(treename); 

  if (!fTree)
    {
      std::cout << "KLFitter::InterfaceGoTopTree::ConnectTree(). Tree not found." << std::endl; 
      return 0; 
    }

  // set branch addresses
  fTree->SetBranchAddress("Event_EventNumber",  &Event_EventNumber); 
  fTree->SetBranchAddress("Event_Weight", &fWeight); 

  fTree->SetBranchAddress("Muon_N",  &Muon_N); 
  fTree->SetBranchAddress("Muon_E",  &Muon_E); 
  fTree->SetBranchAddress("Muon_Px", &Muon_Px); 
  fTree->SetBranchAddress("Muon_Py", &Muon_Py); 
  fTree->SetBranchAddress("Muon_Pz", &Muon_Pz); 
  fTree->SetBranchAddress("Muon_Pt", &Muon_Pt); 
  fTree->SetBranchAddress("Muon_Eta", &Muon_Eta);
  fTree->SetBranchAddress("Muon_Phi", &Muon_Phi); 
  fTree->SetBranchAddress("Muon_IsTopInputs", &Muon_IsTopInputs); 

  fTree->SetBranchAddress("Electron_N",  &Electron_N); 
  fTree->SetBranchAddress("Electron_EgE",  &Electron_E); 
  fTree->SetBranchAddress("Electron_Px", &Electron_Px); 
  fTree->SetBranchAddress("Electron_Py", &Electron_Py); 
  fTree->SetBranchAddress("Electron_Pz", &Electron_Pz); 
  fTree->SetBranchAddress("Electron_EgPt", &Electron_Pt); 
  fTree->SetBranchAddress("Electron_EgEta", &Electron_Eta);
  fTree->SetBranchAddress("Electron_DetEta", &Electron_DetEta); 
  fTree->SetBranchAddress("Electron_EgPhi", &Electron_Phi); 
  fTree->SetBranchAddress("Electron_IsTopInputs", &Electron_IsTopInputs); 

  fTree->SetBranchAddress("Jet_N",   &Jet_N); 
  fTree->SetBranchAddress("Jet_E",   &Jet_E); 
  fTree->SetBranchAddress("Jet_Px",  &Jet_Px); 
  fTree->SetBranchAddress("Jet_Py",  &Jet_Py); 
  fTree->SetBranchAddress("Jet_Pz",  &Jet_Pz); 
  fTree->SetBranchAddress("Jet_Pt",  &Jet_Pt); 
  fTree->SetBranchAddress("Jet_Eta", &Jet_Eta);
  fTree->SetBranchAddress("Jet_EMscaleDetEta", &Jet_DetEta);
  fTree->SetBranchAddress("Jet_Phi", &Jet_Phi); 
  fTree->SetBranchAddress("Jet_IsTopInputs", &Jet_IsTopInputs); 

  fTree->SetBranchAddress("Jet_SV0_Weight", &Jet_SV0_Weight); 

  if (fTree->FindBranch("Photon_N"))
    fTree->SetBranchAddress("Photon_N",  &Photon_N); 
  if (fTree->FindBranch("Photon_E"))
    fTree->SetBranchAddress("Photon_E",  &Photon_E); 
  if (fTree->FindBranch("Photon_Px"))
    fTree->SetBranchAddress("Photon_Px", &Photon_Px); 
  if (fTree->FindBranch("Photon_Py"))
    fTree->SetBranchAddress("Photon_Py", &Photon_Py); 
  if (fTree->FindBranch("Photon_Pz"))
    fTree->SetBranchAddress("Photon_Pz", &Photon_Pz); 
  if (fTree->FindBranch("Photon_Pt"))
    fTree->SetBranchAddress("Photon_Pt", &Photon_Pt); 
  if (fTree->FindBranch("Photon_Eta"))
    fTree->SetBranchAddress("Photon_Eta", &Photon_Eta); 
  if (fTree->FindBranch("Photon_Phi"))
    fTree->SetBranchAddress("Photon_Phi", &Photon_Phi); 
  if (fTree->FindBranch("Photon_IsTopInputs"))
    fTree->SetBranchAddress("Photon_IsTopInputs", &Photon_IsTopInputs); 
        
  fTree->SetBranchAddress("MET_Et",  &MET_Et); 
  fTree->SetBranchAddress("MET_Phi", &MET_Phi); 
  fTree->SetBranchAddress("MET_Etx", &MET_Etx); 
  fTree->SetBranchAddress("MET_Ety", &MET_Ety); 

  fTree->SetBranchAddress("Truth_IsProperMCEvent", &Truth_IsProperMCEvent); 
  fTree->SetBranchAddress("TruthPart_N",           &TruthPart_N); 
  fTree->SetBranchAddress("TruthPart_PDG",         &TruthPart_PDG); 


  fTree->SetBranchAddress("TruthPart_NParents", &TruthPart_NParents); 
  fTree->SetBranchAddress("TruthPart_ParentIdx", &TruthPart_ParentIdx); 
  fTree->SetBranchAddress("TruthPart_NChildren", &TruthPart_NChildren); 
  fTree->SetBranchAddress("TruthPart_ChildIdx", &TruthPart_ChildIdx); 
  fTree->SetBranchAddress("TruthPart_E",   &TruthPart_E); 
  fTree->SetBranchAddress("TruthPart_Px",  &TruthPart_Px); 
  fTree->SetBranchAddress("TruthPart_Py",  &TruthPart_Py); 
  fTree->SetBranchAddress("TruthPart_Pz",  &TruthPart_Pz); 
  fTree->SetBranchAddress("TruthPart_Eta", &TruthPart_Eta); 
  fTree->SetBranchAddress("TruthPart_Phi", &TruthPart_Phi); 
  fTree->SetBranchAddress("TruthPart_Pt",  &TruthPart_Pt); 

  fTree->SetBranchAddress("Truth_WplusHad",  &Truth_WplusHad); 
  fTree->SetBranchAddress("Truth_WminusHad",  &Truth_WminusHad); 

  fTree->SetBranchAddress("TruthIdx_Wplus",  &TruthIdx_Wplus); 
  fTree->SetBranchAddress("TruthIdx_b",      &TruthIdx_b); 
  fTree->SetBranchAddress("TruthIdx_bbar",   &TruthIdx_bbar); 
  fTree->SetBranchAddress("TruthIdx_lminus", &TruthIdx_lminus); 
  fTree->SetBranchAddress("TruthIdx_lplus",  &TruthIdx_lplus); 
  fTree->SetBranchAddress("TruthIdx_n",      &TruthIdx_n); 
  fTree->SetBranchAddress("TruthIdx_nbar",   &TruthIdx_nbar); 
  fTree->SetBranchAddress("TruthIdx_t",      &TruthIdx_t); 
  fTree->SetBranchAddress("TruthIdx_Wminus", &TruthIdx_Wminus); 
  fTree->SetBranchAddress("TruthIdx_tbar",   &TruthIdx_tbar); 
                    
  if (fTree->FindBranch("TruthIdx_photon"))
    fTree->SetBranchAddress("TruthIdx_photon", &TruthIdx_photon);
  fTree->SetBranchAddress("TruthIdx_QfromWminus", &TruthIdx_QfromWminus); 
  fTree->SetBranchAddress("TruthIdx_QfromWplus", &TruthIdx_QfromWplus); 
  fTree->SetBranchAddress("TruthIdx_QbarfromWminus", &TruthIdx_QbarfromWminus); 
  fTree->SetBranchAddress("TruthIdx_QbarfromWplus", &TruthIdx_QbarfromWplus); 

  // no error 
  return 1;
}
// --------------------------------------------------------- 
int KLFitter::InterfaceGoTopTree::ConnectChain(TChain * fChain)
{
	if (!this->fChain) this->fChain = fChain;
  
  // set branch addresses
  fChain->SetBranchAddress("Event_EventNumber",  &Event_EventNumber); 
  fChain->SetBranchAddress("Event_Weight", &fWeight); 

  fChain->SetBranchAddress("Muon_N",  &Muon_N); 
  fChain->SetBranchAddress("Muon_E",  &Muon_E); 
  fChain->SetBranchAddress("Muon_Px", &Muon_Px); 
  fChain->SetBranchAddress("Muon_Py", &Muon_Py); 
  fChain->SetBranchAddress("Muon_Pz", &Muon_Pz); 
  fChain->SetBranchAddress("Muon_Pt", &Muon_Pt); 
  fChain->SetBranchAddress("Muon_Eta", &Muon_Eta);
  fChain->SetBranchAddress("Muon_Phi", &Muon_Phi); 
  fChain->SetBranchAddress("Muon_IsTopInputs", &Muon_IsTopInputs); 

  fChain->SetBranchAddress("Electron_N",  &Electron_N); 
  fChain->SetBranchAddress("Electron_EgE",  &Electron_E); 
  fChain->SetBranchAddress("Electron_Px", &Electron_Px); 
  fChain->SetBranchAddress("Electron_Py", &Electron_Py); 
  fChain->SetBranchAddress("Electron_Pz", &Electron_Pz); 
  fChain->SetBranchAddress("Electron_EgPt", &Electron_Pt); 
  fChain->SetBranchAddress("Electron_EgEta", &Electron_Eta); 
  fChain->SetBranchAddress("Electron_DetEta", &Electron_DetEta); 
  fChain->SetBranchAddress("Electron_EgPhi", &Electron_Phi); 
  fChain->SetBranchAddress("Electron_IsTopInputs", &Electron_IsTopInputs); 

  fChain->SetBranchAddress("Jet_N",   &Jet_N); 
  fChain->SetBranchAddress("Jet_E",   &Jet_E); 
  fChain->SetBranchAddress("Jet_Px",  &Jet_Px); 
  fChain->SetBranchAddress("Jet_Py",  &Jet_Py); 
  fChain->SetBranchAddress("Jet_Pz",  &Jet_Pz); 
  fChain->SetBranchAddress("Jet_Pt",  &Jet_Pt); 
  fChain->SetBranchAddress("Jet_Eta", &Jet_Eta); 
  fChain->SetBranchAddress("Jet_EMscaleDetEta", &Jet_DetEta); 
  fChain->SetBranchAddress("Jet_Phi", &Jet_Phi); 
  fChain->SetBranchAddress("Jet_IsTopInputs", &Jet_IsTopInputs); 

  fChain->SetBranchAddress("Jet_SV0_Weight", &Jet_SV0_Weight); 

  if (fChain->FindBranch("Photon_N"))
    fChain->SetBranchAddress("Photon_N",  &Photon_N); 
  if (fChain->FindBranch("Photon_E"))
    fChain->SetBranchAddress("Photon_E",  &Photon_E); 
  if (fChain->FindBranch("Photon_Px"))
    fChain->SetBranchAddress("Photon_Px", &Photon_Px); 
  if (fChain->FindBranch("Photon_Py"))
    fChain->SetBranchAddress("Photon_Py", &Photon_Py); 
  if (fChain->FindBranch("Photon_Pz"))
    fChain->SetBranchAddress("Photon_Pz", &Photon_Pz); 
  if (fChain->FindBranch("Photon_Pt"))
    fChain->SetBranchAddress("Photon_Pt", &Photon_Pt); 
  if (fChain->FindBranch("Photon_Eta"))
    fChain->SetBranchAddress("Photon_Eta", &Photon_Eta); 
  if (fChain->FindBranch("Photon_Phi"))
    fChain->SetBranchAddress("Photon_Phi", &Photon_Phi); 
  if (fChain->FindBranch("Photon_IsTopInputs"))
    fChain->SetBranchAddress("Photon_IsTopInputs", &Photon_IsTopInputs); 
        
  fChain->SetBranchAddress("MET_Et",  &MET_Et); 
  fChain->SetBranchAddress("MET_Phi", &MET_Phi); 
  fChain->SetBranchAddress("MET_Etx", &MET_Etx); 
  fChain->SetBranchAddress("MET_Ety", &MET_Ety); 

  fChain->SetBranchAddress("Truth_IsProperMCEvent", &Truth_IsProperMCEvent); 
  fChain->SetBranchAddress("TruthPart_N",           &TruthPart_N); 
  fChain->SetBranchAddress("TruthPart_PDG",         &TruthPart_PDG); 


  fChain->SetBranchAddress("TruthPart_NParents", &TruthPart_NParents); 
  fChain->SetBranchAddress("TruthPart_ParentIdx", &TruthPart_ParentIdx); 
  fChain->SetBranchAddress("TruthPart_NChildren", &TruthPart_NChildren); 
  fChain->SetBranchAddress("TruthPart_ChildIdx", &TruthPart_ChildIdx); 
  fChain->SetBranchAddress("TruthPart_E",   &TruthPart_E); 
  fChain->SetBranchAddress("TruthPart_Px",  &TruthPart_Px); 
  fChain->SetBranchAddress("TruthPart_Py",  &TruthPart_Py); 
  fChain->SetBranchAddress("TruthPart_Pz",  &TruthPart_Pz); 
  fChain->SetBranchAddress("TruthPart_Eta", &TruthPart_Eta); 
  fChain->SetBranchAddress("TruthPart_Phi", &TruthPart_Phi); 
  fChain->SetBranchAddress("TruthPart_Pt",  &TruthPart_Pt); 

  fChain->SetBranchAddress("Truth_WplusHad",  &Truth_WplusHad); 
  fChain->SetBranchAddress("Truth_WminusHad",  &Truth_WminusHad); 

  fChain->SetBranchAddress("TruthIdx_Wplus",  &TruthIdx_Wplus); 
  fChain->SetBranchAddress("TruthIdx_b",      &TruthIdx_b); 
  fChain->SetBranchAddress("TruthIdx_bbar",   &TruthIdx_bbar); 
  fChain->SetBranchAddress("TruthIdx_lminus", &TruthIdx_lminus); 
  fChain->SetBranchAddress("TruthIdx_lplus",  &TruthIdx_lplus); 
  fChain->SetBranchAddress("TruthIdx_n",      &TruthIdx_n); 
  fChain->SetBranchAddress("TruthIdx_nbar",   &TruthIdx_nbar); 
  fChain->SetBranchAddress("TruthIdx_t",      &TruthIdx_t); 
  fChain->SetBranchAddress("TruthIdx_Wminus", &TruthIdx_Wminus); 
  fChain->SetBranchAddress("TruthIdx_tbar",   &TruthIdx_tbar); 
                    
  if (fChain->FindBranch("TruthIdx_photon"))
    fChain->SetBranchAddress("TruthIdx_photon", &TruthIdx_photon);
  fChain->SetBranchAddress("TruthIdx_QfromWminus", &TruthIdx_QfromWminus); 
  fChain->SetBranchAddress("TruthIdx_QfromWplus", &TruthIdx_QfromWplus); 
  fChain->SetBranchAddress("TruthIdx_QbarfromWminus", &TruthIdx_QbarfromWminus); 
  fChain->SetBranchAddress("TruthIdx_QbarfromWplus", &TruthIdx_QbarfromWplus); 

  // no error 
  return 1;
}

// --------------------------------------------------------- 
int KLFitter::InterfaceGoTopTree::Event(int index)
{

  // check tree 
  if (!fTree && !fChain)
    {
      std::cout << "KLFitter::InterfaceGoTopTree::GetEvent(). Tree not defined." << std::endl; 
      return 0; 
    } 

  if(fTree){
  	// check event number
  	if (index < 0 || index >= fTree->GetEntries())
    	{
      	std::cout << "KLFitter::InterfaceGoTopTree::GetEvent(). Event number negative or too large." << std::endl; 
      	return 0; 
    	} 
   	// get event 
  	fTree->GetEntry(index);
  }
  
  if(fChain){
  	// check event number
  	if (index < 0 || index >= fChain->GetEntries())
    	{
      	std::cout << "KLFitter::InterfaceGoTopTree::GetEvent(). Event number negative or too large." << std::endl; 
      	return 0; 
    	} 
    // get event 
  	fChain->GetEntry(index);
  } 

  // fill particles
  if (!this->FillParticles())
    return 0; 

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceGoTopTree::FillParticles()
{

  // delete old particle container
  if (fParticles) 
    delete fParticles; 

  // delete old truth particles container
  if (fParticlesTruth) 
    delete fParticlesTruth; 

  // create new particle container
  fParticles = new KLFitter::Particles(); 

  // fill jets
  for (int i = 0; i < Jet_N; ++i)
    {
      if ((*Jet_IsTopInputs)[i]){
        TLorentzVector * tlv_tmp = new TLorentzVector(0,0,0,0);
        tlv_tmp->SetPxPyPzE(Jet_Px->at(i), Jet_Py->at(i), Jet_Pz->at(i), Jet_E->at(i));
        fParticles->AddParticle(tlv_tmp, Jet_DetEta->at(i), KLFitter::Particles::kParton,"",Jet_SV0_Weight->at(i));
        delete tlv_tmp;
      }
    }
  // fill electrons
  for (int i = 0; i < Electron_N; ++i)
    {
      if ((*Electron_IsTopInputs)[i]){
        TLorentzVector * tlv_tmp = new TLorentzVector(0,0,0,0);
        tlv_tmp->SetPtEtaPhiE(Electron_Pt->at(i), Electron_Eta->at(i), Electron_Phi->at(i), Electron_E->at(i));
        fParticles->AddParticle(tlv_tmp, Electron_DetEta->at(i), KLFitter::Particles::kElectron); 
        delete tlv_tmp;
      } 
    }

  // fill muons
  for (int i = 0; i < Muon_N; ++i)
    {
      if ((*Muon_IsTopInputs)[i]){
        TLorentzVector * tlv_tmp = new TLorentzVector(0,0,0,0);
        tlv_tmp->SetPxPyPzE(Muon_Px->at(i), Muon_Py->at(i), Muon_Pz->at(i), Muon_E->at(i));
        fParticles->AddParticle(tlv_tmp, Muon_Eta->at(i), KLFitter::Particles::kMuon);
        delete tlv_tmp;
      }
    }

  // fill photons
  for (int i = 0; i < Photon_N; ++i)
    { 
      if ((*Photon_IsTopInputs)[i]){
        TLorentzVector * tlv_tmp = new TLorentzVector(0,0,0,0);
        tlv_tmp->SetPxPyPzE(Photon_Px->at(i), Photon_Py->at(i), Photon_Pz->at(i), Photon_E->at(i));
        fParticles->AddParticle(tlv_tmp, Photon_Eta->at(i), KLFitter::Particles::kPhoton);
      } 
    }

  // check if input is Signal MC
  if (!fFlagWriteSignalMCTruth)
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
  TLorentzVector * tlv_tmp = new TLorentzVector(0,0,0,0);
  tlv_tmp->SetPxPyPzE(TruthPart_Px->at(index_bhad),
                      TruthPart_Py->at(index_bhad), 
                      TruthPart_Pz->at(index_bhad), 
                      TruthPart_E->at(index_bhad));
  fParticlesTruth->AddParticle(tlv_tmp, KLFitter::Particles::kParton, "hadronic b quark");
  tlv_tmp->SetPxPyPzE(TruthPart_Px->at(index_blep), 
                      TruthPart_Py->at(index_blep), 
                      TruthPart_Pz->at(index_blep), 
                      TruthPart_E->at(index_blep));
  fParticlesTruth->AddParticle(tlv_tmp, KLFitter::Particles::kParton, "leptonic b quark");
  tlv_tmp->SetPxPyPzE(TruthPart_Px->at(index_q1), 
                      TruthPart_Py->at(index_q1), 
                      TruthPart_Pz->at(index_q1), 
                      TruthPart_E->at(index_q1));
  fParticlesTruth->AddParticle(tlv_tmp, KLFitter::Particles::kParton, "light quark 1");
  tlv_tmp->SetPxPyPzE(TruthPart_Px->at(index_q2), 
                      TruthPart_Py->at(index_q2), 
                      TruthPart_Pz->at(index_q2), 
                      TruthPart_E->at(index_q2));
  fParticlesTruth->AddParticle(tlv_tmp, KLFitter::Particles::kParton, "light quark 2"); 

  if (index_l!=-1 && abs(TruthPart_PDG->at(index_l)) == 11){
    tlv_tmp->SetPxPyPzE(TruthPart_Px->at(index_l), 
                        TruthPart_Py->at(index_l), 
                        TruthPart_Pz->at(index_l), 
                        TruthPart_E->at(index_l));
    fParticlesTruth->AddParticle(tlv_tmp, KLFitter::Particles::kElectron, "electron");
  } 
  else if (index_l!=-1 && abs(TruthPart_PDG->at(index_l)) == 13){
    tlv_tmp->SetPxPyPzE(TruthPart_Px->at(index_l),                                                                                             TruthPart_Py->at(index_l), 
                        TruthPart_Pz->at(index_l), 
                        TruthPart_E->at(index_l));
    fParticlesTruth->AddParticle(tlv_tmp, KLFitter::Particles::kMuon, "muon");
  } 
  else if (index_l!=-1 && abs(TruthPart_PDG->at(index_l)) == 15){
    tlv_tmp->SetPxPyPzE(TruthPart_Px->at(index_l), 
                        TruthPart_Py->at(index_l), 
                        TruthPart_Pz->at(index_l), 
                        TruthPart_E->at(index_l));
    fParticlesTruth->AddParticle(tlv_tmp, KLFitter::Particles::kTau, "tau");
  } 
  if (index_photon >= 0) {
    tlv_tmp->SetPxPyPzE(TruthPart_Px->at(index_photon), 
                        TruthPart_Py->at(index_photon), 
                        TruthPart_Pz->at(index_photon), 
                        TruthPart_E->at(index_photon));
    fParticlesTruth->AddParticle(tlv_tmp, KLFitter::Particles::kPhoton, "photon"); 
  }
  if (index_nu!=-1){
    tlv_tmp->SetPxPyPzE(TruthPart_Px->at(index_nu), 
                        TruthPart_Py->at(index_nu), 
                        TruthPart_Pz->at(index_nu), 
                        TruthPart_E->at(index_nu));
    fParticlesTruth->AddParticle(tlv_tmp, KLFitter::Particles::kNeutrino, "neutrino");
  } 
  if (index_tophad >= 0) {
    tlv_tmp->SetPxPyPzE(TruthPart_Px->at(index_tophad), 
                        TruthPart_Py->at(index_tophad), 
                        TruthPart_Pz->at(index_tophad), 
                        TruthPart_E->at(index_tophad));
    fParticlesTruth->AddParticle(tlv_tmp, KLFitter::Particles::kParton, "hadronic top quark"); 
  }

  if (index_toplep >= 0) {
    tlv_tmp->SetPxPyPzE(TruthPart_Px->at(index_toplep), 
                        TruthPart_Py->at(index_toplep), 
                        TruthPart_Pz->at(index_toplep), 
                        TruthPart_E->at(index_toplep));
    fParticlesTruth->AddParticle(tlv_tmp, KLFitter::Particles::kParton, "leptonic top quark"); 
  }
  //free memory
  delete tlv_tmp;

  // no error 
  return 1;

}

// --------------------------------------------------------- 
