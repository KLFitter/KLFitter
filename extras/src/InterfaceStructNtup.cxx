#include "InterfaceStructNtup.h" 
#include "Particles.h"

#include <TLorentzVector.h>

#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>

#include <iostream> 
#include <cmath>

// --------------------------------------------------------- 
KLFitter::InterfaceStructNtup::InterfaceStructNtup()
{
  fTree = 0; 

  event_entry = 0;

  GoodMuonVec_N = 0;
  GoodJetVec_N = 0;
  GoodElectronVec_N = 0;  	

  for(int i=0; i<13; ++i){ 
    GoodMuonVec_E[i] = 0.0;  
    GoodMuonVec_px[i] = 0.0;  
    GoodMuonVec_py[i] = 0.0;  
    GoodMuonVec_pz[i] = 0.0;  
    GoodMuonVec_p_T[i] = 0.0;  
    GoodMuonVec_eta[i] = 0.0;  
    GoodMuonVec_phi[i] = 0.0;  
    
    GoodElectronVec_E[i] = 0.0;  
    GoodElectronVec_px[i] = 0.0;  
    GoodElectronVec_py[i] = 0.0;  
    GoodElectronVec_pz[i] = 0.0;  
    GoodElectronVec_p_T[i] = 0.0;  
    GoodElectronVec_eta[i] = 0.0;  
    GoodElectronVec_phi[i] = 0.0;  

    GoodJetVec_E[i] = 0.0;  
    GoodJetVec_px[i] = 0.0;  
    GoodJetVec_py[i] = 0.0;  
    GoodJetVec_pz[i] = 0.0;  
    GoodJetVec_p_T[i] = 0.0;  
    GoodJetVec_eta[i] = 0.0;  
    GoodJetVec_phi[i] = 0.0;  
    GoodJetVec_weight_SV0[i] = 0.0;
  }  
    EtMiss = 0; 
    PhiMiss = 0; 
    PxMiss = 0; 
    PyMiss = 0; 
    
  // needed due to incomplete dictionaries in ROOT (reading vector<bool> from TTree won't work without)
  gROOT->ProcessLine("#include <vector>");
  gROOT->ProcessLine(".L loader.C+");
}

// --------------------------------------------------------- 
KLFitter::InterfaceStructNtup::~InterfaceStructNtup()
{
}

// --------------------------------------------------------- 
int KLFitter::InterfaceStructNtup::NEvents()
{
  if (!fTree)
    return 0; 
        
  else
    return fTree->GetEntries(); 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceStructNtup::OpenRootFile(const char* filename, Option_t* opt)
{
  // define error code 
  int err = 1; 

  // open file 
  err *= KLFitter::InterfaceRoot::OpenRootFile(filename, opt); 

  // connect Root tree 
  err *= this ->ConnectTree("tree"); 

  // return error code 
  return err; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceStructNtup::ConnectTree(const char* treename)
{
  // check if file exists 
  if (!fRootFile)
    {
      std::cout << "KLFitter::InterfaceStructNtup::ConnectTree(). No Root file defined." << std::endl; 
      return 0; 
    } 

  // check if file is open 
  if (!fRootFile->IsOpen())
    { 
      std::cout << "KLFitter::InterfaceStructNtup::ConnectTree(). Root file not open."<< std::endl; 
      return 0; 
    }

  // get tree from file 
  fTree = (TTree *) fRootFile->Get(treename); 

  if (!fTree)
    {
      std::cout << "KLFitter::InterfaceStructNtup::ConnectTree(). Tree not found." << std::endl; 
      return 0; 
    }

  // set branch addresses
  fTree->SetBranchAddress("event_entry",  &event_entry); 
  fTree->SetBranchAddress("eventWeight", &fWeight); 

  fTree->SetBranchAddress("GoodMuonVec_N",  &GoodMuonVec_N); 
  fTree->SetBranchAddress("GoodMuonVec_E",  &GoodMuonVec_E); 
  fTree->SetBranchAddress("GoodMuonVec_px", &GoodMuonVec_px); 
  fTree->SetBranchAddress("GoodMuonVec_py", &GoodMuonVec_py); 
  fTree->SetBranchAddress("GoodMuonVec_pz", &GoodMuonVec_pz); 
  fTree->SetBranchAddress("GoodMuonVec_p_T", &GoodMuonVec_p_T); 
  fTree->SetBranchAddress("GoodMuonVec_eta", &GoodMuonVec_eta); 
  fTree->SetBranchAddress("GoodMuonVec_phi", &GoodMuonVec_phi); 
 

  fTree->SetBranchAddress("GoodElectronVec_N",  &GoodElectronVec_N); 
  fTree->SetBranchAddress("GoodElectronVec_E",  &GoodElectronVec_E); 
  fTree->SetBranchAddress("GoodElectronVec_px", &GoodElectronVec_px); 
  fTree->SetBranchAddress("GoodElectronVec_py", &GoodElectronVec_py); 
  fTree->SetBranchAddress("GoodElectronVec_pz", &GoodElectronVec_pz); 
  fTree->SetBranchAddress("GoodElectronVec_p_T", &GoodElectronVec_p_T); 
  fTree->SetBranchAddress("GoodElectronVec_eta", &GoodElectronVec_eta); 
  fTree->SetBranchAddress("GoodElectronVec_phi", &GoodElectronVec_phi); 
  
  fTree->SetBranchAddress("GoodJetVec_N",   &GoodJetVec_N); 
  fTree->SetBranchAddress("GoodJetVec_E",   &GoodJetVec_E); 
  fTree->SetBranchAddress("GoodJetVec_px",  &GoodJetVec_px); 
  fTree->SetBranchAddress("GoodJetVec_py",  &GoodJetVec_py); 
  fTree->SetBranchAddress("GoodJetVec_pz",  &GoodJetVec_pz); 
  fTree->SetBranchAddress("GoodJetVec_p_T",  &GoodJetVec_p_T); 
  fTree->SetBranchAddress("GoodJetVec_eta", &GoodJetVec_eta); 
  fTree->SetBranchAddress("GoodJetVec_phi", &GoodJetVec_phi); 
  
  fTree->SetBranchAddress("GoodJetVec_weight_SV0", &GoodJetVec_weight_SV0); 

  fTree->SetBranchAddress("EtMiss",  &EtMiss); 
  fTree->SetBranchAddress("PhiMiss", &PhiMiss); 
  fTree->SetBranchAddress("PxMiss", &PxMiss); 
  fTree->SetBranchAddress("PyMiss", &PyMiss); 

  // no error 
  return 1;
}

// --------------------------------------------------------- 
int KLFitter::InterfaceStructNtup::Event(int index)
{

  // check tree 
  if (!fTree)
    {
      std::cout << "KLFitter::InterfaceStructNtup::GetEvent(). Tree not defined." << std::endl; 
      return 0; 
    } 

  // check event number 
  if (index < 0 || index >= fTree->GetEntries())
    {
      std::cout << "KLFitter::InterfaceStructNtup::GetEvent(). Event number negative or too large." << std::endl; 
      return 0; 
    } 

  // get event 
  fTree->GetEntry(index); 

  // fill particles 
  if (!this->FillParticles())
    return 0; 

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceStructNtup::FillParticles()
{
  // delete old particle container
  if (fParticles) 
    delete fParticles; 

  // create new particle container
  fParticles = new KLFitter::Particles(); 

  // fill Jets
  for (int i = 0; i < GoodJetVec_N; ++i)
    fParticles->AddParticle(new TLorentzVector(GoodJetVec_px[i] / 1000., GoodJetVec_py[i] / 1000., GoodJetVec_pz[i] / 1000., GoodJetVec_E[i] / 1000.), KLFitter::Particles::kParton); 

  for (int i = 0; i < GoodElectronVec_N; ++i)
    {
      //                  if(GoodElectronVec_Iso[i] / 1000.==false)
      //        continue;

      fParticles->AddParticle(new TLorentzVector(GoodElectronVec_px[i] / 1000., GoodElectronVec_py[i] / 1000., GoodElectronVec_pz[i] / 1000., GoodElectronVec_E[i] / 1000.), KLFitter::Particles::kElectron); 
    }

  // fill Muons
  for (int i = 0; i < GoodMuonVec_N; ++i)
    {
      //if(GoodMuonVec_Iso[i] / 1000.==false)
      //        continue;

      fParticles->AddParticle(new TLorentzVector(GoodMuonVec_px[i] / 1000., GoodMuonVec_py[i] / 1000., GoodMuonVec_pz[i] / 1000., GoodMuonVec_E[i] / 1000.), KLFitter::Particles::kMuon); 
    }

  // no error 
  return 1;
}

// --------------------------------------------------------- 
