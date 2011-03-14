#include "InterfaceDummy.h" 
#include <iostream> 

// --------------------------------------------------------- 
KLFitter::InterfaceDummy::InterfaceDummy()
{
  fTree = 0; 

  fBhad_E = 0; 
  fBhad_px = 0; 
  fBhad_py = 0; 
  fBhad_pz = 0; 

  fQup_E = 0; 
  fQup_px = 0; 
  fQup_py = 0; 
  fQup_pz = 0; 

  fQdown_E = 0; 
  fQdown_px = 0; 
  fQdown_py = 0; 
  fQdown_pz = 0; 

  fBlep_E = 0; 
  fBlep_px = 0; 
  fBlep_py = 0; 
  fBlep_pz = 0; 

  fLepton_E = 0; 
  fLepton_px = 0; 
  fLepton_py = 0; 
  fLepton_pz = 0; 

  fPhoton_E = 0; 
  fPhoton_px = 0; 
  fPhoton_py = 0; 
  fPhoton_pz = 0; 

  MET_Et = 0; 
  MET_Etx = 0; 
  MET_Ety = 0; 
  MET_Phi = 0; 

  fTrue_Bhad_E = 0; 
  fTrue_Bhad_px = 0; 
  fTrue_Bhad_py = 0; 
  fTrue_Bhad_pz = 0; 

  fTrue_Qup_E = 0; 
  fTrue_Qup_px = 0; 
  fTrue_Qup_py = 0; 
  fTrue_Qup_pz = 0; 

  fTrue_Qdown_E = 0; 
  fTrue_Qdown_px = 0; 
  fTrue_Qdown_py = 0; 
  fTrue_Qdown_pz = 0; 

  fTrue_Blep_E = 0; 
  fTrue_Blep_px = 0; 
  fTrue_Blep_py = 0; 
  fTrue_Blep_pz = 0; 

  fTrue_Lepton_E = 0; 
  fTrue_Lepton_px = 0; 
  fTrue_Lepton_py = 0; 
  fTrue_Lepton_pz = 0; 

  fTrue_Photon_E = 0; 
  fTrue_Photon_px = 0; 
  fTrue_Photon_py = 0; 
  fTrue_Photon_pz = 0; 

  fTrue_Neutrino_E = 0; 
  fTrue_Neutrino_px = 0; 
  fTrue_Neutrino_py = 0; 
  fTrue_Neutrino_pz = 0; 

}

// --------------------------------------------------------- 
KLFitter::InterfaceDummy::~InterfaceDummy()
{
}

// --------------------------------------------------------- 
int KLFitter::InterfaceDummy::NEvents()
{
  if (!fTree)
    return 0; 
        
  else
    return fTree -> GetEntries(); 
}


// --------------------------------------------------------- 
int KLFitter::InterfaceDummy::OpenRootFile(const char * filename, Option_t * opt)
{
  // define error code 
  int err = 1; 

  // open file 
  err *= KLFitter::InterfaceRoot::OpenRootFile(filename, opt); 

  // connect Root tree 
  err *= this  -> ConnectTree("fTree"); 

  // return error code 
  return err; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceDummy::ConnectTree(const char * treename)
{
  // check if file exists 
  if (!fRootFile)
    {
      std::cout << "KLFitter::InterfaceDummy::ConnectTree(). No Root file defined." << std::endl; 
      return 0; 
    } 

  // check if file is open 
  if (!fRootFile -> IsOpen())
    { 
      std::cout << "KLFitter::InterfaceDummy::ConnectTree(). Root file not open."<< std::endl; 
      return 0; 
    }

  // get tree from file 
  fTree = (TTree *) fRootFile -> Get(treename); 

  if (!fTree)
    {
      std::cout << "KLFitter::InterfaceDummy::ConnectTree(). Tree not found." << std::endl; 
      return 0; 
    }

  // set branch addresses
  fTree -> SetBranchAddress("bhad_E",  &fBhad_E); 
  fTree -> SetBranchAddress("bhad_px", &fBhad_px); 
  fTree -> SetBranchAddress("bhad_py", &fBhad_py); 
  fTree -> SetBranchAddress("bhad_pz", &fBhad_pz); 

  fTree -> SetBranchAddress("blep_E",  &fBlep_E); 
  fTree -> SetBranchAddress("blep_px", &fBlep_px); 
  fTree -> SetBranchAddress("blep_py", &fBlep_py); 
  fTree -> SetBranchAddress("blep_pz", &fBlep_pz); 

  fTree -> SetBranchAddress("qup_E",  &fQup_E); 
  fTree -> SetBranchAddress("qup_px", &fQup_px); 
  fTree -> SetBranchAddress("qup_py", &fQup_py); 
  fTree -> SetBranchAddress("qup_pz", &fQup_pz); 

  fTree -> SetBranchAddress("qdown_E",  &fQdown_E); 
  fTree -> SetBranchAddress("qdown_px", &fQdown_px); 
  fTree -> SetBranchAddress("qdown_py", &fQdown_py); 
  fTree -> SetBranchAddress("qdown_pz", &fQdown_pz); 

  fTree -> SetBranchAddress("lcharged_E",  &fLepton_E); 
  fTree -> SetBranchAddress("lcharged_px", &fLepton_px); 
  fTree -> SetBranchAddress("lcharged_py", &fLepton_py); 
  fTree -> SetBranchAddress("lcharged_pz", &fLepton_pz); 

  fTree -> SetBranchAddress("photon_E",  &fPhoton_E); 
  fTree -> SetBranchAddress("photon_px", &fPhoton_px); 
  fTree -> SetBranchAddress("photon_py", &fPhoton_py); 
  fTree -> SetBranchAddress("photon_pz", &fPhoton_pz); 

  fTree -> SetBranchAddress("MET_Et",  &MET_Et); 
  fTree -> SetBranchAddress("MET_Phi", &MET_Phi); 
  fTree -> SetBranchAddress("MET_Etx", &MET_Etx); 
  fTree -> SetBranchAddress("MET_Ety", &MET_Ety); 

  fTree -> SetBranchAddress("true_bhad_E",  &fTrue_Bhad_E); 
  fTree -> SetBranchAddress("true_bhad_px", &fTrue_Bhad_px); 
  fTree -> SetBranchAddress("true_bhad_py", &fTrue_Bhad_py); 
  fTree -> SetBranchAddress("true_bhad_pz", &fTrue_Bhad_pz); 

  fTree -> SetBranchAddress("true_blep_E",  &fTrue_Blep_E); 
  fTree -> SetBranchAddress("true_blep_px", &fTrue_Blep_px); 
  fTree -> SetBranchAddress("true_blep_py", &fTrue_Blep_py); 
  fTree -> SetBranchAddress("true_blep_pz", &fTrue_Blep_pz); 

  fTree -> SetBranchAddress("true_qup_E",  &fTrue_Qup_E); 
  fTree -> SetBranchAddress("true_qup_px", &fTrue_Qup_px); 
  fTree -> SetBranchAddress("true_qup_py", &fTrue_Qup_py); 
  fTree -> SetBranchAddress("true_qup_pz", &fTrue_Qup_pz); 

  fTree -> SetBranchAddress("true_qdown_E",  &fTrue_Qdown_E); 
  fTree -> SetBranchAddress("true_qdown_px", &fTrue_Qdown_px); 
  fTree -> SetBranchAddress("true_qdown_py", &fTrue_Qdown_py); 
  fTree -> SetBranchAddress("true_qdown_pz", &fTrue_Qdown_pz); 

  fTree -> SetBranchAddress("true_lcharged_E",  &fTrue_Lepton_E); 
  fTree -> SetBranchAddress("true_lcharged_px", &fTrue_Lepton_px); 
  fTree -> SetBranchAddress("true_lcharged_py", &fTrue_Lepton_py); 
  fTree -> SetBranchAddress("true_lcharged_pz", &fTrue_Lepton_pz); 

  fTree -> SetBranchAddress("true_photon_E",  &fTrue_Photon_E); 
  fTree -> SetBranchAddress("true_photon_px", &fTrue_Photon_px); 
  fTree -> SetBranchAddress("true_photon_py", &fTrue_Photon_py); 
  fTree -> SetBranchAddress("true_photon_pz", &fTrue_Photon_pz); 

  fTree -> SetBranchAddress("true_lneutral_E",  &fTrue_Neutrino_E); 
  fTree -> SetBranchAddress("true_lneutral_px", &fTrue_Neutrino_px); 
  fTree -> SetBranchAddress("true_lneutral_py", &fTrue_Neutrino_py); 
  fTree -> SetBranchAddress("true_lneutral_pz", &fTrue_Neutrino_pz); 

  // no error 
  return 1;
}

// --------------------------------------------------------- 
int KLFitter::InterfaceDummy::Event(int index)
{
  // check tree 
  if (!fTree)
    {
      std::cout << "KLFitter::InterfaceDummy::GetEvent(). Tree not defined." << std::endl; 
      return 0; 
    } 

  // check event number 
  if (index < 0 || index >= fTree -> GetEntries())
    {
      std::cout << "KLFitter::InterfaceDummy::GetEvent(). Event number negative or too large." << std::endl; 
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
int KLFitter::InterfaceDummy::FillParticles()
{
  // delete old particles 
  if (fParticles)
    delete fParticles; 

  // delete old truth particles 
  if (fParticlesTruth)
    delete fParticlesTruth; 

  TLorentzVector * lv = 0; 

  // create new particle container
  fParticles = new KLFitter::Particles(); 

  // create Lorentz-vectors and add to list of particles 
  fParticles -> AddParticle(lv = new TLorentzVector(fBhad_px, fBhad_py, fBhad_pz, fBhad_E), KLFitter::Particles::kParton, "jet 1"); 
  fParticles -> AddParticle(lv = new TLorentzVector(fBlep_px, fBlep_py, fBlep_pz, fBlep_E), KLFitter::Particles::kParton, "jet 2"); 
  fParticles -> AddParticle(lv = new TLorentzVector(fQup_px, fQup_py, fQup_pz, fQup_E), KLFitter::Particles::kParton, "jet 3"); 
  fParticles -> AddParticle(lv = new TLorentzVector(fQdown_px, fQdown_py, fQdown_pz, fQdown_E), KLFitter::Particles::kParton, "jet 4"); 
  fParticles -> AddParticle(lv = new TLorentzVector(fLepton_px, fLepton_py, fLepton_pz, fLepton_E), KLFitter::Particles::kElectron, "electron"); 
  fParticles -> AddParticle(lv = new TLorentzVector(fPhoton_px, fPhoton_py, fPhoton_pz, fPhoton_E), KLFitter::Particles::kPhoton, "photon"); 

  // create new truth particles container 
  fParticlesTruth = new KLFitter::Particles(); 

  // create Lorentz-vectors and add to list of particles 
  fParticlesTruth -> AddParticle(lv = new TLorentzVector(fTrue_Bhad_px, fTrue_Bhad_py, fTrue_Bhad_pz, fTrue_Bhad_E), KLFitter::Particles::kParton, "hadronic b quark"); 
  fParticlesTruth -> AddParticle(lv = new TLorentzVector(fTrue_Blep_px, fTrue_Blep_py, fTrue_Blep_pz, fTrue_Blep_E), KLFitter::Particles::kParton, "leptonic b quark"); 
  fParticlesTruth -> AddParticle(lv = new TLorentzVector(fTrue_Qup_px, fTrue_Qup_py, fTrue_Qup_pz, fTrue_Qup_E), KLFitter::Particles::kParton, "up-type quark"); 
  fParticlesTruth -> AddParticle(lv = new TLorentzVector(fTrue_Qdown_px, fTrue_Qdown_py, fTrue_Qdown_pz, fTrue_Qdown_E), KLFitter::Particles::kParton, "down-type quark"); 
  fParticlesTruth -> AddParticle(lv = new TLorentzVector(fTrue_Lepton_px, fTrue_Lepton_py, fTrue_Lepton_pz, fTrue_Lepton_E), KLFitter::Particles::kElectron, "electron"); 
  fParticlesTruth -> AddParticle(lv = new TLorentzVector(fTrue_Photon_px, fTrue_Photon_py, fTrue_Photon_pz, fTrue_Photon_E), KLFitter::Particles::kPhoton, "photon"); 
  fParticlesTruth -> AddParticle(lv = new TLorentzVector(fTrue_Neutrino_px, fTrue_Neutrino_py, fTrue_Neutrino_pz, fTrue_Neutrino_E), KLFitter::Particles::kNeutrino, "neutrino"); 

  // no error 
  return 1;
}

// --------------------------------------------------------- 

