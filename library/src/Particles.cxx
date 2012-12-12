#include "Particles.h" 

#include <iostream> 

// --------------------------------------------------------- 
KLFitter::Particles::Particles() :
  fPartons(new std::vector <TLorentzVector *>(0)),  
  fElectrons(new std::vector <TLorentzVector *>(0)),  
  fMuons(new std::vector <TLorentzVector *>(0)),  
  fTaus(new std::vector <TLorentzVector *>(0)),  
  fNeutrinos(new std::vector <TLorentzVector *>(0)),  
  fBosons(new std::vector <TLorentzVector *>(0)),  
  fPhotons(new std::vector <TLorentzVector *>(0)),  

  fNamePartons(new std::vector <std::string>(0)),  
  fNameElectrons(new std::vector <std::string>(0)),  
  fNameMuons(new std::vector <std::string>(0)),  
  fNameTaus(new std::vector <std::string>(0)),  
  fNameNeutrinos(new std::vector <std::string>(0)),  
  fNameBosons(new std::vector <std::string>(0)),  
  fNamePhotons(new std::vector <std::string>(0)),  

  fJetIndex(new std::vector <int>(0)),  
  fElectronIndex(new std::vector <int>(0)),  
  fMuonIndex(new std::vector <int>(0)),  
  fPhotonIndex(new std::vector <int>(0)),  

  fTrueFlavor(new std::vector<TrueFlavorType>(0)),
  fIsBTagged(new std::vector<bool>(0)),
  fBTaggingEfficiency(new std::vector<double>(0)),
  fBTaggingRejection(new std::vector<double>(0)),

  fBTagWeight(new std::vector<double>(0)),
  fBTagWeightSet(new std::vector<bool>(0)),

  fElectronDetEta(new std::vector<double>(0)),
  fMuonDetEta(new std::vector<double>(0)),
  fJetDetEta(new std::vector<double>(0)),
  fPhotonDetEta(new std::vector<double>(0)),
  fElectronCharge(new std::vector<float>(0)),
  fMuonCharge(new std::vector<float>(0))
{
}

// --------------------------------------------------------- 
KLFitter::Particles::~Particles()
{
  while (!fPartons->empty()) {
    TLorentzVector* lv = fPartons->front();
    fPartons->erase( fPartons->begin() ); 
    if (lv) 
      delete lv;
  }
  delete fPartons; 

  while (!fElectrons->empty()) {
    TLorentzVector* lv = fElectrons->front();
    fElectrons->erase( fElectrons->begin() ); 
    if (lv)
      delete lv; 
  }
  delete fElectrons; 

  while (!fMuons->empty()) {
    TLorentzVector* lv = fMuons->front();
    fMuons->erase( fMuons->begin() ); 
    if (lv)
      delete lv; 
  }
  delete fMuons; 

  while (!fTaus->empty()) {
    TLorentzVector* lv = fTaus->front();
    fTaus->erase( fTaus->begin() ); 
    if (lv)
      delete lv; 
  }
  delete fTaus; 

  while (!fNeutrinos->empty()) {
    TLorentzVector* lv = fNeutrinos->front();
    fNeutrinos->erase( fNeutrinos->begin() ); 
    if (lv)
      delete lv; 
  }
  delete fNeutrinos; 

  while (!fBosons->empty()) {
    TLorentzVector* lv = fBosons->front();
    fBosons->erase( fBosons->begin() ); 
    if (lv) 
      delete lv; 
  }
  delete fBosons; 

  while (!fPhotons->empty()) {
    TLorentzVector* lv = fPhotons->front();
    fPhotons->erase( fPhotons->begin() ); 
    if (lv)
      delete lv; 
  }
  delete fPhotons; 

  if (fNamePartons)
    delete fNamePartons; 

  if (fNameElectrons)
    delete fNameElectrons; 

  if (fNameMuons)
    delete fNameMuons; 

  if (fNameTaus)
    delete fNameTaus; 

  if (fNameNeutrinos)
    delete fNameNeutrinos; 

  if (fNameBosons)
    delete fNameBosons;

  if (fNamePhotons)
    delete fNamePhotons; 

  if (fJetIndex)
    delete fJetIndex;

  if (fElectronIndex)
    delete fElectronIndex;

  if (fMuonIndex)
    delete fMuonIndex;

  if (fPhotonIndex)
    delete fPhotonIndex;

  if (fTrueFlavor)
    delete fTrueFlavor;

  if (fIsBTagged)
    delete fIsBTagged;

  if (fBTaggingEfficiency)
    delete fBTaggingEfficiency;

  if (fBTaggingRejection)
    delete fBTaggingRejection;

  if (fBTagWeight)
    delete fBTagWeight;

  if (fBTagWeightSet)
    delete fBTagWeightSet;

  if (fElectronDetEta)
    delete fElectronDetEta;

  if (fMuonDetEta)
    delete fMuonDetEta;

  if (fJetDetEta)
    delete fJetDetEta;

  if (fPhotonDetEta)
    delete fPhotonDetEta;

  if (fElectronCharge)
    delete fElectronCharge;

  if (fMuonCharge)
    delete fMuonCharge;
}

// --------------------------------------------------------- 
int KLFitter::Particles::AddParticle(TLorentzVector * particle, double DetEta, float LepCharge, KLFitter::Particles::ParticleType ptype, std::string name, int measuredindex)
{
  // get particle container
  std::vector <TLorentzVector *>* container = ParticleContainer(ptype); 

  //std::string name = "";
  //int measuredindex = -1;

  // check if container exists
  if (!container)
    {
      std::cout << "KLFitter::Particles::AddParticle(). Container does not exist." << std::endl; 
      return 0; 
    }

  // check name 
  if (name == "")
    name = Form("particle_%i", NParticles()); 

  // get index and type 
  TLorentzVector* vect = 0; 
  int index = 0; 
  KLFitter::Particles::ParticleType temptype = kParton; 

  // check if particle with name exists already 
  if (!FindParticle(name, vect, index, temptype)) {

    // add particle 
    // create pointer copy of particle content which is owend by Particles
    TLorentzVector * cparticle = new TLorentzVector(particle->Px(), particle->Py(), particle->Pz(), particle->E());
    container->push_back(cparticle); 
    ParticleNameContainer(ptype)->push_back(name); 
    if (ptype == KLFitter::Particles::kElectron){ 
      fElectronIndex->push_back(measuredindex);
      fElectronDetEta->push_back(DetEta);
      fElectronCharge->push_back(LepCharge);
    }
    else if (ptype == KLFitter::Particles::kMuon){ 
      fMuonIndex->push_back(measuredindex);
      fMuonDetEta->push_back(DetEta);
      fMuonCharge->push_back(LepCharge);
    } 
    else if (ptype == KLFitter::Particles::kPhoton){
      fPhotonIndex->push_back(measuredindex);
      fPhotonDetEta->push_back(DetEta);
    } 
  }
  else {
    std::cout << "KLFitter::Particles::AddParticle(). Particle with the name " << name << " exists already." << std::endl; 
    return 0; 
  }

  if (particle->M() < -1.e-3) {
    std::cout << "KLFitter::Particles::AddParticle(). WARNING : A particle with negative mass " << particle->M() << " of type " << ptype << " was added." << std::endl;
    return 1;
  }

  // no error
  return 1;
}

// --------------------------------------------------------- 
int KLFitter::Particles::AddParticle(TLorentzVector * particle, double DetEta, KLFitter::Particles::ParticleType ptype, std::string name, int measuredindex, bool isBtagged, double bTagEff, double bTagRej, TrueFlavorType trueflav, double btagweight)
{
  // get particle container
  std::vector <TLorentzVector *>* container = ParticleContainer(ptype); 

  // check if container exists
  if (!container)
    {
      std::cout << "KLFitter::Particles::AddParticle(). Container does not exist." << std::endl; 
      return 0; 
    }

  // check name 
  if (name == "")
    name = Form("particle_%i", NParticles()); 

  // get index and type 
  TLorentzVector* vect = 0; 
  int index = 0; 
  KLFitter::Particles::ParticleType temptype = kParton; 

  // check if particle with name exists already 
  if (!FindParticle(name, vect, index, temptype)) {

    // add particle 
    // create pointer copy of particle content which is owend by Particles
    TLorentzVector * cparticle = new TLorentzVector(particle->Px(), particle->Py(), particle->Pz(), particle->E());
    container->push_back(cparticle); 
    ParticleNameContainer(ptype)->push_back(name); 
    if (ptype == KLFitter::Particles::kParton) {
      fTrueFlavor->push_back(trueflav); 
      fIsBTagged->push_back(isBtagged); 
      fBTaggingEfficiency->push_back(bTagEff);
      fBTaggingRejection->push_back(bTagRej);
      fJetIndex->push_back(measuredindex);
      fJetDetEta->push_back(DetEta); 
      fBTagWeight->push_back(btagweight);
      if (btagweight != 999) {
	fBTagWeightSet->push_back(true);
	}
	else {
		fBTagWeightSet->push_back(false);
	}
    }
    else if (ptype == KLFitter::Particles::kElectron){ 
      fElectronIndex->push_back(measuredindex);
      fElectronDetEta->push_back(DetEta);
    }
    else if (ptype == KLFitter::Particles::kMuon){ 
      fMuonIndex->push_back(measuredindex);
      fMuonDetEta->push_back(DetEta);
    } 
    else if (ptype == KLFitter::Particles::kPhoton){
      fPhotonIndex->push_back(measuredindex);
      fPhotonDetEta->push_back(DetEta);
    } 
  }
  else {
    std::cout << "KLFitter::Particles::AddParticle(). Particle with the name " << name << " exists already." << std::endl; 
    return 0; 
  }

  if (particle->M() < -1.e-3) {
    std::cout << "KLFitter::Particles::AddParticle(). WARNING : A particle with negative mass " << particle->M() << " of type " << ptype << " was added." << std::endl;
    return 1;
  }

  // no error
  return 1;
}



// --------------------------------------------------------- 
int KLFitter::Particles::AddParticle(TLorentzVector * particle, KLFitter::Particles::ParticleType ptype, std::string name, int measuredindex, bool isBtagged, double bTagEff, double bTagRej, TrueFlavorType trueflav, double btagweight)
{
  //set default DetEta
  double DetEta=-999;
 
  this->AddParticle(particle, DetEta, ptype, name, measuredindex, isBtagged, bTagEff, bTagRej, trueflav, btagweight);
  
  // no error
  return 1;
        
}

// --------------------------------------------------------- 
int KLFitter::Particles::AddParticle(TLorentzVector * particle, KLFitter::Particles::ParticleType ptype, std::string name, int measuredindex, TrueFlavorType trueflav, double btagweight)
{
  //set default DetEta
  double DetEta=-999;
 
  this->AddParticle(particle, DetEta, ptype, name, measuredindex, false, -1., -1., trueflav, btagweight);

  // no error
  return 1;
}

// --------------------------------------------------------- 

int KLFitter::Particles::RemoveParticle(int index, KLFitter::Particles::ParticleType ptype)
{
  // check container and index
  if (!CheckIndex(ParticleContainer(ptype), index))
    return 0; 
        
  // remove particle 
  TLorentzVector* lv = (*ParticleContainer(ptype))[index];
  ParticleContainer(ptype)->erase( ParticleContainer(ptype)->begin() + index); 
  delete lv;
  ParticleNameContainer(ptype)->erase( ParticleNameContainer(ptype)->begin() + index); 

  // no error 
  return 1;
}

// --------------------------------------------------------- 
int KLFitter::Particles::RemoveParticle(std::string name)
{
  // get index and type 
  TLorentzVector* vect = 0; 
  int index = 0; 
  KLFitter::Particles::ParticleType ptype = kParton; 
        
  // remove particle 
  if (FindParticle(name, vect, index, ptype))
    return RemoveParticle(index, ptype); 
  else {
    std::cout << "KLFitter::Particles::RemoveParticles(). Could not find particle with name " << name << "." << std::endl; 
    return 0; 
  }
}


// --------------------------------------------------------- 
TLorentzVector* KLFitter::Particles::Particle(std::string name)
{
  TLorentzVector* particle = 0; 
  int index = 0;
  KLFitter::Particles::ParticleType ptype = kParton; 

  // find particle 
  if (!FindParticle(name, particle, index, ptype)) { 
    std::cout << "KLFitter::Particles::Particle(). Could not find particles." << std::endl; 
    return 0; 
  }

  // return 4-vector 
  return particle; 
}

// --------------------------------------------------------- 
TLorentzVector* KLFitter::Particles::Particle(int index, KLFitter::Particles::ParticleType ptype)
{
  // get particle container
  std::vector <TLorentzVector*>* container = ParticleContainer(ptype);

  if (index < 0 || index > NParticles(ptype)) {
    std::cout << "KLFitter::Particles::Particle(). Index out of range." << std::endl; 
    return 0; 
  }

  // return pointer 
  return (*container)[index];  
}

// --------------------------------------------------------- 
int KLFitter::Particles::FindParticle(std::string name, TLorentzVector* &particle, int &index, KLFitter::Particles::ParticleType &ptype)
{
        
  // loop over all partons 
  unsigned int npartons = fNamePartons->size();
  for (unsigned int i = 0; i < npartons; ++i)
    if (name == (*fNamePartons)[i]) {
      particle = (*fPartons)[i]; 
      index = i; 
      ptype = KLFitter::Particles::kParton; 
      return 1; 
    }

  // loop over all electrons 
  unsigned int nelectrons = fNameElectrons->size();
  for (unsigned int i = 0; i < nelectrons; ++i)
    if (name == (*fNameElectrons)[i]) {
      particle = (*fElectrons)[i]; 
      index = i; 
      ptype = KLFitter::Particles::kElectron; 
      return 1;
    }

  // loop over all muons 
  unsigned int nmuons = fNameMuons->size();
  for (unsigned int i = 0; i < nmuons; ++i)
    if (name == (*fNameMuons)[i]) {
      particle = (*fMuons)[i]; 
      index = i; 
      ptype = KLFitter::Particles::kMuon; 
      return 1;
    }

  // loop over all taus 
  unsigned int ntaus = fNameTaus->size();
  for (unsigned int i = 0; i < ntaus; ++i)
    if (name == (*fNameTaus)[i]) {
      particle = (*fTaus)[i]; 
      index = i; 
      ptype = KLFitter::Particles::kTau; 
      return 1;
    }

  // loop over all neutrinos
  unsigned int nneutrinos = fNameNeutrinos->size(); 
  for (unsigned int i = 0; i < nneutrinos; ++i)
    if (name == (*fNameNeutrinos)[i]) {
      particle = (*fNeutrinos)[i]; 
      index = i; 
      ptype = KLFitter::Particles::kNeutrino; 
      return 1;
    }

  // loop over all bosons 
  unsigned int nbosons = fNameBosons->size();
  for (unsigned int i = 0; i < nbosons; ++i)
    if (name == (*fNameBosons)[i]) {
      particle = (*fBosons)[i]; 
      index = i; 
      ptype = KLFitter::Particles::kBoson; 
      return 1;
    }

  // loop over all photons 
  unsigned int nphotons = fNamePhotons->size();
  for (unsigned int i = 0; i < nphotons; ++i)
    if (name == (*fNamePhotons)[i]) {
      particle = (*fPhotons)[i]; 
      index = i; 
      ptype = KLFitter::Particles::kPhoton; 
      return 1;
    }

  // particle not found
  return 0; 
}

// --------------------------------------------------------- 
TLorentzVector* KLFitter::Particles::Parton(int index)
{
  // no check on index range for CPU-time reasons
  return (*fPartons)[index];

  /*
  // check index 
  if (index < 0 || index >= NPartons())
  {
  std::cout << "KLFitter::Particles::Parton(). Index out of range." << std::endl; 
  return 0; 
  }

  // return pointer 
  return fPartons->at(index); 
  */
}

// --------------------------------------------------------- 
TLorentzVector* KLFitter::Particles::Electron(int index)
{
  // no check on index range for CPU-time reasons
  return (*fElectrons)[index];

  /*
  // check index 
  if (index < 0 || index >= NElectrons())
  {
  std::cout << "KLFitter::Particles::Electron(). Index out of range." << std::endl; 
  return 0; 
  }

  // return pointer 
  return fElectrons->at(index); 
  */
}

// --------------------------------------------------------- 
TLorentzVector* KLFitter::Particles::Muon(int index)
{
  // no check on index range for CPU-time reasons
  return (*fMuons)[index];

  /*
  // check index 
  if (index < 0 || index >= NMuons())
  {
  std::cout << "KLFitter::Particles::Muon(). Index out of range." << std::endl; 
  return 0; 
  }

  // return pointer 
  return fMuons->at(index); 
  */
}

// --------------------------------------------------------- 
TLorentzVector* KLFitter::Particles::Tau(int index)
{
  // no check on index range for CPU-time reasons
  return (*fTaus)[index];

  /*
  // check index 
  if (index < 0 || index >= NTaus())
  {
  std::cout << "KLFitter::Particles::Tau(). Index out of range." << std::endl; 
  return 0; 
  }

  // return pointer 
  return fTaus->at(index); 
  */
}

// --------------------------------------------------------- 
TLorentzVector* KLFitter::Particles::Boson(int index)
{
  // no check on index range for CPU-time reasons
  return (*fBosons)[index];

  /*
  // check index 
  if (index < 0 || index >= NBosons())
  {
  std::cout << "KLFitter::Particles::Boson(). Index out of range." << std::endl; 
  return 0; 
  }

  // return pointer 
  return fBosons->at(index); 
  */

}

// --------------------------------------------------------- 

TLorentzVector* KLFitter::Particles::Neutrino(int index)
{
  // no check on index range for CPU-time reasons
  return (*fNeutrinos)[index];

  /*
  // check index 
  if (index < 0 || index >= NNeutrinos())
  {
  std::cout << "KLFitter::Particles::Neutrino(). Index out of range." << std::endl; 
  return 0; 
  }

  // return pointer 
  return fNeutrinos->at(index); 
  */
}

// --------------------------------------------------------- 
TLorentzVector* KLFitter::Particles::Photon(int index)
{
  // no check on index range for CPU-time reasons
  return (*fPhotons)[index];

  /*
  // check index 
  if (index < 0 || index >= NPhotons())
  {
  std::cout << "KLFitter::Particles::Photon(). Index out of range." << std::endl; 
  return 0; 
  }

  // return pointer 
  return fPhotons->at(index); 
  */
}

// --------------------------------------------------------- 
int  KLFitter::Particles::NParticles(KLFitter::Particles::ParticleType ptype)
{
  return int(ParticleContainer(ptype)->size()); 
}

// --------------------------------------------------------- 
std::string KLFitter::Particles::NameParticle(int index, KLFitter::Particles::ParticleType ptype)
{
  // get particle container
  std::vector <TLorentzVector *>* container = ParticleContainer(ptype); 
        
  // check container and index
  if (!CheckIndex(container, index))
    return ""; 

  // return name 
  return (*ParticleNameContainer(ptype))[index]; 

}

// --------------------------------------------------------- 
int KLFitter::Particles::CheckIndex(std::vector <TLorentzVector *>* container, int index)
{
  // check container
  if (!container)
    {
      std::cout << "KLFitter::Particles::CheckIndex(). Container does not exist." << std::endl; 
      return 0; 
    }

  // check index
  if (index < 0 || index >= int(container->size()))
    {
      std::cout << "KLFitter::Particles::CheckIndex(). Index out of range." << std::endl; 
      return 0; 
    }

  // no error
  return 1; 
}

// --------------------------------------------------------- 
std::vector <TLorentzVector *>* KLFitter::Particles::ParticleContainer(KLFitter::Particles::ParticleType ptype)
{
  // return particle container 
  switch(ptype)
    {
    case KLFitter::Particles::kParton:
      return fPartons;
      break; 
    case KLFitter::Particles::kElectron:
      return fElectrons;
      break; 
    case KLFitter::Particles::kMuon:
      return fMuons;
      break; 
    case KLFitter::Particles::kTau:
      return fTaus;
      break; 
    case KLFitter::Particles::kNeutrino:
      return fNeutrinos;
      break; 
    case KLFitter::Particles::kBoson:
      return fBosons;
      break; 
    case KLFitter::Particles::kPhoton:
      return fPhotons;
      break; 
    }

  // or null pointer 
  std::cout << "KLFitter::Particles::ParticleContainer(). Particle type unknown." << std::endl;
  return 0; 
}

// --------------------------------------------------------- 
std::vector <std::string>* KLFitter::Particles::ParticleNameContainer(KLFitter::Particles::ParticleType ptype)
{
  // return container 
  if (ptype == KLFitter::Particles::kParton)
    return fNamePartons; 
        
  else if (ptype == KLFitter::Particles::kElectron)
    return fNameElectrons; 

  else if (ptype == KLFitter::Particles::kMuon)
    return fNameMuons; 

  else if (ptype == KLFitter::Particles::kTau)
    return fNameTaus; 

  else if (ptype == KLFitter::Particles::kBoson)
    return fNameBosons; 

  else if (ptype == KLFitter::Particles::kNeutrino)
    return fNameNeutrinos; 

  else if (ptype == KLFitter::Particles::kPhoton)
    return fNamePhotons; 

  // or null pointer 
  else
    {
      std::cout << "KLFitter::Particles::ParticleNameContainer(). Particle type not known." << std::endl; 
      return 0; 
    }
}

// --------------------------------------------------------- 
double KLFitter::Particles::DetEta(int index, KLFitter::Particles::ParticleType ptype)
{
  if (index < 0 || index > NParticles(ptype)) {
    std::cout << "KLFitter::Particles::DetEta(). Index out of range." << std::endl; 
    return 0; 
  }

  if (ptype == KLFitter::Particles::kParton)
    return (*fJetDetEta)[index];
  else if (ptype == KLFitter::Particles::kElectron) 
    return (*fElectronDetEta)[index];
  else if (ptype == KLFitter::Particles::kMuon)
    return (*fMuonDetEta)[index];
  else if (ptype == KLFitter::Particles::kPhoton)
    return (*fPhotonDetEta)[index];

  // return error value
  return -100;
}
// --------------------------------------------------------- 
float KLFitter::Particles::LeptonCharge(int index, KLFitter::Particles::ParticleType ptype)
{
  if (index < 0 || index > NParticles(ptype)) {
    std::cout << "KLFitter::Particles::LepCharge(). Index out of range." << std::endl; 
    return 0; 
  }

  
  if (ptype == KLFitter::Particles::kElectron){
    if (fElectronCharge->size()==0)
      return -9;
    else
      return (*fElectronCharge)[index];
  }
  else if (ptype == KLFitter::Particles::kMuon){
    if (fMuonCharge->size()==0)
      return -9;
    else
      return (*fMuonCharge)[index];
  }
  else
    std::cout << "KLFitter::Particles::LepCharge NO LEPTON TYPE!" << std::endl;
  
  // return error value
  return -9;
}
// --------------------------------------------------------- 
int KLFitter::Particles::JetIndex(int index)
{
  // no check on index range for CPU-time reasons
  return (*fJetIndex)[index]; 

  /*
  // check index 
  if (index < 0 || index >= NPartons())
  {
  std::cout << "KLFitter::Particles::JetIndex(). Index out of range." << std::endl; 
  return 0; 
  }

  // return jet index
  return fJetIndex->at(index); 
  */
}

// --------------------------------------------------------- 
int KLFitter::Particles::ElectronIndex(int index)
{
  // no check on index range for CPU-time reasons
  return (*fElectronIndex)[index]; 

  /*
  // check index 
  if (index < 0 || index >= NElectrons())
  {
  std::cout << "KLFitter::Particles::ElectronIndex(). Index out of range." << std::endl; 
  return 0; 
  }

  // return electron index
  return fElectronIndex->at(index); 
  */
}

// --------------------------------------------------------- 
int KLFitter::Particles::MuonIndex(int index)
{
  // no check on index range for CPU-time reasons
  return (*fMuonIndex)[index]; 

  /*
  // check index 
  if (index < 0 || index >= NMuons())
  {
  std::cout << "KLFitter::Particles::MuonIndex(). Index out of range." << std::endl; 
  return 0; 
  }

  // return muon index
  return fMuonIndex->at(index); 
  */
}

// --------------------------------------------------------- 
int KLFitter::Particles::PhotonIndex(int index)
{
  // no check on index range for CPU-time reasons
  return (*fPhotonIndex)[index]; 

  /*
  // check index 
  if (index < 0 || index >= NPhotons())
  {
  std::cout << "KLFitter::Particles::PhotonIndex(). Index out of range." << std::endl; 
  return 0; 
  }

  // return photon index
  return fPhotonIndex->at(index); 
  */
}

// --------------------------------------------------------- 
int KLFitter::Particles::SetIsBTagged(int index, bool isBTagged)
{
  // check index
  if (index < 0 || index >= int(fIsBTagged->size()))
    {
      std::cout << " KLFitter::SetIsBTagged(). Index out of range." << std::endl; 
      return 0; 
    }

  (*fIsBTagged)[index] = isBTagged; 

  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::Particles::SetBTagWeight(int index, double btagweight)
{
  // check index
  if (index < 0 || index >= int(fBTagWeight->size()))
    {
      std::cout << " KLFitter::SetBTagWeight(). Index out of range." << std::endl; 
      return 0; 
    }

  (*fBTagWeight)[index] = btagweight; 
  SetBTagWeightSet(index, true);

  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::Particles::SetBTagWeightSet(int index, bool btagweightset)
{
  // check index
  if (index < 0 || index >= int(fBTagWeightSet->size()))
    {
      std::cout << " KLFitter::SetBTagWeightSet(). Index out of range." << std::endl; 
      return 0; 
    }

  (*fBTagWeightSet)[index] = btagweightset; 

  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::Particles::SetBTaggingEfficiency(int index, double btagEff)
{
  // check index
  if (index < 0 || index >= int(fBTaggingEfficiency->size()))
    {
      std::cout << " KLFitter::SetBTaggingEfficiency(). Index out of range." << std::endl; 
      return 0; 
    }

  (*fBTaggingEfficiency)[index] = btagEff; 

  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::Particles::SetBTaggingRejection(int index, double btagRej)
{
  // check index
  if (index < 0 || index >= int(fBTaggingRejection->size()))
    {
      std::cout << " KLFitter::SetBTaggingRejection(). Index out of range." << std::endl; 
      return 0; 
    }

  (*fBTaggingRejection)[index] = btagRej; 

  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::Particles::NBTags()
{
  unsigned int n = fIsBTagged->size(); 
  int sum = 0; 

  for (unsigned int i = 0; i < n; ++i) {
    if ((*fIsBTagged)[i])
      sum++; 
  }

  return sum; 
}

// --------------------------------------------------------- 
