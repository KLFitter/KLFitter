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

  fBTaggingProbability(new std::vector <double>(0)),
  fFlavorTag(new std::vector<double>(0)),

  fElectronDetEta(new std::vector<double>(0)),
  fMuonDetEta(new std::vector<double>(0)),
  fJetDetEta(new std::vector<double>(0)),
  fPhotonDetEta(new std::vector<double>(0))
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

  if (fBTaggingProbability)
    delete fBTaggingProbability; 

  if (fFlavorTag)
    delete fFlavorTag;

  if (fElectronDetEta)
    delete fElectronDetEta;

  if (fMuonDetEta)
    delete fMuonDetEta;

  if (fJetDetEta)
    delete fJetDetEta;

  if (fPhotonDetEta)
    delete fPhotonDetEta;
}
// --------------------------------------------------------- 
int KLFitter::Particles::AddParticle(TLorentzVector* particle, KLFitter::Particles::ParticleType ptype, std::string name, double btagprob, double flavortag, int measuredindex)
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
      fBTaggingProbability->push_back(btagprob); 
      fFlavorTag->push_back(flavortag); 
      fJetIndex->push_back(measuredindex);      
    }
    else if (ptype == KLFitter::Particles::kElectron){ 
      fElectronIndex->push_back(measuredindex);
    }
    else if (ptype == KLFitter::Particles::kMuon){ 
      fMuonIndex->push_back(measuredindex);
     } 
    else if (ptype == KLFitter::Particles::kPhoton){
      fPhotonIndex->push_back(measuredindex);
     } 
  }
  else {
    std::cout << "KLFitter::Particles::AddParticle(). Particle with the name " << name << " exists already." << std::endl; 
    return 0; 
  }

  // no error
  return 1;
        
}

// --------------------------------------------------------- 
int KLFitter::Particles::AddParticle(TLorentzVector* particle, double DetEta, KLFitter::Particles::ParticleType ptype, std::string name, double btagprob, double flavortag, int measuredindex)
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
      fBTaggingProbability->push_back(btagprob); 
      fFlavorTag->push_back(flavortag); 
      fJetIndex->push_back(measuredindex);
      fJetDetEta->push_back(DetEta); 
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
  TLorentzVector* lv = ParticleContainer(ptype)->at(index);
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
  return container->at(index);  
}

// --------------------------------------------------------- 
int KLFitter::Particles::FindParticle(std::string name, TLorentzVector* &particle, int &index, KLFitter::Particles::ParticleType &ptype)
{
        
  // loop over all partons 
  unsigned int npartons = fNamePartons->size();
  for (unsigned int i = 0; i < npartons; ++i)
    if (name == fNamePartons->at(i)) {
      particle = fPartons->at(i); 
      index = i; 
      ptype = KLFitter::Particles::kParton; 
      return 1; 
    }

  // loop over all electrons 
  unsigned int nelectrons = fNameElectrons->size();
  for (unsigned int i = 0; i < nelectrons; ++i)
    if (name == fNameElectrons->at(i)) {
      particle = fElectrons->at(i); 
      index = i; 
      ptype = KLFitter::Particles::kElectron; 
      return 1;
    }

  // loop over all muons 
  unsigned int nmuons = fNameMuons->size();
  for (unsigned int i = 0; i < nmuons; ++i)
    if (name == fNameMuons->at(i)) {
      particle = fMuons->at(i); 
      index = i; 
      ptype = KLFitter::Particles::kMuon; 
      return 1;
    }

  // loop over all taus 
  unsigned int ntaus = fNameTaus->size();
  for (unsigned int i = 0; i < ntaus; ++i)
    if (name == fNameTaus->at(i)) {
      particle = fTaus->at(i); 
      index = i; 
      ptype = KLFitter::Particles::kTau; 
      return 1;
    }

  // loop over all neutrinos
  unsigned int nneutrinos = fNameNeutrinos->size(); 
  for (unsigned int i = 0; i < nneutrinos; ++i)
    if (name == fNameNeutrinos->at(i)) {
      particle = fNeutrinos->at(i); 
      index = i; 
      ptype = KLFitter::Particles::kNeutrino; 
      return 1;
    }

  // loop over all bosons 
  unsigned int nbosons = fNameBosons->size();
  for (unsigned int i = 0; i < nbosons; ++i)
    if (name == fNameBosons->at(i)) {
      particle = fBosons->at(i); 
      index = i; 
      ptype = KLFitter::Particles::kBoson; 
      return 1;
    }

  // loop over all photons 
  unsigned int nphotons = fNamePhotons->size();
  for (unsigned int i = 0; i < nphotons; ++i)
    if (name == fNamePhotons->at(i)) {
      particle = fPhotons->at(i); 
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
  return fPartons->at(index);

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
  return fElectrons->at(index);

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
  return fMuons->at(index);

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
  return fTaus->at(index);

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
  return fBosons->at(index);

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
  return fNeutrinos->at(index);

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
  return fPhotons->at(index);

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
  return ParticleNameContainer(ptype)->at(index); 

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
double KLFitter::Particles::BTaggingProbability(int index)
{
  // no check on index range for CPU-time reasons
  return fBTaggingProbability->at(index); 

  /*
  // check index 
  if (index < 0 || index >= NPartons())
  {
  std::cout << "KLFitter::Particles::BTaggingProbability(). Index out of range." << std::endl; 
  return 0; 
  }

  // return b-tagging probability
  return fBTaggingProbability->at(index); 
  */
}

// --------------------------------------------------------- 
double KLFitter::Particles::FlavorTag(int index)
{
  // no check on index range for CPU-time reasons
  return fFlavorTag->at(index); 

  /*
  // check index 
  if (index < 0 || index >= NPartons())
  {
  std::cout << "KLFitter::Particles::FlavorTag(). Index out of range." << std::endl; 
  return 0; 
  }

  // return b-tagging probability
  return fFlavorTag->at(index); 
  */
}
// --------------------------------------------------------- 
double KLFitter::Particles::DetEta(int index, KLFitter::Particles::ParticleType ptype)
{
  if (index < 0 || index > NParticles(ptype)) {
    std::cout << "KLFitter::Particles::DetEta(). Index out of range." << std::endl; 
    return 0; 
  }

  if (ptype == KLFitter::Particles::kParton)
    return fJetDetEta->at(index);
  else if (ptype == KLFitter::Particles::kElectron) 
    return fElectronDetEta->at(index);
  else if (ptype == KLFitter::Particles::kMuon)
    return fMuonDetEta->at(index);
  else if (ptype == KLFitter::Particles::kPhoton)
    return fPhotonDetEta->at(index);

  // return error value
  return -100;
}
// --------------------------------------------------------- 
int KLFitter::Particles::JetIndex(int index)
{
  // no check on index range for CPU-time reasons
  return fJetIndex->at(index); 

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
  return fElectronIndex->at(index); 

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
  return fMuonIndex->at(index); 

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
  return fPhotonIndex->at(index); 

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
int KLFitter::Particles::SetFlavorTag(int index, double tag)
{
  // check index
  if (index < 0 || index >= int(fFlavorTag->size()))
    {
      std::cout << " KLFitter::Combinatorics::SetFlavorTag(). Index out of range." << std::endl; 
      return 0; 
    }

  (*fFlavorTag)[index] = tag; 

  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::Particles::SetBTaggingProbability(int index, double prob)
{
  // check index
  if (index < 0 || index >= int(fFlavorTag->size()))
    {
      std::cout << " KLFitter::Combinatorics::SetBTaggingProbability(). Index out of range." << std::endl; 
      return 0; 
    }

  (*fBTaggingProbability)[index] = prob; 

  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::Particles::NBTags()
{
  unsigned int n = fFlavorTag->size(); 
  int sum = 0; 

  for (unsigned int i = 0; i < n; ++i) {
    if (fFlavorTag->at(i) == 1)
      sum++; 
  }

  return sum; 
}

// --------------------------------------------------------- 
