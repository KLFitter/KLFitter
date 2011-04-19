#include "MatchingTool.h" 

#include <iostream> 


// --------------------------------------------------------- 
KLFitter::MatchingTool::MatchingTool(KLFitter::Particles ** particles, KLFitter::Particles ** particlestruth)
{
  fParticles = particles;
  fParticlesTruth = particlestruth;
  fDeltaRPartons = 0.3; 
  fDeltaRElectrons = 0.1; 
  fDeltaRMuons = 0.1; 
  fDeltaRPhotons = 0.1; 
  fMatchedPartons = new std::vector< std::vector<int> >(0);
  fMatchedElectrons = new std::vector< std::vector<int> >(0);
  fMatchedMuons = new std::vector< std::vector<int> >(0);
  fMatchedPhotons = new std::vector< std::vector<int> >(0);
}

// --------------------------------------------------------- 
KLFitter::MatchingTool::~MatchingTool()
{
  if (fMatchedPartons)
    delete fMatchedPartons;
        
  if (fMatchedElectrons)
    delete fMatchedElectrons;
        
  if (fMatchedMuons)
    delete fMatchedMuons;

  if (fMatchedPhotons)
    delete fMatchedPhotons;
}

// --------------------------------------------------------- 
void KLFitter::MatchingTool::SetDeltaRPartons(double dR)
{
  if (dR < 0)
    dR = 0.3; 
  else
    fDeltaRPartons = dR; 
}

// --------------------------------------------------------- 
void KLFitter::MatchingTool::SetDeltaRElectrons(double dR)
{
  if (dR < 0)
    dR = 0.1; 
  else
    fDeltaRElectrons = dR; 
}

// --------------------------------------------------------- 
void KLFitter::MatchingTool::SetDeltaRMuons(double dR)
{
  if (dR < 0)
    dR = 0.1; 
  else
    fDeltaRMuons = dR; 
}

// --------------------------------------------------------- 
void KLFitter::MatchingTool::SetDeltaRPhotons(double dR)
{
  if (dR < 0)
    dR = 0.1; 
  else
    fDeltaRPhotons = dR; 
}

// ---------------------------------------------------------
int KLFitter::MatchingTool::Initialize(KLFitter::Particles::ParticleType ptype)
{
  if (!fParticlesTruth)
    {
      std::cout << "KLFitter::MatchingTool::Initialize(). Pointer to truth particle container not defined." << std::endl; 
      return 0; 
    }

  if (!(*fParticlesTruth))
    {
      std::cout << "KLFitter::MatchingTool::Initialize(). Truth particle container not defined." << std::endl; 
      return 0; 
    }

  if (!fParticles)
    {
      std::cout << "KLFitter::MatchingTool::Initialize(). Pointer to particle container not defined." << std::endl; 
      return 0; 
    }

  if (!(*fParticles))
    {
      std::cout << "KLFitter::MatchingTool::Initialize(). Particle container not defined." << std::endl; 
      return 0; 
    }

  int n_truth = -1; 
  int n_reco = -1; 
  if (ptype == KLFitter::Particles::kParton)
    {
      // get number of particles 
      n_truth = (*fParticlesTruth) -> NPartons(); 
      n_reco = (*fParticles) -> NPartons(); 
                        
      // clear matrix with matching status
      fMatchedPartons -> clear(); 

      // initialze matrix with status of matching
      for(int j=0; j<n_truth; ++j)
        {
          fMatchedPartons -> push_back(std::vector<int>(0));
          for(int i=0; i<n_reco; ++i)
            fMatchedPartons -> at(j).push_back(-1);
        }
    }

  else if (ptype == KLFitter::Particles::kElectron)
    {
      // get number of particles 
      n_truth = (*fParticlesTruth) -> NElectrons(); 
      n_reco = (*fParticles) -> NElectrons(); 
                        
      // clear matrix with matching status
      fMatchedElectrons -> clear(); 

      // initialze matrix with status of matching
      for(int j=0; j<n_truth; ++j)
        {
          fMatchedElectrons -> push_back(std::vector<int>(0));
          for(int i=0; i<n_reco; ++i)
            fMatchedElectrons -> at(j).push_back(-1);
        }
    }

  else if (ptype == KLFitter::Particles::kMuon)
    {
      // get number of particles 
      n_truth = (*fParticlesTruth) -> NMuons(); 
      n_reco = (*fParticles) -> NMuons(); 
                        
      // clear matrix with matching status
      fMatchedMuons -> clear(); 

      // initialze matrix with status of matching
      for(int j=0; j<n_truth; ++j)
        {
          fMatchedMuons -> push_back(std::vector<int>(0));
          for(int i=0; i<n_reco; ++i)
            fMatchedMuons -> at(j).push_back(-1);
        }
    }

  else if (ptype == KLFitter::Particles::kPhoton)
    {
      // get number of particles 
      n_truth = (*fParticlesTruth) -> NPhotons(); 
      n_reco = (*fParticles) -> NPhotons(); 
                        
      // clear matrix with matching status
      fMatchedPhotons -> clear(); 

      // initialze matrix with status of matching
      for(int j=0; j<n_truth; ++j)
        {
          fMatchedPhotons -> push_back(std::vector<int>(0));
          for(int i=0; i<n_reco; ++i)
            fMatchedPhotons -> at(j).push_back(-1);
        }
    }

  else
    {
      std::cout << "KLFitter::MatchingTool::Initialize(). Unknown particle type." << std::endl; 
      return 0; 
    }

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::MatchingTool::MatchVectors(int index1, int index2, KLFitter::Particles::ParticleType ptype)
{
  // get Lorentz vectors
  TLorentzVector * vect1 = (*fParticles) -> Particle(index1, ptype);
  TLorentzVector * vect2 = (*fParticlesTruth) -> Particle(index2, ptype);

  double dR = 0.; 
  if (ptype == KLFitter::Particles::kParton)
    dR = fDeltaRPartons; 
  else if (ptype == KLFitter::Particles::kElectron)
    dR = fDeltaRElectrons; 
  else if (ptype == KLFitter::Particles::kMuon)
    dR = fDeltaRMuons; 
  else if (ptype == KLFitter::Particles::kPhoton)
    dR = fDeltaRPhotons; 

  // check if Lorentz vectors exist 
  if (vect1 == 0 || vect2 == 0)
    {
      std::cout << "KLFitter::MatchingTool::MatchVectors(). Could not find Lorentz vectors." << std::endl; 
      return -1; 
    }
  else
    return MatchDeltaR(vect1, vect2, dR); 
}

// --------------------------------------------------------- 
int KLFitter::MatchingTool::MatchDeltaR(TLorentzVector * vect1, TLorentzVector * vect2, double dr)
{
  // get deltaR
  double MaxDeltaR = dr;
  
  // get eta and phi for the objects
  double recophi = 0.;
  double recoeta = 0.;
  recophi = vect1->Phi();
  recoeta = vect1->Eta();
  double truephi = 0.;
  double trueeta = 0.;
  truephi = vect2->Phi();
  trueeta = vect2->Eta();
  
  // calculate deltaR and make sure about boundaries
  double dEta = recoeta - trueeta;
  double dPhi = recophi - truephi; 
  const double dPi = 3.1415926535897932384626433832795;
  while (dPhi > dPi) {dPhi -= 2*dPi;}
  while (dPhi < -dPi) {dPhi += 2*dPi;}
  
  double dR = sqrt(dEta*dEta + dPhi*dPhi);  

  // check if deltaR is within our range
  if (dR < MaxDeltaR)
    return 1;
  //else
  return 0;
}

// --------------------------------------------------------- 
int KLFitter::MatchingTool::MatchTruth(int index, KLFitter::Particles::ParticleType ptype)
{
  // get number of reco particles
  int nreco = (*fParticles) -> NParticles(ptype); 

  // loop over all reco particles
  for (int i = 0; i < nreco; ++i)
    {
      // match truth and reco particle
      int stat = this -> MatchVectors(i, index, ptype); 

      // write into list of matching status 
      this -> SetMatchingStatus(index, i, ptype, stat); 
    }

  // no error
  return 1;
}

// --------------------------------------------------------- 
int KLFitter::MatchingTool::MatchTruthAll(KLFitter::Particles::ParticleType ptype)
{
  // initialize
  this -> Initialize(ptype); 

  // get number of reco particles
  int ntruth = (*fParticlesTruth) -> NParticles(ptype); 

  // loop over all truth particles
  for (int i = 0; i < ntruth; ++i)
    {
      // match truth and reco particle
      if (!this -> MatchTruth(i, ptype))
        return 0; 
    }

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::MatchingTool::SetMatchingStatus(int indextruth, int indexreco, KLFitter::Particles::ParticleType ptype, int stat)
{
  if (ptype == KLFitter::Particles::kParton)
    {
      if (indextruth >= 0 && indextruth < (*fParticlesTruth) -> NPartons() && 
          indexreco >= 0 && indexreco < (*fParticles) -> NPartons())
        fMatchedPartons -> at(indextruth).at(indexreco) = stat; 
      else
        {
          std::cout << "KLFitter::MatchingTool::SetMatchingStatus(). Parton index out of range." << std::endl; 
          return 0; 
        }
    }
  else if (ptype == KLFitter::Particles::kElectron)
    {
      if (indextruth >= 0 && indextruth < (*fParticlesTruth) -> NElectrons() && 
          indexreco >= 0 && indexreco < (*fParticles) -> NElectrons())
        fMatchedElectrons -> at(indextruth).at(indexreco) = stat; 
      else
        {
          std::cout << "KLFitter::MatchingTool::SetMatchingStatus(). Electron index out of range." << std::endl; 
          return 0; 
        }
    }
  else if (ptype == KLFitter::Particles::kMuon)
    {
      if (indextruth >= 0 && indextruth < (*fParticlesTruth) -> NMuons() && 
          indexreco >= 0 && indexreco < (*fParticles) -> NMuons())
        fMatchedMuons -> at(indextruth).at(indexreco) = stat; 
      else
        {
          std::cout << "KLFitter::MatchingTool::SetMatchingStatus(). Muon index out of range." << std::endl; 
          return 0; 
        }
    }
  else if (ptype == KLFitter::Particles::kPhoton)
    {
      if (indextruth >= 0 && indextruth < (*fParticlesTruth) -> NPhotons() && 
          indexreco >= 0 && indexreco < (*fParticles) -> NPhotons())
        fMatchedPhotons -> at(indextruth).at(indexreco) = stat; 
      else
        {
          std::cout << "KLFitter::MatchingTool::SetMatchingStatus(). Photon index out of range." << std::endl; 
          return 0; 
        }
    }
  else
    {
      std::cout << "KLFitter::MatchingTool::SetMatchingStatus(). Unknown particle type." << std::endl; 
      return 0; 
    }

  // no error 
  return 1;
}

// --------------------------------------------------------- 
int KLFitter::MatchingTool::MatchingStatus(int indextruth, int indexreco, KLFitter::Particles::ParticleType ptype)
{
  if (ptype == KLFitter::Particles::kParton)
    {
      if (indextruth >= 0 && indextruth < (*fParticlesTruth) -> NPartons() && 
          indexreco >= 0 && indexreco < (*fParticles) -> NPartons())
        return fMatchedPartons -> at(indextruth).at(indexreco); 
      else
        {
          std::cout << "KLFitter::MatchingTool::IsMatched(). Parton index out of range." << std::endl; 
          return 0; 
        }
    }
  else if (ptype == KLFitter::Particles::kElectron)
    {
      if (indextruth >= 0 && indextruth < (*fParticlesTruth) -> NElectrons() && 
          indexreco >= 0 && indexreco < (*fParticles) -> NElectrons())
        return fMatchedElectrons -> at(indextruth).at(indexreco); 
      else
        {
          std::cout << "KLFitter::MatchingTool::IsMatched(). Electron index out of range." << std::endl; 
          return 0; 
        }
    }
  else if (ptype == KLFitter::Particles::kMuon)
    {
      if (indextruth >= 0 && indextruth < (*fParticlesTruth) -> NMuons() && 
          indexreco >= 0 && indexreco < (*fParticles) -> NMuons())
        return fMatchedMuons -> at(indextruth).at(indexreco); 
      else
        {
          std::cout << "KLFitter::MatchingTool::IsMatched(). Muon index out of range." << std::endl; 
          return 0; 
        }
    }
  else if (ptype == KLFitter::Particles::kPhoton)
    {
      if (indextruth >= 0 && indextruth < (*fParticlesTruth) -> NPhotons() && 
          indexreco >= 0 && indexreco < (*fParticles) -> NPhotons())
        return fMatchedPhotons -> at(indextruth).at(indexreco); 
      else
        {
          std::cout << "KLFitter::MatchingTool::IsMatched(). Photon index out of range." << std::endl; 
          return 0; 
        }
    }
  else
    {
      std::cout << "KLFitter::MatchingTool::IsMatched(). Unknown particle type." << std::endl; 
      return 0; 
    }

  // no error 
  return 1;
}

// --------------------------------------------------------- 
int KLFitter::MatchingTool::NMatchedTruth(int index,  KLFitter::Particles::ParticleType ptype)
{
  int nmatches = 0;
  int njets = (*fParticles) -> NPartons();   
  int nelectrons = (*fParticles) -> NElectrons(); 
  int nmuons = (*fParticles) -> NMuons(); 
  int nphotons = (*fParticles) -> NPhotons(); 

  int npartonstruth = (*fParticlesTruth) -> NPartons(); 
  int nelectronstruth = (*fParticlesTruth) -> NElectrons(); 
  int nmuonstruth = (*fParticlesTruth) -> NMuons(); 
  int nphotonstruth = (*fParticlesTruth) -> NPhotons(); 

  // return size of our vector to check number of matches - parton
  if(ptype == KLFitter::Particles::kParton)
    { 
      if (index < 0 || index >= npartonstruth)
        return -1;
                        
      for(int i = 0; i < njets; ++i)
        if ((fMatchedPartons-> at(index)).at(i) > 0)
          nmatches++;
                        
      // return number of matches
      return nmatches; 
    }

  // return size of our vector to check number of matches - electrons
  if(ptype == KLFitter::Particles::kElectron)
    { 
      if (index < 0 || index >= nelectronstruth)
        return -1;

      for(int i = 0; i < nelectrons; ++i)
        if ((fMatchedElectrons-> at(index)).at(i) > 0)
          nmatches++;

      // return number of matches
      return nmatches; 
    }

  // return size of our vector to check number of matches - muons
  if(ptype == KLFitter::Particles::kMuon)
    { 
      if (index < 0 || index >= nmuonstruth)
        return -1;

      for(int i = 0; i < nmuons; ++i)
        if ((fMatchedMuons-> at(index)).at(i) > 0)
          nmatches++;

      // return number of matches
      return nmatches; 
    }

  // return size of our vector to check number of matches - photons
  if(ptype == KLFitter::Particles::kPhoton)
    { 

      if (index < 0 || index >= nphotonstruth)
        return -1;

      for(int i = 0; i < nphotons; ++i) {
        if ((fMatchedPhotons-> at(index)).at(i) > 0)
          nmatches++;

        // return number of matches
        return nmatches; 
      }
    }
  else
    {
      std::cout << "KLFitter::MatchingTool::NMatchedTruth(): Particle type unknown." << std::endl;
      return -1;
    }
  return -1;
}

// --------------------------------------------------------- 
// int KLFitter::MatchingTool::NMatchedTruth(const char * name)
// {
//   KLFitter::Particles::ParticleType ptype;
//   int index;

//   // find the particle according to the name
//   fParticlesTruth -> FindParticle(name, index, ptype);

//   // run the function NMatchedTruth to find the number of matches
//   int size = NMatchedTruth(index, ptype); 
  
//      // return the number of reco particles matched to the truth particle
//   return size;
// }

// --------------------------------------------------------- 
std::vector<int> KLFitter::MatchingTool::ListMatchedTruth(int index,  KLFitter::Particles::ParticleType ptype)
{
  int njets = (*fParticlesTruth) -> NPartons();   
  int neles = (*fParticlesTruth) -> NElectrons();
  int nmuons = (*fParticlesTruth) -> NMuons();
  int nphotons = (*fParticlesTruth) -> NPhotons();
  // check type
  if(ptype == KLFitter::Particles::kParton){ 
    // check for out of range index
    if(index > njets-1 ){
      std::vector<int> errvec[1];
      errvec->push_back(-1);
      return *(errvec);   
    } 
    // return vector linking truth particle
    return fMatchedPartons->at(index);
  }
  else if(ptype == KLFitter::Particles::kElectron){
    // check for out of range index
    if(index > neles-1){
      std::vector<int> errvec[1];
      errvec->push_back(-1);
      return *(errvec);  
    }
    // return vector linking truth particle
    return fMatchedElectrons->at(index);
  }
  else if(ptype == KLFitter::Particles::kMuon){
    // check for out of range index
    if(index > nmuons-1){
      std::vector<int> errvec[1];
      errvec->push_back(-1);
      return *(errvec);  
    }
    // return vector linking truth particle
    return fMatchedMuons->at(index);
  }
  else if(ptype == KLFitter::Particles::kPhoton){
    // check for out of range index
    if(index > nphotons-1){
      std::vector<int> errvec[1];
      errvec->push_back(-1);
      return *(errvec);  
    }
    // return vector linking truth particle
    return fMatchedPhotons->at(index);
  }  
  else{
    // error - particle type incorrect - return vector with -1
    std::vector<int> errvec[1];
    errvec->push_back(-1);
    return *(errvec);
  }
}

// --------------------------------------------------------- 
// std::vector<int> KLFitter::MatchingTool::ListMatchedTruth(const char * name)
// {
//   KLFitter::Particles::ParticleType ptype;
//   int index;
//   // find the particle according to the name  
//   fParticlesTruth -> FindParticle(name ,index, ptype);
//   // run the function ListMatchedTruth to get the vector on matches
//   std::vector<int> outvec = ListMatchedTruth(index,ptype);

//   return outvec;
// }
