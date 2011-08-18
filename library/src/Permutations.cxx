#include "Permutations.h" 
#include <iostream> 
#include <set>

// --------------------------------------------------------- 
KLFitter::Permutations::Permutations(KLFitter::Particles ** p, KLFitter::Particles ** pp)
{
  fParticles = p; 
  fParticlesPermuted = pp; 
  fParticlesTable = new std::vector <KLFitter::Particles *>(0); 
  fPermutationTable = new std::vector < std::vector<int> *>(0); 
  fPermutationIndex = -1; 
  fTablePartons = 0; 
  fTableElectrons = 0; 
  fTableMuons = 0; 
  fTablePhotons = 0; 
}

// --------------------------------------------------------- 
KLFitter::Permutations::~Permutations()
{
  // delete particles and permutations tables 
  Reset(); 
}

// --------------------------------------------------------- 
int KLFitter::Permutations::SetPermutation(int index)
{
  // check if permutation table exists 
  if (!fParticlesTable)
    {
      std::cout << "KLFitter::Permutations::SetPermutation(). Table does not exist yet." << std::endl; 
      return 0; 
    }
        
  // check index 
  if (index < 0 || index >= NPermutations())
    {
      std::cout << "KLFitter::Permutations::SetPermutation(). Index out of range." << std::endl; 
      return 0; 
    }

  if (!fParticlesPermuted)
    {
      std::cout << "KLFitter::Permutations::SetPermutation(). Pointer to permuted particles not available." << std::endl; 
      return 0; 
    }

  // set permutation 
  (*fParticlesPermuted) = (*fParticlesTable)[index]; 

  // set permutation index
  fPermutationIndex = index;

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::Permutations::CreatePermutations()
{
  // reset existing particle and permuation tables
  Reset(); 

  // create new table of particles 
  if (!fParticlesTable)
    fParticlesTable = new std::vector <KLFitter::Particles *>(0); 

  // create new table of permutations 
  if (!fPermutationTable)
    fPermutationTable = new std::vector < std::vector<int> *>(0); 

  // check particles
  CheckParticles(); 

  // get number of objects per category 
  int npartons   = (*fParticles)->NPartons(); 
  int nelectrons = (*fParticles)->NElectrons();
  int nmuons     = (*fParticles)->NMuons(); 
  int nphotons     = (*fParticles)->NPhotons(); 

  // create table for parton, electron, muon and photons permutations 
  fTablePartons = new std::vector < std::vector<int> * >(0); 
  CreateSubTable(npartons, fTablePartons); 

  fTableElectrons = new std::vector < std::vector<int> * >(0); 
  CreateSubTable(nelectrons, fTableElectrons); 

  fTableMuons = new std::vector < std::vector<int> * >(0); 
  CreateSubTable(nmuons, fTableMuons); 

  fTablePhotons = new std::vector < std::vector<int> * >(0); 
  CreateSubTable(nphotons, fTablePhotons); 

  // get number of possible permutations for each category 
  int npermpartons   = fTablePartons->size() <= 0 ? 1 : fTablePartons->size(); 
  int npermelectrons = fTableElectrons->size() <= 0 ? 1 : fTableElectrons->size(); 
  int npermmuons     = fTableMuons->size() <= 0 ? 1 : fTableMuons->size(); 
  int npermphotons     = fTablePhotons->size() <= 0 ? 1 : fTablePhotons->size(); 
  int npermoverall   = npartons + nelectrons + nmuons + nphotons;
  
  // loop over all parton permutations
  for (int ipermparton = 0; ipermparton < npermpartons; ++ipermparton)
    {
      // loop over all electron permutations 
      for (int ipermelectron = 0; ipermelectron < npermelectrons; ++ipermelectron)
        {
          // loop over all muon permutations
          for (int ipermmuon = 0; ipermmuon < npermmuons; ++ipermmuon)
            {   
              // loop over all photon permutations
              for (int ipermphoton = 0; ipermphoton < npermphotons; ++ipermphoton)
                {       
                  // create new particles 
                  KLFitter::Particles * particles = new KLFitter::Particles(); 

                  // create new permutation 
                  std::vector <int> * permutation = new std::vector <int>(npermoverall); 

                  // loop over all partons 
                  for (int i = 0; i < npartons; ++i)
                    {
                      // get index 
                      int index = (*(*fTablePartons)[ipermparton])[i];

                      // add parton 
                      particles->AddParticle((*fParticles)->Parton(index),
                                             (*fParticles)->DetEta(index, KLFitter::Particles::kParton),
                                             KLFitter::Particles::kParton,
                                             (*fParticles)->NameParticle(index, KLFitter::Particles::kParton),
                                             (*fParticles)->JetIndex(index),
                                             (*fParticles)->IsBTagged(index),
                                             (*fParticles)->BTaggingEfficiency(index),
                                             (*fParticles)->BTaggingRejection(index),
                                             (*fParticles)->TrueFlavor(index));

                      // set permutation 
                      (*permutation)[i] = index;
                    }

                  // loop over all electrons 
                  for (int i = 0; i < nelectrons; ++i)
                    {
                      // get index 
                      int index = (*(*fTableElectrons)[ipermelectron])[i];

                      // add electron 
                      particles->AddParticle((*fParticles)->Electron(index),
                                             (*fParticles)->DetEta(index, KLFitter::Particles::kElectron),
                                             KLFitter::Particles::kElectron,
                                             (*fParticles)->NameParticle(index, KLFitter::Particles::kElectron),
                                             (*fParticles)->ElectronIndex(index));

                      // set permutation 
                      (*permutation)[npartons + i] = index;
                    }
                                                        
                  // loop over all muons 
                  for (int i = 0; i < nmuons; ++i)
                    {
                      // get index 
                      int index = (*(*fTableMuons)[ipermmuon])[i]; 

                      // add muon 
                      particles->AddParticle((*fParticles)->Muon(index),
                                             (*fParticles)->DetEta(index, KLFitter::Particles::kMuon),
                                             KLFitter::Particles::kMuon,
                                             (*fParticles)->NameParticle(index, KLFitter::Particles::kMuon),
                                             (*fParticles)->MuonIndex(index));

                      // set permutation 
                      (*permutation)[npartons + nelectrons + i] = index; 
                    }

                  // loop over all photons 
                  for (int i = 0; i < nphotons; ++i)
                    {
                      // get index 
                      int index = (*(*fTablePhotons)[ipermphoton])[i]; 

                      // add photon 
                      particles->AddParticle((*fParticles)->Photon(index),
                                             (*fParticles)->DetEta(index, KLFitter::Particles::kPhoton),
                                             KLFitter::Particles::kPhoton,
                                             (*fParticles)->NameParticle(index, KLFitter::Particles::kPhoton),
                                             (*fParticles)->PhotonIndex(index));

                      // set permutation 
                      (*permutation)[npartons + nelectrons + nmuons + i] = index; 
                    }

                  // add particles to table 
                  fParticlesTable->push_back(particles); 
                                                        
                  // add permutation to table 
                  fPermutationTable->push_back(permutation); 
                }
            }
        }
    }                                                   
        
  // no error 
  return 1; 

}

// --------------------------------------------------------- 
int KLFitter::Permutations::Reset()
{
  // delete table of particles 
  if (fParticlesTable) 
    {
      while (!fParticlesTable->empty())
        {
          KLFitter::Particles * p = fParticlesTable->front(); 
          fParticlesTable->erase( fParticlesTable->begin() ); 
          delete p; 
        }                               
      delete fParticlesTable; 
      fParticlesTable = 0;
    }

  // delete permutation table 
  if (fPermutationTable)
    {
      while (!fPermutationTable->empty())
        {
          std::vector <int> * p = fPermutationTable->front(); 
          fPermutationTable->erase( fPermutationTable->begin() ); 
          delete p; 
        }                               
      delete fPermutationTable; 
      fPermutationTable = 0; 
    }

  // free memory 
  if (fTablePartons) {
    while (!fTablePartons->empty())
      {
        std::vector<int> * l = fTablePartons->front(); 
        fTablePartons->erase( fTablePartons->begin() ); 
        l->clear(); 
        delete l; 
      }                                                 
    delete fTablePartons; 
    fTablePartons = 0;
  }

  if (fTableElectrons) {
    while (!fTableElectrons->empty())
      {
        std::vector<int> * l = fTableElectrons->front(); 
        fTableElectrons->erase( fTableElectrons->begin() ); 
        l->clear(); 
        delete l; 
      }                                                 
    delete fTableElectrons; 
    fTableElectrons = 0;
  }

  if (fTableMuons) {
    while (!fTableMuons->empty())
      {
        std::vector<int> * l = fTableMuons->front(); 
        fTableMuons->erase( fTableMuons->begin() ); 
        l->clear(); 
        delete l; 
      }                                                 
    delete fTableMuons;
    fTableMuons = 0; 
  }

  if (fTablePhotons) {
    while (!fTablePhotons->empty())
      {
        std::vector<int> * l = fTablePhotons->front(); 
        fTablePhotons->erase( fTablePhotons->begin() ); 
        l->clear(); 
        delete l; 
      }                                                 
    delete fTablePhotons;
    fTablePhotons = 0; 
  }

  // no error 
  return 1; 
}


// --------------------------------------------------------- 
int KLFitter::Permutations::CreateSubTable(int Nobj, std::vector < std::vector<int> * > * table)
{
  // initialize permutations (original and temp) 
  std::vector<int> op(Nobj, -1); 
  std::vector<int> tp(Nobj+1, -1); 

  // calculate number of permutations
  int Nperm = 1; 
  for (int i = Nobj; i > 0; --i)
    Nperm *= i; 
        
  // fill permutations
  for (int i = 0; i < Nobj; ++i)
    {
      op[i] = i; 
      tp[i] = i; 
    }
  tp[Nobj] = Nobj; 

  // add original permutation to table 
  std::vector<int> * tempp = new std::vector<int>(0); 
  *tempp = op; 
  table->push_back(tempp); 

  int i, j, temp; 
  i = 1; 
  while (i < Nobj)
    {
      (tp[i])--; 
      j = (Nobj-1) - i % 2 * tp[i]; 
      i = (Nobj-1) - i; 

      // swap elements i and j 
      temp = op[j]; 
      op[j] = op[i]; 
      op[i] = temp; 

      // add new permutation 
      tempp = new std::vector<int>(Nobj, -1);
      for (int k = 0; k < Nobj; ++k)
        (*tempp)[k] = op[k]; 
      table->push_back(tempp); 

      i = 1; 
      while (!tp[i])
        {
          tp[i] = i; 
          i++;
        }
    }

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::Permutations::InvariantParticlePermutations(KLFitter::Particles::ParticleType ptype, std::vector<int> indexVector)
{
  // check if particles are defined
  if (!CheckParticles())
    return 0; 

  // check indices 
  std::set<int> indexSet;
  for (unsigned int i = 0, I = indexVector.size(); i < I; i++)
    {
      if (indexSet.find(indexVector[i]) != indexSet.end())
        {
          std::cout << "KLFitter::Permutations::InvariantParticlePermutations(). Indices have to be different." << std::endl; 
          return 0; 
        }
      else
        indexSet.insert(indexVector[i]);
    }

  for (unsigned int i = 0, I = indexVector.size(); i < I; i++)
    {
      int index = indexVector[i];
      if (index < 0 || index >= (*fParticles)->NParticles(ptype))
        {
          std::cout << "KLFitter::Permutations::InvariantParticlePermutations(). Index out of range." << std::endl; 
          return 0; 
        }
    }

  // swap indices 
  indexVector.clear();
  std::set<int>::iterator it_indexSetBegin = indexSet.begin();
  std::set<int>::iterator it_indexSetEnd   = indexSet.end();
  for (; it_indexSetBegin != it_indexSetEnd; it_indexSetBegin++)
    indexVector.push_back(*it_indexSetBegin);

  // check particles table 
  if (!fParticlesTable)
    {
      std::cout << "KLFitter::Permutations::InvariantParticlePermutations(). Table does not exist yet." << std::endl; 
      return 0; 
    }
 
  // check permutation table
  if (!fPermutationTable)
    {
      std::cout << "KLFitter::Permutations::InvariantParticlePermutations(). Table of parton permutations doesn ot exist." << std::endl; 
      return 0; 
    }
 
  // no error 
  int err = 1; 

  // loop over all permutations (if there are only 2 indices left)
  if (indexVector.size() == 2)
    {
      // get number of permutations 
      int nperm = NPermutations(); 

      for (int iperm = nperm-1; iperm >= 0; --iperm)
        {
          int offset = 0;
          for (KLFitter::Particles::ParticleType itype = KLFitter::Particles::kParton; itype < ptype; ++itype)
            offset += (*fParticles)->NParticles(itype);
          int index1 = indexVector[0] + offset;
          int index2 = indexVector[1] + offset;

          // get permutation 
          std::vector<int> * permutation = (*fPermutationTable)[iperm]; 
                        
          // check indices 
          if ((*permutation)[index1] >= (*permutation)[index2]) 
            {
              std::vector <int> * pt = (*fPermutationTable)[iperm]; 
              fPermutationTable->erase( fPermutationTable->begin() + iperm ); 
              delete pt; 

              KLFitter::Particles * p = (*fParticlesTable)[iperm]; 
              fParticlesTable->erase( fParticlesTable->begin() + iperm ); 
              delete p; 
            }
        }
    }
  else
    {
      // repeat until there are only 2 indices left
      while (indexVector.size() >= 2)
        {
          int index2 = indexVector.back();
          for (unsigned int i = 0, I = indexVector.size()-1; i < I; i++)
            {
              int index1 = indexVector[i];
              std::vector<int> newIndexVector;
              newIndexVector.push_back(index1);
              newIndexVector.push_back(index2);
              err *= InvariantParticlePermutations(ptype, newIndexVector);
            }
          indexVector.erase(--indexVector.end());
        }
    }

  // return error code 
  return err; 
}

// --------------------------------------------------------- 
int KLFitter::Permutations::CheckParticles()
{
  // check if particles are defined
  if (!fParticles)
    {
      std::cout << "KLFitter::Permutations::CheckParticles(). Pointer to particles not defined yet." << std::endl; 
      return 0; 
    }

  if (!(*fParticles))
    {
      std::cout << "KLFitter::Permutations::CheckParticles(). Particles not defined yet." << std::endl; 
      return 0; 
    }
        
  // no error 
  return 1;
}

// --------------------------------------------------------- 

std::vector<std::vector<int>* > *KLFitter::Permutations::PermutationTable()
{
  return fPermutationTable;

  /*
  // copy the table
  std::vector<std::vector<int>* > *permutationTable = new std::vector<std::vector<int>* >(0);
  int npartons   = (*fParticles)->NPartons(); 
  int nelectrons = (*fParticles)->NElectrons(); 
  int nmuons     = (*fParticles)->NMuons(); 
  int nphotons     = (*fParticles)->NPhotons(); 
  for (unsigned int iPerm = 0; iPerm < fPermutationTable->size(); iPerm++)
  {
  std::vector<int> *iVec = new std::vector<int>;
  for (int iParticle = 0; iParticle < npartons + nelectrons + nmuons + nphotons; iParticle++)
  iVec->push_back(fPermutationTable->at(iPerm)->at(iParticle));
  permutationTable->push_back(iVec);
  }

  // return copy
  return permutationTable;
  */
}

// --------------------------------------------------------- 

