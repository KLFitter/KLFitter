#include "SelectionTool.h" 
#include <iostream> 

// --------------------------------------------------------- 
KLFitter::SelectionTool::SelectionTool() :
  fParticlesSelected(0),
  fJetPt(0.0), 
  fJetEta(2.5),
  fJetJVF(0.75),
  fElectronPt(0.0),
  fElectronEta(2.5), 
  fMuonPt(0.0),
  fMuonEta(2.5),
  fPhotonPt(0.0),
  fPhotonEta(2.5),
  fMET(0), 
  fMWT(0), 
	fMET_plus_MWT(0),
  fCounterEvents(0),  
  fCounterJets(0),  
  fCounterBJets(0),  
  fCounterElectrons(0),  
  fCounterMuons(0),  
  fCounterPhotons(0),  
  fCounterMET(0),
  fCounterMWT(0),
  fCounterTriangular(0),
  fCounterSelected(0),
  fMaxNJetsForFit(4)
{
}

// --------------------------------------------------------- 
KLFitter::SelectionTool::~SelectionTool()
{
  if (fParticlesSelected)
    delete fParticlesSelected; 
}

// --------------------------------------------------------- 
int KLFitter::SelectionTool::SelectObjects(KLFitter::Particles * particles)
{
  // check if particles exist
  if (!particles)
    {
      std::cout << "KLFitter::SelectionTool::SelectObjects. Particles do not exist." << std::endl; 
      return 0; 
    }
  
  // clear events
  if (fParticlesSelected) 
    delete fParticlesSelected; 

  bool isDilepton(false);

  fParticlesSelected = new KLFitter::Particles(); 
  
  // reset maps 
  this->ResetMaps(); 
  
  // jet selection
  int npartons = particles->NPartons(); 
  
  // debugKK
  // re-order such that b-tagged jets come first
  
  if (fNBJets.size() > 0) {
    // run over b-jets first
    for (int i = 0; i < npartons; ++i) {
      // check b-tag
      if (particles->BTagWeight(i) < fNBJets.at(0).value) 
	continue;
      
      // check eta region 
      if (TMath::Abs( particles->DetEta(i, KLFitter::Particles::kParton) ) > fJetEta) 
        continue; 
      
      // check pT 
      if (particles->Parton(i)->Pt() < fJetPt)
        continue; 
      
      // add jet 
      TLorentzVector * tlv_tmp = new TLorentzVector(*particles->Parton(i));
      fParticlesSelected->AddParticle( tlv_tmp,
				       //				       new TLorentzVector(*particles->Parton(i)),
                                       particles->DetEta(i, KLFitter::Particles::kParton),
                                       KLFitter::Particles::kParton, 
                                       particles->NameParticle(i, KLFitter::Particles::kParton),
                                       particles->JetIndex(i),
                                       particles->IsBTagged(i),
                                       particles->BTaggingEfficiency(i),
                                       particles->BTaggingRejection(i),
                                       particles->TrueFlavor(i),
                                       particles->BTagWeight(i));
      
      
      
      // add index to map 
      fMapJets.push_back(i); 
      delete tlv_tmp;
    }
  }
  
  // run of non-b-tag jets
  for (int i = 0; i < npartons; ++i)
    {
      if (fNBJets.size() > 0) {
	// check non-b-tag
	if (particles->BTagWeight(i) > fNBJets.at(0).value) 
	  continue;
      }
      
      // check eta region 
      if (TMath::Abs( particles->DetEta(i, KLFitter::Particles::kParton) ) > fJetEta) 
        continue; 
      
      // check pT 
      if (particles->Parton(i)->Pt() < fJetPt)
        continue; 
      
      // add jet 
      TLorentzVector * tlv_tmp = new TLorentzVector(*particles->Parton(i));
      fParticlesSelected->AddParticle( tlv_tmp,
				       //				       new TLorentzVector(*particles->Parton(i)),
                                       particles->DetEta(i, KLFitter::Particles::kParton),
                                       KLFitter::Particles::kParton, 
                                       particles->NameParticle(i, KLFitter::Particles::kParton),
                                       particles->JetIndex(i),
                                       particles->IsBTagged(i),
                                       particles->BTaggingEfficiency(i),
                                       particles->BTaggingRejection(i),
                                       particles->TrueFlavor(i),
                                       particles->BTagWeight(i));
      
      // add index to map 
      fMapJets.push_back(i); 
      delete tlv_tmp;
    }
  
  
  
  
  // electron selection 
  int nelectrons = particles->NElectrons(); 

  if(nelectrons!=0 && particles->LeptonCharge(0, KLFitter::Particles::kElectron)!=-9)
    isDilepton=true;

  //std::cout << "isDilepton? " << isDilepton  << std::endl;
  
  for (int i = 0; i < nelectrons; ++i)
    {
      // check eta region 
      if (TMath::Abs(particles->DetEta(i, KLFitter::Particles::kElectron)) > fElectronEta) 
        continue; 
      
      // check pT 
      if (particles->Electron(i)->Pt() < fElectronPt)
        continue; 
      TLorentzVector * tlv_tmp = new TLorentzVector(*particles->Electron(i));

      if (isDilepton){
	// add electron 
	fParticlesSelected->AddParticle( tlv_tmp, 
					 //				       new TLorentzVector(*particles->Electron(i)),
					 particles->DetEta(i, KLFitter::Particles::kElectron),
					 particles->LeptonCharge(i, KLFitter::Particles::kElectron),
					 KLFitter::Particles::kElectron, 
					 particles->NameParticle(i, KLFitter::Particles::kElectron),
					 particles->ElectronIndex(i));
      }
      else {
	// add electron 
	fParticlesSelected->AddParticle( tlv_tmp, 
					 //				       new TLorentzVector(*particles->Electron(i)),
					 particles->DetEta(i, KLFitter::Particles::kElectron),
					 KLFitter::Particles::kElectron, 
					 particles->NameParticle(i, KLFitter::Particles::kElectron),
					 particles->ElectronIndex(i));
      }


      // add index to map 
      fMapElectrons.push_back(i); 
      delete tlv_tmp;
    }
  
  // muon selection 
  int nmuons = particles->NMuons(); 

  if(nmuons!=0 && particles->LeptonCharge(0, KLFitter::Particles::kMuon)!=-9)
    isDilepton=true;

  //std::cout << "isDilepton? " << isDilepton  << std::endl;
  
  for (int i = 0; i < nmuons; ++i)
    {
      // check eta region 
      if (TMath::Abs(particles->DetEta(i, KLFitter::Particles::kMuon)) > fMuonEta) 
        continue; 
      
      // check pT 
      if (particles->Muon(i)->Pt() < fMuonPt)
        continue; 
      TLorentzVector * tlv_tmp = new TLorentzVector(*particles->Muon(i));

      if (isDilepton) {
	// add muon 
	fParticlesSelected->AddParticle( tlv_tmp,  
					 //				       new TLorentzVector(*particles->Muon(i)),
					 particles->DetEta(i, KLFitter::Particles::kMuon),
					 particles->LeptonCharge(i, KLFitter::Particles::kMuon),
					 KLFitter::Particles::kMuon, 
					 particles->NameParticle(i, KLFitter::Particles::kMuon),
					 particles->MuonIndex(i));
      }
      else {
	// add muon 
	fParticlesSelected->AddParticle( tlv_tmp,  
					 //				       new TLorentzVector(*particles->Muon(i)),
					 particles->DetEta(i, KLFitter::Particles::kMuon),
					 KLFitter::Particles::kMuon, 
					 particles->NameParticle(i, KLFitter::Particles::kMuon),
					 particles->MuonIndex(i));
      }


      // add index to map 
      fMapMuons.push_back(i); 
      delete tlv_tmp;
    }
  
  // photon selection 
  int nphotons = particles->NPhotons(); 
  
  for (int i = 0; i < nphotons; ++i)
    {
      // check eta region 
      if (TMath::Abs(particles->DetEta(i, KLFitter::Particles::kPhoton)) > fPhotonEta) 
        continue; 

      // check pT 
      if (particles->Photon(i)->Pt() < fPhotonPt)
        continue; 
      TLorentzVector * tlv_tmp = new TLorentzVector(*particles->Photon(i));
      // add photon 
      fParticlesSelected->AddParticle( tlv_tmp,   
				       //				       new TLorentzVector(*particles->Photon(i)),
                                       particles->DetEta(i, KLFitter::Particles::kPhoton),
                                       KLFitter::Particles::kPhoton, 
                                       particles->NameParticle(i, KLFitter::Particles::kPhoton),
                                       particles->PhotonIndex(i));
      
      // add index to map 
      fMapPhotons.push_back(i); 
      delete tlv_tmp;
    }



  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::SelectionTool::SelectEvent(KLFitter::Particles * particles, double MET, double MWT)
{
  // get the object selection pt cuts from the event cuts
  fElectronPt = ObjectPtCut(fNElectronsPt);
  fMuonPt = ObjectPtCut(fNMuonsPt);
  fJetPt = ObjectPtCut(fNJetsPt);
  fPhotonPt = ObjectPtCut(fNPhotonsPt);
  
  
  // select objects
  if (!this->SelectObjects(particles))
    return 0; 
  
  // check if particles exist
  if (!fParticlesSelected)
    {
      std::cout << "KLFitter::SelectionTool::SelectEvent. Particles do not exist." << std::endl; 
      return 0; 
    }

  // increase coutner 
  fCounterEvents++; 
  
  // ------------------
  // electron selection 
  // ------------------

  if (int(fNElectronsPt.size()) > 0) {
    // counting variables
    int nelectrons = fParticlesSelected->NElectrons(); 
    unsigned int nelectroncuts = fNElectronsPt.size(); 
    
    std::vector<int> nelectronspt; 
    nelectronspt.assign(nelectroncuts, 0); 
    
    // loop over electrons
    for (int i = 0; i < nelectrons; ++i) {
      // get pt of electron
      double pt = fParticlesSelected->Electron(i)->Pt(); 
      
      // loop over all cuts and count
      for (unsigned int j = 0; j < nelectroncuts; ++j) {
        // increase counter if pt larger than cut value 
        if ( pt > fNElectronsPt[j].value)
          nelectronspt[j]++; 
      }
    }
    
    // check electron cuts
    for (unsigned int i = 0; i < nelectroncuts; ++i) {
      if (nelectronspt[i] != fNElectronsPt[i].n)
        return 0; 
    }
  }
  
  // increase counter 
  fCounterElectrons++; 
    
  // --------------
  // muon selection 
  // --------------

  if (int(fNMuonsPt.size()) > 0) {
    // counting variables
    int nmuons = fParticlesSelected->NMuons(); 
    int nmuoncuts = int(fNMuonsPt.size()); 
    
    std::vector<int> nmuonspt; 
    nmuonspt.assign(nmuoncuts, 0); 
    
    // loop over muons
    for (int i = 0; i < nmuons; ++i) {
      // get pt of muon
      double pt = fParticlesSelected->Muon(i)->Pt(); 
      
      // loop over all cuts and count
      for (int j = 0; j < nmuoncuts; ++j) {
        // increase counter if pt larger than cut value 
        if ( pt > fNMuonsPt[j].value)
          nmuonspt[j]++; 
      }
    }
    
    // check muon cuts
    for (int i = 0; i < nmuoncuts; ++i) {
      if (nmuonspt[i] != fNMuonsPt[i].n)
        return 0; 
    }
  }

  // increase counter 
  fCounterMuons++; 
  
  
  //------------
  // jet selection 
  //------------
  
  int nbjets = 0;
  
  if (int(fNJetsPt.size()) > 0)
    {
      // counting variables
      int njets = fParticlesSelected->NPartons(); 
      int njetcuts = int(fNJetsPt.size()); 	   
      std::vector<int> njetspt; 
      njetspt.assign(njetcuts, 0); 
      
      // loop over jets
      for (int i = 0; i < njets; ++i)
        {
          // get pt of jet
          double pt = fParticlesSelected->Parton(i)->Pt();                                       
	  double tag = fParticlesSelected->BTagWeight(i);
          // loop over all cuts and count
          for (int j = 0; j < njetcuts; ++j)
            {
              // increase counter if pt larger than cut value 
              if ( pt > fNJetsPt.at(j).value) {
                njetspt[j]++; 
		if (fNBJets.size()>0) {
		  if (tag > fNBJets.at(j).value)
		    nbjets++;
		}
	      }
            }
        }
      
      // check jet cuts
      for (int i = 0; i < njetcuts; ++i)
        {
          if (njetspt.at(i) < fNJetsPt.at(i).n)
            return 0;
	  
          if (fNJetsPt.at(i).dn >= 0 && njetspt.at(i) - fNJetsPt.at(i).n  > fNJetsPt.at(i).dn)
            return 0;
        }

      RemoveAdditionalParticles(int(fMaxNJetsForFit), KLFitter::Particles::kParton);
			
			// increase counter 
      fCounterJets++; 
      
    }
  else {
    // increase counter 
    fCounterJets++; 
  }
  
  // --------------
  // photon selection 
  // --------------

  if (int(fNPhotonsPt.size()) > 0)
    {
      // counting variables
      int nphotons = fParticlesSelected->NPhotons(); 
      int nphotoncuts = int(fNPhotonsPt.size()); 

      std::vector<int> nphotonspt; 
      nphotonspt.assign(nphotoncuts, 0); 

      // loop over photons
      for (int i = 0; i < nphotons; ++i)
        {
          // get pt of photon
          double pt = fParticlesSelected->Photon(i)->Pt(); 
                                        
          // loop over all cuts and count
          for (int j = 0; j < nphotoncuts; ++j)
            {
              // increase counter if pt larger than cut value 
              if ( pt > fNPhotonsPt.at(j).value)
                nphotonspt[j]++; 
            }
        }
      
      // check photon cuts
      for (int i = 0; i < nphotoncuts; ++i)
        {
          if (nphotonspt.at(i) != fNPhotonsPt.at(i).n)
            return 0; 
        }
    }

  // increase counter 
  fCounterPhotons++; 

  // MET selection
  if (MET < fMET)
    return 0; 
        
  // increase counter 
  fCounterMET++; 
  
  // MWT selection
  if (MWT < fMWT)
    return 0; 
  
  // increase counter 
  fCounterMWT++; 
  
  // triangular cut
  if ( (MWT+MET) < fMET_plus_MWT)
    return 0;
  
  // increase counter 
  fCounterTriangular++; 
  
  // check jet cuts
  if (fNBJets.size()>0) {
    
    if (nbjets < fNBJets.at(0).n)
      return 0;
    
    if (fNBJets.at(0).dn >= 0 && nbjets - fNBJets.at(0).n  > fNBJets.at(0).dn)
      return 0;
  }
  
  // increase counter 
  fCounterBJets++; 
  
  
  // increase counter 
  fCounterSelected++; 
    
  // event passed
  return 1;
}

// --------------------------------------------------------- 
int KLFitter::SelectionTool::RequireNBJets(double weight, int n, int dn){

  // add cut to set of cuts
  KLFitter::SelectionTool::Cut cut; 
  cut.value = weight; 
  cut.n = n; 
  cut.dn = dn; 
  fNBJets.push_back(cut);
  
  // no errors 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::SelectionTool::RequireNJetsPt(double pt, int n, int dn)
{
  if (pt < 0) 
    {
      std::cout << "KLFitter::SelectionTool::RequireNJetsPt. Pt < 0 does not make sense." << std::endl; 
      return 0; 
    }
  
  if (n < 0)
    {
      std::cout << "KLFitter::SelectionTool::RequireNJetsPt. n < 0 does not make sense." << std::endl; 
      return 0; 
    }

  // add cut to set of cuts
  KLFitter::SelectionTool::Cut cut; 
  cut.value = pt; 
  cut.n = n; 
  cut.dn = dn; 
  fNJetsPt.push_back(cut); 
  
  // no errors 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::SelectionTool::RequireNElectronsPt(double pt, int n)
{
  if (pt < 0) 
    {
      std::cout << "KLFitter::SelectionTool::RequireNElectronsPt. Pt < 0 does not make sense." << std::endl; 
      return 0; 
    }

  if (n < 0)
    {
      std::cout << "KLFitter::SelectionTool::RequireNElectronsPt. n < 0 does not make sense." << std::endl; 
      return 0; 
    }

  // add cut to set of cuts
  KLFitter::SelectionTool::Cut cut; 
  cut.value = pt; 
  cut.n = n; 
  cut.dn = 0; 
  fNElectronsPt.push_back(cut); 

  // no errors 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::SelectionTool::RequireNMuonsPt(double pt, int n)
{
  if (pt < 0) 
    {
      std::cout << "KLFitter::SelectionTool::RequireNMuonsPt. Pt < 0 does not make sense." << std::endl; 
      return 0; 
    }

  if (n < 0)
    {
      std::cout << "KLFitter::SelectionTool::RequireNMuonsPt. n < 0 does not make sense." << std::endl; 
      return 0; 
    }

  // add cut to set of cuts
  KLFitter::SelectionTool::Cut cut; 
  cut.value = pt; 
  cut.n = n; 
  cut.dn = 0; 
  fNMuonsPt.push_back(cut); 

  // no errors 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::SelectionTool::RequireNPhotonsPt(double pt, int n)
{
  if (pt < 0) 
    {
      std::cout << "KLFitter::SelectionTool::RequireNPhotonsPt. Pt < 0 does not make sense." << std::endl; 
      return 0; 
    }

  if (n < 0)
    {
      std::cout << "KLFitter::SelectionTool::RequireNPhotonsPt. n < 0 does not make sense." << std::endl; 
      return 0; 
    }

  // add cut to set of cuts
  KLFitter::SelectionTool::Cut cut; 
  cut.value = pt; 
  cut.n = n; 
  cut.dn = 0; 
  fNPhotonsPt.push_back(cut); 

  // no errors 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::SelectionTool::RequireMET(double met)
{
  if (met < 0)
    {
      std::cout << "KLFitter::SelectionTool::RequireMET. Missing ET < 0 does not make sense." << std::endl; 
      return 0; 
    }

  fMET = met; 

  // no errors
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::SelectionTool::RequireMWT(double mwt)
{
  if (mwt < 0)
    {
      std::cout << "KLFitter::SelectionTool::RequireMWT. MWT < 0 does not make sense." << std::endl; 
      return 0; 
    }

  fMWT = mwt; 

  // no errors
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::SelectionTool::RequireTriangular(double met_plus_mwt)
{
  if (met_plus_mwt < 0)
    {
      std::cout << "KLFitter::SelectionTool::RequireTriangular. MET+MWT < 0 does not make sense." << std::endl; 
      return 0; 
    }

  fMET_plus_MWT = met_plus_mwt; 

  // no errors
  return 1; 
}

// --------------------------------------------------------- 
void KLFitter::SelectionTool::ResetMaps()
{
  fMapJets.clear(); 
  fMapElectrons.clear(); 
  fMapMuons.clear();    
  fMapPhotons.clear();  
}

// --------------------------------------------------------- 
void KLFitter::SelectionTool::ResetCounter()
{
  fCounterEvents = 0; 
  fCounterJets = 0; 
  fCounterElectrons = 0; 
  fCounterMuons = 0; 
  fCounterPhotons = 0; 
  fCounterMET = 0;      
  fCounterMWT = 0;      
  fCounterTriangular = 0;      
  fCounterSelected = 0;      
}

// --------------------------------------------------------- 
double KLFitter::SelectionTool::ObjectPtCut(std::vector<Cut> const& cuts)
{
  double tmp_pt = -1.0;
  for (unsigned int iCut = 0; iCut < cuts.size(); iCut++)
    {
      if (tmp_pt > cuts.at(iCut).value || tmp_pt < 0)
        {
          tmp_pt = cuts.at(iCut).value;
        }
    }
  return tmp_pt;
}

// --------------------------------------------------------- 
void KLFitter::SelectionTool::RemoveAdditionalParticles(int n, KLFitter::Particles::ParticleType type)
{
  // remove particles - starting with lowest pT !
  while (fParticlesSelected->NParticles(type) > n)
    {
      double pTmin = 0.;
      int remove = -1;
      for (int j = 0; j < fParticlesSelected->NParticles(type); ++j)
        {
          double pt = fParticlesSelected->Particle(j, type)->Pt();
          if (j == 0 || pt < pTmin)
            {
              pTmin = pt;
              remove = j;
            }
        }
      fParticlesSelected->RemoveParticle(remove, type); 
      std::vector<int>::iterator it = fMapJets.begin();
      for (int i = 0; i < remove; i++)
        it++;
      //                        fMapJets.erase(it);
    }
}

// --------------------------------------------------------- 
