#include "SelectionTool.h" 
#include <iostream> 


// --------------------------------------------------------- 
KLFitter::SelectionTool::SelectionTool() :
	fParticlesSelected(new KLFitter::Particles()),
	//	fJetPt(20.0), 
	fJetPt(0.0), 
	fJetEta(2.5),
	//	fElectronPt(20.0),
	fElectronPt(0.0),
	fElectronEta(2.5), 
	fMuonPt(20.0),
	fMuonEta(2.5),
	//	fPhotonPt(10.0),
	fPhotonPt(0.0),
	fPhotonEta(2.5),
	fMET(0), 
	fCounterEvents(0),  
	fCounterJets(0),  
	fCounterElectrons(0),  
	fCounterMuons(0),  
	fCounterPhotons(0),  
	fCounterMET(0)
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
	fParticlesSelected = new KLFitter::Particles(); 

	// reset maps 
	this -> ResetMaps(); 

	// jet selection
	int npartons = particles -> NPartons(); 

	for (int i = 0; i < npartons; ++i)
		{
			// check eta region 
			if (TMath::Abs(particles -> Parton(i) -> Eta()) > fJetEta) 
				continue; 

			// check pT 
			if (particles -> Parton(i) -> Pt() < fJetPt)
				continue; 

			// add jet 
			fParticlesSelected -> AddParticle( particles -> Parton(i),
																				 KLFitter::Particles::kParton, 
																				 particles -> NameParticle(i, KLFitter::Particles::kParton) ); 

			// add flavor tag
			fParticlesSelected->SetFlavorTag( fParticlesSelected->NPartons()-1, particles->FlavorTag(i)); 

			// add b-tagging probability
			fParticlesSelected->SetBTaggingProbability( fParticlesSelected->NPartons()-1, particles->BTaggingProbability(i)); 
			// add index to map 
			fMapJets.push_back(i); 
		}

	// electron selection 
	int nelectrons = particles -> NElectrons(); 

	for (int i = 0; i < nelectrons; ++i)
		{
			// check eta region 
			if (TMath::Abs(particles -> Electron(i) -> Eta()) > fElectronEta) 
				continue; 

			// check pT 
			if (particles -> Electron(i) -> Pt() < fElectronPt)
				continue; 

			// add electron 
			fParticlesSelected -> AddParticle( particles -> Electron(i),
																				 KLFitter::Particles::kElectron, 
																				 particles -> NameParticle(i, KLFitter::Particles::kElectron) ); 
			// add index to map 
			fMapElectrons.push_back(i); 
		}

	// muon selection 
	int nmuons = particles -> NMuons(); 

	for (int i = 0; i < nmuons; ++i)
		{
			// check eta region 
			if (TMath::Abs(particles -> Muon(i) -> Eta()) > fMuonEta) 
				continue; 

			// check pT 
			if (particles -> Muon(i) -> Pt() < fMuonPt)
				continue; 

			// add muon 
			fParticlesSelected -> AddParticle( particles -> Muon(i),
																				 KLFitter::Particles::kMuon, 
																				 particles -> NameParticle(i, KLFitter::Particles::kMuon) ); 

			// add index to map 
			fMapMuons.push_back(i); 
		}

	// photon selection 
	int nphotons = particles -> NPhotons(); 

	for (int i = 0; i < nphotons; ++i)
		{
			// check eta region 
			if (TMath::Abs(particles -> Photon(i) -> Eta()) > fPhotonEta) 
				continue; 

			// check pT 
			if (particles -> Photon(i) -> Pt() < fPhotonPt)
				continue; 

			// add photon 
			fParticlesSelected -> AddParticle( particles -> Photon(i),
																				 KLFitter::Particles::kPhoton, 
																				 particles -> NameParticle(i, KLFitter::Particles::kPhoton) ); 

			// add index to map 
			fMapPhotons.push_back(i); 
		}

	// no error 
	return 1; 
}

// --------------------------------------------------------- 
int KLFitter::SelectionTool::SelectEvent(KLFitter::Particles * particles, double MET)
{
	// select objects
	if (!this -> SelectObjects(particles))
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

	if (int(fNElectronsPt.size()) > 0)
		{
			// counting variables
			int nelectrons = fParticlesSelected -> NElectrons(); 
			int nelectroncuts = int(fNElectronsPt.size()); 
			
			std::vector<int> nelectronspt; 
			nelectronspt.assign(nelectroncuts, 0); 

			// loop over electrons
			for (int i = 0; i < nelectrons; ++i)
				{
					// get pt of electron
					double pt = fParticlesSelected -> Electron(i) -> Pt(); 
					
					// loop over all cuts and count
					for (int j = 0; j < nelectroncuts; ++j)
						{
							// increase counter if pt larger than cut value 
							if ( pt > fNElectronsPt.at(j).value)
								nelectronspt[j]++; 
						}
				}

			// check electron cuts
			for (int i = 0; i < nelectroncuts; ++i)
				{
					if (fNElectronsPt.at(i).dn > 0 && nelectronspt.at(i) < fNElectronsPt.at(i).n)
						return 0; 
					else if (fNElectronsPt.at(i).dn <= 0 && nelectronspt.at(i) != fNElectronsPt.at(i).n)
						return 0; 
				}
		}

	// increase counter 
	fCounterElectrons++; 

	// --------------
	// muon selection 
	// --------------

	if (int(fNMuonsPt.size()) > 0)
		{
			// counting variables
			int nmuons = fParticlesSelected -> NMuons(); 
			int nmuoncuts = int(fNMuonsPt.size()); 

			std::vector<int> nmuonspt; 
			nmuonspt.assign(nmuoncuts, 0); 

			// loop over muons
			for (int i = 0; i < nmuons; ++i)
				{
					// get pt of muon
					double pt = fParticlesSelected -> Muon(i) -> Pt(); 
					
					// loop over all cuts and count
					for (int j = 0; j < nmuoncuts; ++j)
						{
							// increase counter if pt larger than cut value 
							if ( pt > fNMuonsPt.at(j).value)
								nmuonspt[j]++; 
						}
				}
			
			// check muon cuts
			for (int i = 0; i < nmuoncuts; ++i)
				{
					if (fNMuonsPt.at(i).dn > 0 && nmuonspt.at(i) < fNMuonsPt.at(i).n)
						return 0; 
					else if (fNMuonsPt.at(i).dn <= 0 && nmuonspt.at(i) != fNMuonsPt.at(i).n)
						return 0; 
				}
		}

	// increase counter 
	fCounterMuons++; 
			
	//------------
	// jet selection 
	//------------

	if (int(fNJetsPt.size()) > 0)
		{
			// counting variables
			int njets = fParticlesSelected -> NPartons(); 
			int njetcuts = int(fNJetsPt.size()); 
			
			std::vector<int> njetspt; 
			njetspt.assign(njetcuts, 0); 
			
			// loop over jets
			for (int i = 0; i < njets; ++i)
				{
					// get pt of jet
					double pt = fParticlesSelected -> Parton(i) -> Pt(); 
					
					// loop over all cuts and count
					for (int j = 0; j < njetcuts; ++j)
						{
							// increase counter if pt larger than cut value 
							if ( pt > fNJetsPt.at(j).value)
								njetspt[j]++; 
						}
				}
			
			// check jet cuts
			for (int i = 0; i < njetcuts; ++i)
				{
					if (fNJetsPt.at(i).dn > 0)
						{
							if (njetspt.at(i) < fNJetsPt.at(i).n)
								return 0; 
						}
					else
						{
							if (njetspt.at(i) < fNJetsPt.at(i).n)
								return 0; 

							// remove particles 
							// - starting with lowest pT !
							else if (njetspt.at(i) > fNJetsPt.at(i).n)
								{
									while (fParticlesSelected -> NPartons() > fNJetsPt.at(i).n)
										{
											double pTmin = 0.;
											int jetToRemove = -1;
											for (int j = 0; j < fParticlesSelected -> NPartons(); ++j)
												{
													double pt = fParticlesSelected -> Parton(j) -> Pt();
													if (j == 0 || pt < pTmin)
														{
															pTmin = pt;
															jetToRemove = j;
														}
												}
											fParticlesSelected -> RemoveParticle(jetToRemove, KLFitter::Particles::kParton); 
										}
								}
						}
				}
		}

	// increase counter 
	fCounterJets++; 
	
	// --------------
	// photon selection 
	// --------------

	if (int(fNPhotonsPt.size()) > 0)
		{
			// counting variables
			int nphotons = fParticlesSelected -> NPhotons(); 
			int nphotoncuts = int(fNPhotonsPt.size()); 

			std::vector<int> nphotonspt; 
			nphotonspt.assign(nphotoncuts, 0); 

			// loop over photons
			for (int i = 0; i < nphotons; ++i)
				{
					// get pt of photon
					double pt = fParticlesSelected -> Photon(i) -> Pt(); 
					
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
					if (fNPhotonsPt.at(i).dn > 0 && nphotonspt.at(i) < fNPhotonsPt.at(i).n)
						return 0; 
					else if (fNPhotonsPt.at(i).dn <= 0 && nphotonspt.at(i) != fNPhotonsPt.at(i).n)
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
	
	// event passed
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
int KLFitter::SelectionTool::RequireNElectronsPt(double pt, int n, int dn)
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
	cut.dn = dn; 
	fNElectronsPt.push_back(cut); 

	// no errors 
	return 1; 
}

// --------------------------------------------------------- 
int KLFitter::SelectionTool::RequireNMuonsPt(double pt, int n, int dn)
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
	cut.dn = dn; 
	fNMuonsPt.push_back(cut); 

	// no errors 
	return 1; 
}

// --------------------------------------------------------- 
int KLFitter::SelectionTool::RequireNPhotonsPt(double pt, int n, int dn)
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
	cut.dn = dn; 
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
}

// --------------------------------------------------------- 

