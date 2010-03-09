#include "LikelihoodTopLeptonJets.h" 
#include <iostream> 
#include <BAT/BCMath.h> 

// --------------------------------------------------------- 
KLFitter::LikelihoodTopLeptonJets::LikelihoodTopLeptonJets() : KLFitter::LikelihoodBase::LikelihoodBase()
{
	// define model particles 
	this->DefineModelParticles(); 

	// define parameters 
	this->DefineParameters(); 

	// initialize missing Et
	ETmiss_x = 0.;
	ETmiss_y = 0.;

	// initialize flags
	fFlagTopMassFixed = false; 
	fFlagUseJetMass = false; 

	// initialize lepton type 
	fTypeLepton = 1; 
}

// --------------------------------------------------------- 
KLFitter::LikelihoodTopLeptonJets::~LikelihoodTopLeptonJets()
{
	while (!fLorentzVectorContainer->empty())
		{
			TLorentzVector * lv = fLorentzVectorContainer->front(); 
			fLorentzVectorContainer->erase(fLorentzVectorContainer ->begin() ); 
			delete lv; 
		}
	delete fLorentzVectorContainer; 

	if (fParticlesModel)
		delete fParticlesModel; 
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTopLeptonJets::SetET_miss_XY(double etx, double ety)
{
	// set missing ET x and y component
	ETmiss_x = etx;
	ETmiss_y = ety;

	// no error
	return 1;
}

// --------------------------------------------------------- 
void KLFitter::LikelihoodTopLeptonJets::SetLeptonType(int leptontype)
{
	if (leptontype < 1 || leptontype > 2)
		{
			std::cout << " Warnung: lepton type not defined. Set electron as lepton type." << std::endl; 
			fTypeLepton = 1;
		}
	else
		fTypeLepton = leptontype; 

	// define model particles 
	this->DefineModelParticles(); 
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTopLeptonJets::DefineModelParticles()
{
	// check if model particles and lorentz vector container exist and delete
	this->DeleteModelParticles();

	// create the particles of the model 
	fParticlesModel = new KLFitter::Particles(); 
	fLorentzVectorContainer = new std::vector <TLorentzVector *>(0); 
	TLorentzVector * lv = 0; 

	// add model particles 
	fParticlesModel->AddParticle(lv = new TLorentzVector(),    // 4-vector 
															 KLFitter::Particles::kParton, // type 
															 "hadronic b quark",           // name 
															 1.0,                          // b-jet (truth) 
															 0.0,                          // not tagged (meas)
															 0);                           // index of corresponding particle 
	fLorentzVectorContainer->push_back(lv); 

	fParticlesModel->AddParticle(lv = new TLorentzVector(),
															 KLFitter::Particles::kParton, 
															 "leptonic b quark",
															 1.0, 
															 0.0,                          // not tagged (meas)
															 1);                           // index of corresponding particle 
	fLorentzVectorContainer->push_back(lv); 	

	fParticlesModel->AddParticle(lv = new TLorentzVector(),
															 KLFitter::Particles::kParton,
															 "light quark 1",
															 0.0,                          // light jet (truth)
															 0.0,                          // not tagged (meas)
															 2);                           // index of corresponding particle 
	fLorentzVectorContainer->push_back(lv); 	

	fParticlesModel->AddParticle(lv = new TLorentzVector(),
																 KLFitter::Particles::kParton,
																 "light quark 2",
															 0.0, 
															 0.0,                          // not tagged (meas)
															 3);                           // index of corresponding particle
	fLorentzVectorContainer->push_back(lv); 
	
	if (fTypeLepton == 1)
		{
			fParticlesModel->AddParticle(lv = new TLorentzVector(),
																		 KLFitter::Particles::kElectron,
																		 "electron"); 
			fLorentzVectorContainer->push_back(lv); 
		}
	else if (fTypeLepton == 2)
		{
			fParticlesModel->AddParticle(lv = new TLorentzVector(),
																		 KLFitter::Particles::kMuon,
																		 "muon"); 
			fLorentzVectorContainer->push_back(lv); 
		}

	fParticlesModel->AddParticle(lv = new TLorentzVector(),
																 KLFitter::Particles::kNeutrino, 
																 "neutrino"); 
	fLorentzVectorContainer->push_back(lv); 

	fParticlesModel->AddParticle(lv = new TLorentzVector(),
																 KLFitter::Particles::kBoson, 
																 "hadronic W"); 
	fLorentzVectorContainer->push_back(lv); 

	fParticlesModel->AddParticle(lv = new TLorentzVector(),
																 KLFitter::Particles::kBoson,
																 "leptonic W"); 
	fLorentzVectorContainer->push_back(lv); 

	fParticlesModel->AddParticle(lv = new TLorentzVector(),
																 KLFitter::Particles::kParton,
																 "hadronic top",
																 2.0); 
	fLorentzVectorContainer->push_back(lv); 

	fParticlesModel->AddParticle(lv = new TLorentzVector(),
																 KLFitter::Particles::kParton,
																 "leptonic top",
																 2.0); 
	fLorentzVectorContainer->push_back(lv); 

	// no error 
	return 1;
}

// --------------------------------------------------------- 
void KLFitter::LikelihoodTopLeptonJets::DefineParameters()
{
	// add parameters of model
	this->AddParameter("energy hadronic b",       fPhysicsConstants->MassBottom(), 1000.0); // par 0
	this->AddParameter("energy leptonic b",       fPhysicsConstants->MassBottom(), 1000.0); // par 1
	this->AddParameter("energy light quark 1",    0.0, 1000.0);                               // par 2
	this->AddParameter("energy light quark 2",    0.0, 1000.0);                               // par 3
	this->AddParameter("energy lepton",           0.0, 1000.0);                               // par 4
	this->AddParameter("p_z neutrino",        -1000.0, 1000.0);                               // par 5
	this->AddParameter("missPx",              -1000.0, 1000.0);                               // par 6
	this->AddParameter("missPy",              -1000.0, 1000.0);                               // par 7
	this->AddParameter("eta hadronic b",       -2.5, 2.5);                                    // par 8
	this->AddParameter("eta leptonic b",       -2.5, 2.5);                                    // par 9
	this->AddParameter("eta light quark 1",    -2.5, 2.5);                                    // par 10
	this->AddParameter("eta light quark 2",    -2.5, 2.5);                                    // par 11
	this->AddParameter("phi hadronic b",       -TMath::Pi(), TMath::Pi());                    // par 12
	this->AddParameter("phi leptonic b",       -TMath::Pi(), TMath::Pi());                    // par 13
	this->AddParameter("phi light quark 1",    -TMath::Pi(), TMath::Pi());                    // par 14
	this->AddParameter("phi light quark 2",    -TMath::Pi(), TMath::Pi());                    // par 15
	this->AddParameter("top mass",             100.0, 1000.0);                                // par 16
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTopLeptonJets::CalculateLorentzVectors(std::vector <double> parameters)
{
	// variables
  double E, m, px, py, pz, pabs, theta, scale; 
	TLorentzVector * vect; 

	// hadronic b quark 
	E = parameters.at(0); 
	if (fFlagUseJetMass)
		m = TMath::Max(0., (*fParticlesPermuted)->Parton(0)->M()); 
	else
		m = fPhysicsConstants->MassBottom(); 
	pabs = sqrt(E*E - m*m);
	theta = atan(exp(-parameters.at(8)))*2;
	px = pabs * cos(parameters.at(12)) * sin(theta);
	py = pabs * sin(parameters.at(12)) * sin(theta);
	pz = pabs * cos(theta);
	fParticlesModel->Parton(0)->SetPxPyPzE(px, py, pz, E);

	// leptonic b quark 
	E = parameters.at(1); 
	if (fFlagUseJetMass)
		m = TMath::Max(0., (*fParticlesPermuted)->Parton(1)->M());
	else
		m = fPhysicsConstants->MassBottom(); 
	pabs = sqrt(E*E - m*m);
	theta = atan(exp(-parameters.at(9)))*2;
	px = pabs * cos(parameters.at(13)) * sin(theta);
	py = pabs * sin(parameters.at(13)) * sin(theta);
	pz = pabs * cos(theta);
	fParticlesModel->Parton(1)->SetPxPyPzE(px, py, pz, E);

	// light quark 1 
	E = parameters.at(2); 
	if (fFlagUseJetMass)
		m = TMath::Max(0., (*fParticlesPermuted)->Parton(2)->M());
	else
		m = 0; 
	pabs = sqrt(E*E - m*m);
	theta = atan(exp(-parameters.at(10)))*2;
	px = pabs * cos(parameters.at(14)) * sin(theta);
	py = pabs * sin(parameters.at(14)) * sin(theta);
	pz = pabs * cos(theta);
	fParticlesModel->Parton(2)->SetPxPyPzE(px, py, pz, E);

	// light quark 2 
	E = parameters.at(3); 
	if (fFlagUseJetMass)
		m = TMath::Max(0., (*fParticlesPermuted)->Parton(3)->M());
	else
		m = 0; 
	pabs = sqrt(E*E - m*m);
	theta = atan(exp(-parameters.at(11)))*2;
	px = pabs * cos(parameters.at(15)) * sin(theta);
	py = pabs * sin(parameters.at(15)) * sin(theta);
	pz = pabs * cos(theta);
	fParticlesModel->Parton(3)->SetPxPyPzE(px, py, pz, E);

	// lepton
	if (fTypeLepton == 1)
		{
			vect = (*fParticlesPermuted)->Electron(0); 
			E = parameters.at(4); 
			px = vect->Px(); 
			py = vect->Py(); 
			pz = vect->Pz(); 
			pabs = sqrt(px*px + py*py + pz*pz); 
			scale = E / vect->E(); 
			fParticlesModel->Electron(0)->SetPxPyPzE(scale * px, scale * py, scale * pz, E); 
		}
	else if (fTypeLepton == 2)
		{
			vect = (*fParticlesPermuted)->Muon(0); 
			E = parameters.at(4); 
			px = vect->Px(); 
			py = vect->Py(); 
			pz = vect->Pz(); 
			pabs = sqrt(px*px + py*py + pz*pz); 
			scale = E / vect->E(); 
			fParticlesModel->Muon(0)->SetPxPyPzE(scale * px, scale * py, scale * pz, E); 
		}

	// neutrino 
	px = parameters.at(6); // yes: (px, py, pz) = (par6, par7, par5) - this order (!)
	py = parameters.at(7);
	pz = parameters.at(5);
	E = sqrt(px*px + py*py + pz*pz);
	fParticlesModel->Neutrino(0)->SetPxPyPzE(px, py, pz, E); 

	// composite particles 

	// hadronic W 
	*(fParticlesModel->Boson(0)) = *(fParticlesModel->Parton(2)) + *(fParticlesModel->Parton(3)); 
	
	// leptonic W 
	if (fTypeLepton == 1)
		*(fParticlesModel->Boson(1)) = *(fParticlesModel->Electron(0)) + *(fParticlesModel->Neutrino(0)); 
	else if (fTypeLepton == 2)
		*(fParticlesModel->Boson(1)) = *(fParticlesModel->Muon(0)) + *(fParticlesModel->Neutrino(0)); 

	// hadronic top 
	*(fParticlesModel->Parton(4)) = *(fParticlesModel->Boson(0)) + *(fParticlesModel->Parton(0)); 

	// leptonic top 
	*(fParticlesModel->Parton(5)) = *(fParticlesModel->Boson(1)) + *(fParticlesModel->Parton(1)); 

	// no error 
	return 1; 
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTopLeptonJets::Initialize()
{
	// error code 
	int err = 1; 

	// adjust parameter ranges 
	err *= this->AdjustParameterRanges(); 

	// set initial values
	// (only for Markov chains - initial parameters for other minimisation methods are set in Fitter.cxx)
	this->SetInitialParameters( this->GetInitialParameters() );	

	// return error code 
	return err; 
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTopLeptonJets::RemoveInvariantParticlePermutations()
{
	// error code 
	int err = 1; 

	// remove the permutation from the second and the third jet
	KLFitter::Particles::ParticleType ptype = KLFitter::Particles::kParton;
	std::vector<int> indexVector_Jets;
	indexVector_Jets.push_back(2);
	indexVector_Jets.push_back(3);
	err *= (*fPermutations)->InvariantParticlePermutations(ptype, indexVector_Jets); 
	
	//remove invariant jet permutations of notevent jets
	KLFitter::Particles * particles = (*fPermutations)->Particles();
	indexVector_Jets.clear();
	for (int iPartons = 4; iPartons < particles->NPartons(); iPartons++)
		indexVector_Jets.push_back(iPartons);
	err *= (*fPermutations)->InvariantParticlePermutations(ptype, indexVector_Jets);

	// remove the permutation from the muons
	ptype = KLFitter::Particles::kMuon;
	std::vector<int> indexVector_Muons;
	for (int iMuon = 0; iMuon < particles->NMuons(); iMuon++)
		indexVector_Muons.push_back(iMuon);
	err *= (*fPermutations)->InvariantParticlePermutations(ptype, indexVector_Muons); 

	// return error code 
	return err; 
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTopLeptonJets::AdjustParameterRanges()
{
	// adjust limits 
	double nsigmas_jet = 7.0; 
	double nsigmas_lepton = 2.0; 

	// parameter 0: energy of hadronic b quark 
	double E = (*fParticlesPermuted)->Parton(0)->E(); 
	double m = fPhysicsConstants->MassBottom(); 
	if (fFlagUseJetMass)
		m = TMath::Max(0.0, (*fParticlesPermuted)->Parton(0)->M()); 
	double Emin = TMath::Max(m, E - nsigmas_jet * sqrt(E)); 
	double Emax  = E + nsigmas_jet * sqrt(E); 
	this->SetParameterRange(0, Emin, Emax); 

	// parameter 1: energy of leptonic b quark 
	E = (*fParticlesPermuted)->Parton(1)->E(); 
	m = fPhysicsConstants->MassBottom(); 
	if (fFlagUseJetMass)
		m = TMath::Max(0.0, (*fParticlesPermuted)->Parton(1)->M()); 
	Emin = TMath::Max(m, E - nsigmas_jet * sqrt(E)); 
	Emax  = E + nsigmas_jet * sqrt(E); 
	this->SetParameterRange(1, Emin, Emax); 

	// parameter 2: energy of light quark 1 
	E = (*fParticlesPermuted)->Parton(2)->E(); 
	m = 0.001;
	if (fFlagUseJetMass)
		m = TMath::Max(0.0, (*fParticlesPermuted)->Parton(2)->M()); 
	Emin = TMath::Max(m, E - nsigmas_jet * sqrt(E)); 
	Emax  = E + nsigmas_jet * sqrt(E); 
	this->SetParameterRange(2, Emin, Emax); 

	// parameter 3: energy of light quark2 
	E = (*fParticlesPermuted)->Parton(3)->E(); 
	m = 0.001;
	if (fFlagUseJetMass)
		m = TMath::Max(0.0, (*fParticlesPermuted)->Parton(3)->M()); 
	Emin = TMath::Max(m, E - nsigmas_jet * sqrt(E)); 
	Emax  = E + nsigmas_jet * sqrt(E); 
	this->SetParameterRange(3, Emin, Emax); 

	// noteKK: range for muons has to be adjusted (no sqrt(E) behaviour) 

	// parameter 2: energy of lepton
	if (fTypeLepton == 1)
		E = (*fParticlesPermuted)->Electron(0)->E(); 
	else if (fTypeLepton == 2)
		E = (*fParticlesPermuted)->Muon(0)->E(); 
	Emin = TMath::Max(0.001, E - nsigmas_lepton * sqrt(E)); 
	Emax  = E + nsigmas_lepton * sqrt(E); 
	this->SetParameterRange(4, Emin, Emax); 

	// note: this is hard-coded in the momement 

	// missing ET 
	this->SetParameterRange(6, ETmiss_x-100.0, ETmiss_x+100);
	this->SetParameterRange(7, ETmiss_y-100.0, ETmiss_y+100);

	// eta
	double eta = (*fParticlesPermuted)->Parton(0)->Eta(); 
	double etamin = TMath::Max(-2.5, eta - 0.2); 
	double etamax = TMath::Min(2.5, eta + 0.2); 
	this->SetParameterRange(8, etamin, etamax); 

	eta = (*fParticlesPermuted)->Parton(1)->Eta(); 
	etamin = TMath::Max(-2.5, eta - 0.2); 
	etamax = TMath::Min(2.5, eta + 0.2); 
	this->SetParameterRange(9, etamin, etamax); 

	eta = (*fParticlesPermuted)->Parton(2)->Eta(); 
	etamin = TMath::Max(-2.5, eta - 0.2); 
	etamax = TMath::Min(2.5, eta + 0.2); 
	this->SetParameterRange(10, etamin, etamax); 

	eta = (*fParticlesPermuted)->Parton(3)->Eta(); 
	etamin = TMath::Max(-2.5, eta - 0.2); 
	etamax = TMath::Min(2.5, eta + 0.2); 
	this->SetParameterRange(11, etamin, etamax); 

	// phi 
	double phi =  eta = (*fParticlesPermuted)->Parton(0)->Phi(); 
	double phimin = phi - 0.1;
	double phimax = phi + 0.1;
	this->SetParameterRange(12, phimin, phimax); 

	phi =  eta = (*fParticlesPermuted)->Parton(1)->Phi(); 
	phimin = phi - 0.1;
	phimax = phi + 0.1;
	this->SetParameterRange(13, phimin, phimax); 

	phi =  eta = (*fParticlesPermuted)->Parton(2)->Phi(); 
	phimin = phi - 0.1;
	phimax = phi + 0.1;
	this->SetParameterRange(14, phimin, phimax); 

	phi =  eta = (*fParticlesPermuted)->Parton(3)->Phi(); 
	phimin = phi - 0.1;
	phimax = phi + 0.1;
	this->SetParameterRange(15, phimin, phimax); 

	// top mass 
	if (fFlagTopMassFixed)
		this->SetParameterRange(16, this->fPhysicsConstants->MassTop(), this->fPhysicsConstants->MassTop()); 

	// no error 
	return 1; 
}

// --------------------------------------------------------- 
double KLFitter::LikelihoodTopLeptonJets::LogLikelihood(std::vector <double> parameters)
{
//	// check if W mass is within range
//	TLorentzVector Whad; 
//	double px, py, pz, E; 
//	px = (*fParticlesPermuted)->Parton(2)->Px() + (*fParticlesPermuted)->Parton(3)->Px(); 
//	py = (*fParticlesPermuted)->Parton(2)->Py() + (*fParticlesPermuted)->Parton(3)->Py(); 
//	pz = (*fParticlesPermuted)->Parton(2)->Pz() + (*fParticlesPermuted)->Parton(3)->Pz(); 
//	E  = (*fParticlesPermuted)->Parton(2)->E()  + (*fParticlesPermuted)->Parton(3)->E(); 
//	Whad.SetPxPyPzE(px, py, pz, E); 
//
//	if ( fabs(Whad.M() - fPhysicsConstants->MassW()) > 40.0 )
//		return -500; 

	// calculate 4-vectors 
	this->CalculateLorentzVectors(parameters); 

	// define log of likelihood 
	double logprob = 0.; 

	// jet energy resolution terms 
	logprob += log( (*fDetector)->ResEnergyBJet( parameters.at(8) )->p( parameters.at(0), (*fParticlesPermuted)->Parton(0)->E()) ); 
	logprob += log( (*fDetector)->ResEnergyBJet( parameters.at(9) )->p( parameters.at(1), (*fParticlesPermuted)->Parton(1)->E()) ); 
	logprob += log( (*fDetector)->ResEnergyLightJet( parameters.at(10) )->p( parameters.at(2), (*fParticlesPermuted)->Parton(2)->E()) ); 
	logprob += log( (*fDetector)->ResEnergyLightJet( parameters.at(11) )->p( parameters.at(3), (*fParticlesPermuted)->Parton(3)->E()) ); 

	// lepton energy resolution terms 
	if (fTypeLepton == 1)
		logprob += log( (*fDetector)->ResEnergyElectron( (*fParticlesPermuted)->Electron(0)->Eta() )->p( parameters.at(4), (*fParticlesPermuted)->Electron(0)->E()) ); 
	else if (fTypeLepton == 2) {
     //    logprob += log( (*fDetector)->ResEnergyMuon( (*fParticlesPermuted)->Muon(0)->Eta() )->p( parameters.at(4),		(*fParticlesPermuted)->Muon(0)->E()) );
     double pt_fit = parameters.at(4) * sin((*fParticlesPermuted)->Muon(0)->Theta());
     double pt_reco = (*fParticlesPermuted)->Muon(0)->Pt();
     double pt2 = pt_fit * pt_fit;
     // (pt)^2 == 0 cannot happen if event selection does not allow for eta = +/- infinity
     logprob += log( (*fDetector)->ResEnergyMuon( (*fParticlesPermuted)->Muon(0)->Eta() )->p( pt_fit/pt2, pt_reco/pt2 ));
	}

	//	ETmiss_x and ETmiss_y
	logprob += log( (*fDetector)->ResMissingET()->p( parameters.at(6), ETmiss_x) ); 
	logprob += log( (*fDetector)->ResMissingET()->p( parameters.at(7), ETmiss_y) ); 

	// eta resolution 
	logprob += log( (*fDetector)->ResEtaBJet(parameters.at(8))->p(parameters.at(8), (*fParticlesPermuted)->Parton(0)->Eta()) ); 
	logprob += log( (*fDetector)->ResEtaBJet(parameters.at(9))->p(parameters.at(9), (*fParticlesPermuted)->Parton(1)->Eta()) ); 
	logprob += log( (*fDetector)->ResEtaLightJet(parameters.at(10))->p(parameters.at(10), (*fParticlesPermuted)->Parton(2)->Eta()) );
	logprob += log( (*fDetector)->ResEtaLightJet(parameters.at(11))->p(parameters.at(11), (*fParticlesPermuted)->Parton(3)->Eta()) ); 	

 	// check phi variable
 	if (parameters.at(12) < -TMath::Pi())
 		parameters[12] += TMath::TwoPi(); 
 	if (parameters.at(12) > TMath::Pi())
 		parameters[12] -= TMath::TwoPi(); 

 	if (parameters.at(13) < -TMath::Pi())
 		parameters[13] += TMath::TwoPi(); 
 	if (parameters.at(13) > TMath::Pi())
 		parameters[13] -= TMath::TwoPi(); 

 	if (parameters.at(14) < -TMath::Pi())
 		parameters[14] += TMath::TwoPi(); 
 	if (parameters.at(14) > TMath::Pi())
 		parameters[14] -= TMath::TwoPi(); 

 	if (parameters.at(15) < -TMath::Pi())
 		parameters[15] += TMath::TwoPi(); 
 	if (parameters.at(15) > TMath::Pi())
 		parameters[15] -= TMath::TwoPi(); 
	
 	// phi resolution (this implementation assumes a symmetric phi resolution function
 	double diff0 = fabs(parameters.at(12)-(*fParticlesPermuted)->Parton(0)->Phi());
 	double diff1 = fabs(parameters.at(13)-(*fParticlesPermuted)->Parton(1)->Phi());
 	double diff2 = fabs(parameters.at(14)-(*fParticlesPermuted)->Parton(2)->Phi());
 	double diff3 = fabs(parameters.at(15)-(*fParticlesPermuted)->Parton(3)->Phi());

	if (diff0>TMath::Pi())
		diff0 -= TMath::TwoPi();
	if (diff1>TMath::Pi())
		diff1 -= TMath::TwoPi();
	if (diff2>TMath::Pi())
		diff2 -= TMath::TwoPi();
	if (diff3>TMath::Pi())
		diff3 -= TMath::TwoPi();

	logprob += log( (*fDetector)->ResPhiBJet(parameters.at(8))->p(0.0, diff0) ); 
	logprob += log( (*fDetector)->ResPhiBJet(parameters.at(9))->p(0.0, diff1) ); 
	logprob += log( (*fDetector)->ResPhiLightJet(parameters.at(10))->p(0.0, diff2) ); 
	logprob += log( (*fDetector)->ResPhiLightJet(parameters.at(11))->p(0.0, diff3) ); 

	// Breit-Wigner of hadronically decaying W-boson
	logprob += BCMath::LogBreitWignerRel( fParticlesModel->Boson(0)->M(), fPhysicsConstants->MassW(), fPhysicsConstants->GammaW()); 

	// Breit-Wigner of leptonically decaying W-boson
	logprob += BCMath::LogBreitWignerRel( fParticlesModel->Boson(1)->M(), fPhysicsConstants->MassW(), fPhysicsConstants->GammaW()); 

	// note: top mass width should be made DEPENDENT on the top mass at a certain point
	//	fPhysicsConstants->SetMassTop(parameters.at(16));

	// Breit-Wigner of hadronically decaying top quark
	logprob += BCMath::LogBreitWignerRel( fParticlesModel->Parton(4)->M(), parameters.at(16), fPhysicsConstants->GammaTop()); 
	
	// Breit-Wigner of leptonically decaying top quark
	logprob += BCMath::LogBreitWignerRel( fParticlesModel->Parton(5)->M(), parameters.at(16), fPhysicsConstants->GammaTop()); 

	// return log of likelihood 
	return logprob; 
}

// --------------------------------------------------------- 
std::vector<double> KLFitter::LikelihoodTopLeptonJets::GetInitialParameters()
{
	std::vector<double> values(this->GetInitialParametersWoNeutrinoPz());

	// check second neutrino solution
	std::vector<double> neutrino_pz_solutions = GetNeutrinoPzSolutions();
	if (int(neutrino_pz_solutions.size()) == 1)
		values[5] = neutrino_pz_solutions.at(0); 
	else if(int(neutrino_pz_solutions.size()) == 2)
		{
			double sol1, sol2; 
			values[5] = neutrino_pz_solutions.at(0); 
			sol1 = this->LogLikelihood(values); 
			values[5] = neutrino_pz_solutions.at(1); 
			sol2 = this->LogLikelihood(values); 

			if (sol1 > sol2)
				values[5] = neutrino_pz_solutions.at(0); 
		}

	return values;
}

// --------------------------------------------------------- 
std::vector<double> KLFitter::LikelihoodTopLeptonJets::GetInitialParametersWoNeutrinoPz()
{
	std::vector<double> values;

	// energies of the quarks
	values.push_back( (*fParticlesPermuted)->Parton(0)->E() );
	values.push_back( (*fParticlesPermuted)->Parton(1)->E() );
	values.push_back( (*fParticlesPermuted)->Parton(2)->E() );
	values.push_back( (*fParticlesPermuted)->Parton(3)->E() );

	// energy of the lepton
	if (fTypeLepton == 1)
		values.push_back( (*fParticlesPermuted)->Electron(0)->E() );
	else if (fTypeLepton == 2)
		values.push_back( (*fParticlesPermuted)->Muon(0)->E() );

	// pz of the neutrino
	values.push_back(0.0); 

	// missing px and py
	values.push_back(ETmiss_x);
	values.push_back(ETmiss_y);

	// eta and phi values of the jets
	values.push_back( (*fParticlesPermuted)->Parton(0)->Eta() );
	values.push_back( (*fParticlesPermuted)->Parton(1)->Eta() );
	values.push_back( (*fParticlesPermuted)->Parton(2)->Eta() );
	values.push_back( (*fParticlesPermuted)->Parton(3)->Eta() );
	values.push_back( (*fParticlesPermuted)->Parton(0)->Phi() );
	values.push_back( (*fParticlesPermuted)->Parton(1)->Phi() );
	values.push_back( (*fParticlesPermuted)->Parton(2)->Phi() );
	values.push_back( (*fParticlesPermuted)->Parton(3)->Phi() );

	// top mass
	double mtop = ( *(*fParticlesPermuted)->Parton(0) + *(*fParticlesPermuted)->Parton(2) + *(*fParticlesPermuted)->Parton(3) ).M(); 
	if (mtop < this->GetParameter(16)->GetLowerLimit())
		mtop = this->GetParameter(16)->GetLowerLimit(); 
	else if (mtop > this->GetParameter(16)->GetUpperLimit())
		mtop = this->GetParameter(16)->GetUpperLimit(); 
	values.push_back( mtop);

	// return the vector
	return values;
}

// --------------------------------------------------------- 
std::vector<double> KLFitter::LikelihoodTopLeptonJets::GetNeutrinoPzSolutions() {
	return this->CalculateNeutrinoPzSolutions();
}

// --------------------------------------------------------- 
std::vector<double> KLFitter::LikelihoodTopLeptonJets::CalculateNeutrinoPzSolutions(TLorentzVector * additionalParticle)
{
	std::vector<double> pz;

	KLFitter::PhysicsConstants constants;
	// electron mass
	double mE = 0.;

	double px_c = 0.0;
	double py_c = 0.0;
	double pz_c = 0.0;
	double Ec = 0.0; 

	if (fTypeLepton == 1)
		{
			px_c = (*fParticlesPermuted)->Electron(0)->Px();
			py_c = (*fParticlesPermuted)->Electron(0)->Py();
			pz_c = (*fParticlesPermuted)->Electron(0)->Pz();
			Ec = (*fParticlesPermuted)->Electron(0)->E();
		}
	else if (fTypeLepton == 2)
		{
			px_c = (*fParticlesPermuted)->Muon(0)->Px();
			py_c = (*fParticlesPermuted)->Muon(0)->Py();
			pz_c = (*fParticlesPermuted)->Muon(0)->Pz();
			Ec = (*fParticlesPermuted)->Muon(0)->E();
		}

	// add additional particle to "charged lepton" 4-vector
	if (additionalParticle) {
		px_c += additionalParticle->Px();
		py_c += additionalParticle->Py();
		pz_c += additionalParticle->Pz();
		Ec += additionalParticle->E();
	}

	double px_nu = ETmiss_x;
	double py_nu = ETmiss_y;
	double alpha = constants.MassW()*constants.MassW() - mE*mE + 2*(px_c*px_nu + py_c*py_nu);

	double a = pz_c*pz_c - Ec*Ec;
	double b = alpha * pz_c;
	double c = - Ec*Ec * (px_nu*px_nu + py_nu*py_nu) + alpha*alpha/4.;

	double discriminant = b*b - 4*a*c;
	if (discriminant < 0.)
		return pz;

	double pz_offset = - b / (2*a);

	double squareRoot = sqrt(discriminant);
	if (squareRoot < 1.e-6)
		pz.push_back(pz_offset);
	else {
		pz.push_back(pz_offset + squareRoot / (2*a));
		pz.push_back(pz_offset - squareRoot / (2*a));
	}

	return pz;
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTopLeptonJets::DeleteModelParticles()
{
	// check if model particles and lorentz vector container exist and delete
	if (fParticlesModel) {
		delete fParticlesModel; 
		fParticlesModel = 0;
	}

	if (fLorentzVectorContainer)
		{
			while (!fLorentzVectorContainer->empty())
				{
					TLorentzVector * lv = fLorentzVectorContainer->front(); 
					fLorentzVectorContainer->erase(fLorentzVectorContainer ->begin() ); 
					delete lv; 
				}
			delete fLorentzVectorContainer; 
			fLorentzVectorContainer = 0;
		}

	// no error
	return 1;
}
// --------------------------------------------------------- 
