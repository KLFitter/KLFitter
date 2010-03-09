#include "LikelihoodTTHElectron.h" 
#include <iostream> 
#include <BAT/BCMath.h> 

// --------------------------------------------------------- 
KLFitter::LikelihoodTTHElectron::LikelihoodTTHElectron()
{
	// define model particles 
	this -> DefineModelParticles(); 

	// define parameters 
	this -> DefineParameters(); 
}

// --------------------------------------------------------- 
KLFitter::LikelihoodTTHElectron::~LikelihoodTTHElectron()
{
	while (!fLorentzVectorContainer -> empty())
		{
			TLorentzVector * lv = fLorentzVectorContainer -> front(); 
			fLorentzVectorContainer -> erase(fLorentzVectorContainer  -> begin() ); 
			delete lv; 
		}
	delete fLorentzVectorContainer; 

	if (fParticlesModel)
		delete fParticlesModel; 
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTTHElectron::DefineModelParticles()
{
	// create the particles of the model 
	fParticlesModel = new KLFitter::Particles(); 

	fLorentzVectorContainer = new std::vector <TLorentzVector *>(0); 

	TLorentzVector * lv = 0; 

	// add model particles 
	fParticlesModel -> AddParticle(lv = new TLorentzVector(),
																 KLFitter::Particles::kParton,
 																 "hadronic b quark",
																 1.0); 
	fLorentzVectorContainer -> push_back(lv); 

	fParticlesModel -> AddParticle(lv = new TLorentzVector(),
																 KLFitter::Particles::kParton, 
																 "leptonic b quark",
																 1.0); 
	fLorentzVectorContainer -> push_back(lv); 	

	fParticlesModel -> AddParticle(lv = new TLorentzVector(),
																 KLFitter::Particles::kParton,
																 "light quark 1",
																 0.0); 
	fLorentzVectorContainer -> push_back(lv); 	

	fParticlesModel -> AddParticle(lv = new TLorentzVector(),
																 KLFitter::Particles::kParton,
																 "light quark 2",
																 0.0); 
	fLorentzVectorContainer -> push_back(lv); 
	
	fParticlesModel -> AddParticle(lv = new TLorentzVector(),
																 KLFitter::Particles::kParton, 
																 "Higgs b quark 1",
																 1.0); 
	fLorentzVectorContainer -> push_back(lv); 	

	fParticlesModel -> AddParticle(lv = new TLorentzVector(),
																 KLFitter::Particles::kParton, 
																 "Higgs b quark 2",
																 1.0); 
	fLorentzVectorContainer -> push_back(lv); 	

	fParticlesModel -> AddParticle(lv = new TLorentzVector(),
																 KLFitter::Particles::kElectron,
																 "electron"); 
	fLorentzVectorContainer -> push_back(lv); 

	fParticlesModel -> AddParticle(lv = new TLorentzVector(),
																 KLFitter::Particles::kNeutrino, 
																 "neutrino"); 
	fLorentzVectorContainer -> push_back(lv); 

	fParticlesModel -> AddParticle(lv = new TLorentzVector(),
																 KLFitter::Particles::kBoson, 
																 "hadronic W"); 
	fLorentzVectorContainer -> push_back(lv); 

	fParticlesModel -> AddParticle(lv = new TLorentzVector(),
																 KLFitter::Particles::kBoson,
																 "leptonic W"); 
	fLorentzVectorContainer -> push_back(lv); 

	fParticlesModel -> AddParticle(lv = new TLorentzVector(),
																 KLFitter::Particles::kParton,
																 "hadronic top",
																 0.0); 
	fLorentzVectorContainer -> push_back(lv); 

	fParticlesModel -> AddParticle(lv = new TLorentzVector(),
																 KLFitter::Particles::kParton,
																 "leptonic top",
																 0.0); 
	fLorentzVectorContainer -> push_back(lv); 

	fParticlesModel -> AddParticle(lv = new TLorentzVector(),
																 KLFitter::Particles::kBoson,
																 "Higgs",
																 0.0); 
	fLorentzVectorContainer -> push_back(lv); 

	// no error 
	return 1;
}

// --------------------------------------------------------- 
void KLFitter::LikelihoodTTHElectron::DefineParameters()
{
	// add parameters of model
	this -> AddParameter("energy hadronic b",       fPhysicsConstants -> MassBottom(), 1000.0); 
	this -> AddParameter("energy leptonic b",       fPhysicsConstants -> MassBottom(), 1000.0); 
	this -> AddParameter("energy light quark 1",    0.0, 1000.0); 
	this -> AddParameter("energy light quark 2",    0.0, 1000.0); 
	this -> AddParameter("energy Higgs b quark 1",  0.0, 1000.0); 
	this -> AddParameter("energy Higgs b quark 2",  0.0, 1000.0); 
	this -> AddParameter("energy electron",         0.0, 1000.0); 
	this -> AddParameter("p_z neutrino",        -1000.0, 1000.0); 
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTTHElectron::CalculateLorentzVectors(std::vector <double> parameters)
{
	// variables
  double E, px, py, pz, pabs, scale; 
	double sumpx = 0; 
	double sumpy = 0; 
	TLorentzVector * vect; 

	// hadronic b quark 
	vect = (*fParticlesPermuted) -> Parton(0); 
	E  = parameters.at(0); 
	px = vect -> Px(); 
	py = vect -> Py(); 
	pz = vect -> Pz(); 
	pabs = sqrt(px*px + py*py + pz*pz); 
	scale = sqrt(E*E - fPhysicsConstants -> MassBottom()*fPhysicsConstants -> MassBottom()) / pabs; 
	fParticlesModel -> Parton(0) -> SetPxPyPzE(scale * px, scale * py, scale * pz, E); 

	sumpx += px; 
	sumpy += py; 

	// leptonic b quark 
	vect = (*fParticlesPermuted) -> Parton(1); 
	E  = parameters.at(1); 
	px = vect -> Px(); 
	py = vect -> Py(); 
	pz = vect -> Pz(); 
	pabs = sqrt(px*px + py*py + pz*pz); 
	scale = sqrt(E*E - fPhysicsConstants -> MassBottom()*fPhysicsConstants -> MassBottom()) / pabs; 
	fParticlesModel -> Parton(1) -> SetPxPyPzE(scale * px, scale * py, scale * pz, E); 

	sumpx += px; 
	sumpy += py; 

	// light quark 1 
	vect = (*fParticlesPermuted) -> Parton(2); 
	E  = parameters.at(2); 
	px = vect -> Px(); 
	py = vect -> Py(); 
	pz = vect -> Pz(); 
	pabs = sqrt(px*px + py*py + pz*pz); 
	scale = E / pabs; 
	fParticlesModel -> Parton(2) -> SetPxPyPzE(scale * px, scale * py, scale * pz, E); 

	sumpx += px; 
	sumpy += py; 

	// light quark 2 
	vect = (*fParticlesPermuted) -> Parton(3); 
	E  = parameters.at(3); 
	px = vect -> Px(); 
	py = vect -> Py(); 
	pz = vect -> Pz(); 
	pabs = sqrt(px*px + py*py + pz*pz); 
	scale = E / pabs;
	fParticlesModel -> Parton(3) -> SetPxPyPzE(scale * px, scale * py, scale * pz, E); 

	sumpx += px; 
	sumpy += py; 

	// Higgs b quark 1 
	vect = (*fParticlesPermuted) -> Parton(4); 
	E  = parameters.at(4); 
	px = vect -> Px(); 
	py = vect -> Py(); 
	pz = vect -> Pz(); 
	pabs = sqrt(px*px + py*py + pz*pz); 
	scale = E / pabs; 
	fParticlesModel -> Parton(4) -> SetPxPyPzE(scale * px, scale * py, scale * pz, E); 

	sumpx += px; 
	sumpy += py; 

	// Higgs b quark 2
	vect = (*fParticlesPermuted) -> Parton(5); 
	E  = parameters.at(5); 
	px = vect -> Px(); 
	py = vect -> Py(); 
	pz = vect -> Pz(); 
	pabs = sqrt(px*px + py*py + pz*pz); 
	scale = E / pabs; 
	fParticlesModel -> Parton(5) -> SetPxPyPzE(scale * px, scale * py, scale * pz, E); 

	sumpx += px; 
	sumpy += py; 

	// electron
	vect = (*fParticlesPermuted) -> Electron(0); 
	E  = parameters.at(6); 
	px = vect -> Px(); 
	py = vect -> Py(); 
	pz = vect -> Pz(); 
	pabs = sqrt(px*px + py*py + pz*pz); 
	scale = E / pabs; 
	fParticlesModel -> Electron(0) -> SetPxPyPzE(scale * px, scale * py, scale * pz, E); 

	sumpx += px; 
	sumpy += py; 

	// neutrino 
	px = - sumpx; 
	py = - sumpy; 
	pz = parameters.at(7); 
	E = sqrt(px*px + py*py + pz*pz);
	fParticlesModel -> Neutrino(0) -> SetPxPyPzE(px, py, pz, E); 

	// composite particles 

	// hadronic W 
	*(fParticlesModel -> Boson(0)) = *(fParticlesModel -> Parton(2)) + *(fParticlesModel -> Parton(3)); 

	// leptonic W 
	*(fParticlesModel -> Boson(1)) = *(fParticlesModel -> Electron(0)) + *(fParticlesModel -> Neutrino(0)); 

	// hadronic top 
	*(fParticlesModel -> Parton(6)) = *(fParticlesModel -> Boson(0)) + *(fParticlesModel -> Parton(0)); 

	// leptonic top 
	*(fParticlesModel -> Parton(7)) = *(fParticlesModel -> Boson(1)) + *(fParticlesModel -> Parton(1)); 

	// Higgs
	*(fParticlesModel -> Boson(2)) = *(fParticlesModel -> Parton(4)) + *(fParticlesModel -> Parton(5)); 

	// no error 
	return 1; 
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTTHElectron::Initialize()
{
	// error code 
	int err = 1; 

	// adjust parameter ranges 
	err *= this -> AdjustParameterRanges(); 

	// return error code 
	return err; 
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTTHElectron::RemoveInvariantParticlePermutations()
{
	// error code 
	int err = 1; 

	// remove the permutation from the second and the third jet
	KLFitter::Particles::ParticleType ptype = KLFitter::Particles::kParton;
	int indexArray_Jets[2] = {2, 3};
	std::vector<int> indexVector_Jets(indexArray_Jets, indexArray_Jets + sizeof(indexArray_Jets)/sizeof(int));
	err *= (*fPermutations) -> InvariantParticlePermutations(ptype, indexVector_Jets); 

	// remove the permutation from the two Higgs b jets
	ptype = KLFitter::Particles::kParton;
	indexArray_Jets[0] = 4;
	indexArray_Jets[1] = 5;
	std::vector<int> indexVector_Jets2(indexArray_Jets, indexArray_Jets + sizeof(indexArray_Jets)/sizeof(int));
	err *= (*fPermutations) -> InvariantParticlePermutations(ptype, indexVector_Jets2); 

	// remove the permutation from the muons
	ptype = KLFitter::Particles::kMuon;
	KLFitter::Particles * particles = (*fPermutations) -> Particles();
	int indexArray_Muons[particles->NMuons()];
	for (int iMuon = 0; iMuon < particles->NMuons(); iMuon++)
		indexArray_Muons[iMuon] = iMuon;
	std::vector<int> indexVector_Muons(indexArray_Muons, indexArray_Muons + sizeof(indexArray_Muons)/sizeof(int));
	err *= (*fPermutations) -> InvariantParticlePermutations(ptype, indexVector_Muons); 

	// return error code 
	return err; 
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTTHElectron::AdjustParameterRanges()
{
	// adjust limits 
	double nsigmas = 7.0; 
	double nsigmas_lepton = 2.0; 

	// parameter 0: energy of hadronic b quark 
	double E = (*fParticlesPermuted) -> Parton(0) -> E(); 
	double Emin = TMath::Max(fPhysicsConstants -> MassBottom(), E - nsigmas * sqrt(E)); 
	double Emax  = E + nsigmas * sqrt(E); 
	this -> SetParameterRange(0, Emin, Emax); 

	// parameter 1: energy of leptonic b quark 
	E = (*fParticlesPermuted) -> Parton(1) -> E(); 
	Emin = TMath::Max(fPhysicsConstants -> MassBottom(), E - nsigmas * sqrt(E)); 
	Emax  = E + nsigmas * sqrt(E); 
	this -> SetParameterRange(1, Emin, Emax); 

	// parameter 2: energy of light quark 1 
	E = (*fParticlesPermuted) -> Parton(2) -> E(); 
	Emin = TMath::Max(0.001, E - nsigmas * sqrt(E)); 
	Emax  = E + nsigmas * sqrt(E); 
	this -> SetParameterRange(2, Emin, Emax); 

	// parameter 3: energy of light quark2 
	E = (*fParticlesPermuted) -> Parton(3) -> E(); 
	Emin = TMath::Max(0.001, E - nsigmas * sqrt(E)); 
	Emax  = E + nsigmas * sqrt(E); 
	this -> SetParameterRange(3, Emin, Emax); 

	// parameter 4: energy of Higgs b quark 1
	E = (*fParticlesPermuted) -> Parton(4) -> E(); 
	Emin = TMath::Max(0.001, E - nsigmas * sqrt(E)); 
	Emax  = E + nsigmas * sqrt(4); 
	this -> SetParameterRange(4, Emin, Emax); 

	// parameter 5: energy of Higgs b quark 2
	E = (*fParticlesPermuted) -> Parton(5) -> E(); 
	Emin = TMath::Max(0.001, E - nsigmas * sqrt(E)); 
	Emax  = E + nsigmas * sqrt(E); 
	this -> SetParameterRange(5, Emin, Emax); 

	// parameter 2: energy of electron
	E = (*fParticlesPermuted) -> Electron(0) -> E(); 
	Emin = TMath::Max(0.001, E - nsigmas_lepton * sqrt(E)); 
	Emax  = E + nsigmas_lepton * sqrt(E); 
	this -> SetParameterRange(6, Emin, Emax); 

	// no error 
	return 1; 
}

// --------------------------------------------------------- 
double KLFitter::LikelihoodTTHElectron::LogLikelihood(std::vector <double> parameters)
{
	// calculate 4-vectors 
	this -> CalculateLorentzVectors(parameters); 

	// define log of likelihood 
	double logprob = 0.; 

	// energy resolution terms 
	logprob += log( (*fDetector) -> ResEnergyBJet( (*fParticlesPermuted) -> Parton(0) -> Eta() ) -> p( parameters.at(0), (*fParticlesPermuted) -> Parton(0) -> E()) ); 
	logprob += log( (*fDetector) -> ResEnergyBJet( (*fParticlesPermuted) -> Parton(1) -> Eta() ) -> p( parameters.at(1), (*fParticlesPermuted) -> Parton(1) -> E()) ); 
	logprob += log( (*fDetector) -> ResEnergyLightJet( (*fParticlesPermuted) -> Parton(2) -> Eta() ) -> p( parameters.at(2), (*fParticlesPermuted) -> Parton(2) -> E()) ); 
	logprob += log( (*fDetector) -> ResEnergyLightJet( (*fParticlesPermuted) -> Parton(3) -> Eta() ) -> p( parameters.at(3), (*fParticlesPermuted) -> Parton(3) -> E()) ); 
	logprob += log( (*fDetector) -> ResEnergyLightJet( (*fParticlesPermuted) -> Parton(4) -> Eta() ) -> p( parameters.at(4), (*fParticlesPermuted) -> Parton(4) -> E()) ); 
	logprob += log( (*fDetector) -> ResEnergyLightJet( (*fParticlesPermuted) -> Parton(5) -> Eta() ) -> p( parameters.at(5), (*fParticlesPermuted) -> Parton(5) -> E()) ); 
	logprob += log( (*fDetector) -> ResEnergyElectron( (*fParticlesPermuted) -> Electron(0) -> Eta() ) -> p( parameters.at(6), (*fParticlesPermuted) -> Electron(0) -> E()) ); 

	// Breit-Wigner of hadronically decaying W-boson
 	logprob += BCMath::LogBreitWignerRel( fParticlesModel -> Boson(0) -> M(), fPhysicsConstants -> MassW(), fPhysicsConstants -> GammaW()); 
	
	// Breit-Wigner of leptonically decaying W-boson
 	logprob += BCMath::LogBreitWignerRel( fParticlesModel -> Boson(1) -> M(), fPhysicsConstants -> MassW(), fPhysicsConstants -> GammaW()); 

	// Breit-Wigner of hadronically decaying top quark
	logprob += BCMath::LogBreitWignerRel( fParticlesModel -> Parton(6) -> M(), fPhysicsConstants -> MassTop(), fPhysicsConstants -> GammaTop()); 
	
	// Breit-Wigner of leptonically decaying top quark
	logprob += BCMath::LogBreitWignerRel( fParticlesModel -> Parton(7) -> M(), fPhysicsConstants -> MassTop(), fPhysicsConstants -> GammaTop()); 

	// return log of likelihood 
	return logprob; 
}

// --------------------------------------------------------- 

