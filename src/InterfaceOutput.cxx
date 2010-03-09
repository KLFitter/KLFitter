#include "InterfaceOutput.h" 
#include <iostream> 
#include <map> 

// --------------------------------------------------------- 
KLFitter::InterfaceOutput::InterfaceOutput()
{
	fFitter = 0; 
	fMatchingTool = 0; 
	fSelectionTool = 0; 
	fParticlesTruth = 0; 
	fParticlesModel = 0; 
	fParticlesMeasured = 0; 
	fParticlesSelected = 0; 
	fTreeTruth = 0; 
	fTreeModel = 0;
 	fTreeMeasured = 0; 
 	fTreeSelected = 0; 
	fTreeMatching = 0; 
	fTreeMap = 0; 
	fTreeVarBestPermutation = 0; 
	fTreeVarLogLikelihood = 0; 
	fTreeVarIntegral = 0;
 	fTreeVarEventProbability = 0; 
	fTreeVarMinuitStatus = 0; 
	fTreeVarEventNumber = 0; 
	fTreeVarParameters = new std::vector<double *>(0); 
	fTreeVarParameterErrors = new std::vector<double *>(0); 
	fTreeVarModel = new std::vector<double *>(0); 
	fTreeVarTruth = new std::vector<double *>(0); 
	fTreeVarMeasured = new std::vector<std::vector<double> *>(0); 
	fTreeVarSelected = new std::vector<std::vector<double> *>(0); 
	fTreeVarNMatchedPartons = new std::vector <int *>(0); 
	fTreeVarNMatchedElectrons = new std::vector <int *>(0); 
	fTreeVarNMatchedMuons = new std::vector <int *>(0); 
	fTreeVarNMatchedPhotons = new std::vector <int *>(0); 
	fTreeVarMatchedPartons = new std::vector <int *>(0); 
	fTreeVarMatchedElectrons = new std::vector <int *>(0); 
	fTreeVarMatchedMuons = new std::vector <int *>(0); 
	fTreeVarMatchedPhotons = new std::vector <int *>(0); 
	fTreeVarMapJets = 0; 
	fTreeVarMapElectrons =  0;
	fTreeVarMapMuons = 0; 
	fTreeVarMapPhotons = 0; 
	fEventWeight = 0.;
	fIsNotClassified = false;
	fIsRadTopProd = false;
	fIsHadTopRadDecay = false;
	fIsLepTopRadDecay = false;
	fIsHadWRadDecay = false;
	fIsLepWRadDecay = false;
}

// --------------------------------------------------------- 
KLFitter::InterfaceOutput::~InterfaceOutput()
{
	if (fTreeVarBestPermutation)
		delete [] fTreeVarBestPermutation; 

	if (fTreeVarLogLikelihood)
		delete [] fTreeVarLogLikelihood; 

	if (fTreeVarIntegral)
		delete [] fTreeVarIntegral; 

	if (fTreeVarEventProbability)
		delete [] fTreeVarEventProbability; 

	if (fTreeVarMinuitStatus)
		delete [] fTreeVarMinuitStatus; 

	while (!fTreeVarParameters -> empty())
		{
			double * d = fTreeVarParameters -> front(); 
			fTreeVarParameters -> erase(fTreeVarParameters -> begin()); 
			delete [] d; 
		}
	delete fTreeVarParameters; 

	while (!fTreeVarParameterErrors -> empty())
		{
			double * d = fTreeVarParameterErrors -> front(); 
			fTreeVarParameterErrors -> erase(fTreeVarParameterErrors -> begin()); 
			delete [] d; 
		}
	delete fTreeVarParameterErrors; 
	
	while (!fTreeVarModel -> empty())
		{
			double * d = fTreeVarModel -> front(); 
			fTreeVarModel -> erase(fTreeVarModel -> begin()); 
			delete [] d; 
		}
	delete fTreeVarModel; 
	
	while (!fTreeVarTruth -> empty())
		{
			double * d = fTreeVarTruth -> front(); 
			fTreeVarTruth -> erase(fTreeVarTruth -> begin()); 
			delete d; 
		}
	delete fTreeVarTruth; 
	
	while (!fTreeVarMeasured -> empty())
		{
			std::vector<double> * d = fTreeVarMeasured -> front(); 
			fTreeVarMeasured -> erase(fTreeVarMeasured -> begin()); 
			delete d; 
		}
	delete fTreeVarMeasured;

	while (!fTreeVarSelected -> empty())
		{
			std::vector<double> * d = fTreeVarSelected -> front(); 
			fTreeVarSelected -> erase(fTreeVarSelected -> begin()); 
			delete d; 
		}
	delete fTreeVarSelected;
	
	while (!fTreeVarNMatchedPartons -> empty())
		{
			int * i = fTreeVarNMatchedPartons -> front(); 
			fTreeVarNMatchedPartons -> erase(fTreeVarNMatchedPartons -> begin()); 
			delete i; 
		}
	delete fTreeVarNMatchedPartons;

	while (!fTreeVarNMatchedElectrons -> empty())
		{
			int * i = fTreeVarNMatchedElectrons -> front(); 
			fTreeVarNMatchedElectrons -> erase(fTreeVarNMatchedElectrons -> begin()); 
			delete i; 
		}
	delete fTreeVarNMatchedElectrons;

	while (!fTreeVarNMatchedMuons -> empty())
		{
			int * i = fTreeVarNMatchedMuons -> front(); 
			fTreeVarNMatchedMuons -> erase(fTreeVarNMatchedMuons -> begin()); 
			delete i; 
		}
	delete fTreeVarNMatchedMuons;

	while (!fTreeVarNMatchedPhotons -> empty())
		{
			int * i = fTreeVarNMatchedPhotons -> front(); 
			fTreeVarNMatchedPhotons -> erase(fTreeVarNMatchedPhotons -> begin()); 
			delete i; 
		}
	delete fTreeVarNMatchedPhotons;

	while (!fTreeVarMatchedPartons -> empty())
		{
			int * i = fTreeVarMatchedPartons -> front(); 
			fTreeVarMatchedPartons -> erase(fTreeVarMatchedPartons -> begin()); 
			delete i; 
		}
	delete fTreeVarMatchedPartons;

	while (!fTreeVarMatchedElectrons -> empty())
		{
			int * i = fTreeVarMatchedElectrons -> front(); 
			fTreeVarMatchedElectrons -> erase(fTreeVarMatchedElectrons -> begin()); 
			delete i; 
		}
	delete fTreeVarMatchedElectrons;

	while (!fTreeVarMatchedMuons -> empty())
		{
			int * i = fTreeVarMatchedMuons -> front(); 
			fTreeVarMatchedMuons -> erase(fTreeVarMatchedMuons -> begin()); 
			delete i; 
		}
	delete fTreeVarMatchedMuons;

	while (!fTreeVarMatchedPhotons -> empty())
		{
			int * i = fTreeVarMatchedPhotons -> front(); 
			fTreeVarMatchedPhotons -> erase(fTreeVarMatchedPhotons -> begin()); 
			delete i; 
		}
	delete fTreeVarMatchedPhotons;

	if (fTreeVarMapJets)
		delete [] fTreeVarMapJets; 

	if (fTreeVarMapElectrons)
		delete [] fTreeVarMapElectrons; 

	if (fTreeVarMapMuons)
		delete [] fTreeVarMapMuons; 

	if (fTreeVarMapPhotons)
		delete [] fTreeVarMapPhotons; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput::SetFitter(KLFitter::Fitter * fitter)
{
	// check if fitter exists 
	if (!fitter)
		{
			std::cout << "KLFitter::InterfaceOutput::SetFitter(). Fitter does not exist." << std::endl; 
			return 0; 
		}

	// set pointers to pointer 
	fParticlesModel = fitter -> Likelihood() -> PParticlesModel(); 
	fParticlesSelected = fitter -> PParticles(); 

	// set fitter
	fFitter = fitter;

	// no error 
	return 1; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput::SetMatchingTool(KLFitter::MatchingTool * matchingtool)
{
	// check if fitter exists 
	if (!matchingtool)
		{
			std::cout << "KLFitter::InterfaceOutput::SetMatchingTool(). Matching tool does not exist." << std::endl; 
			return 0; 
		}

	// set pointer to matching tool 
	fMatchingTool = matchingtool; 

	// no error 
	return 1;
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput::SetSelectionTool(KLFitter::SelectionTool * selectiontool)
{
	// check if fitter exists 
	if (!selectiontool)
		{
			std::cout << "KLFitter::InterfaceOutput::SetSelectionTool(). Selection tool does not exist." << std::endl; 
			return 0; 
		}

	// set pointer to matching tool 
	fSelectionTool = selectiontool; 

	// no error 
	return 1;
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput::OpenRootFile(const char * filename, Option_t * opt)
{
	// define error code 
	int err = 1; 

	// open file 
	err *= KLFitter::InterfaceRoot::OpenRootFile(filename, opt); 

	// return error code 
	return err; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput::CloseRootFile()
{
	// define error code 
	int err = 1; 

	// check if file exists 
	if (!fRootFile)
		return 0; 
	
	// check if file is open 
	if (!fRootFile -> IsOpen())
		return 0; 

	if (fTreeTruth)
		fRootFile -> WriteTObject(fTreeTruth); 
	
	if (fTreeModel)
		fRootFile -> WriteTObject(fTreeModel); 

	if (fTreeMeasured)
		fRootFile -> WriteTObject(fTreeMeasured); 

	if (fTreeSelected)
		fRootFile -> WriteTObject(fTreeSelected); 

	if (fTreeMatching)
		fRootFile -> WriteTObject(fTreeMatching); 

	if (fTreeMap)
		fRootFile -> WriteTObject(fTreeMap); 

	// close file 
	KLFitter::InterfaceRoot::CloseRootFile(); 

	// return error code 
	return err; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput::CreateTreeModel()
{
	// check if particles exist
	if (!fParticlesModel)
		{
			std::cout << "KLFitter::InterfaceOutput::CreateTreeModel(). Particles do not exist." << std::endl; 
			return 0; 
		}

	// delete old tree if necessary 
	if (fTreeModel)
		delete fTreeModel; 

	// get number of permutations 
	int nperm = fFitter -> Permutations() -> NPermutations(); 

	// create new tree 
	fTreeModel = new TTree("TreeModel", "TreeModel"); 

	// reset variables
	fTreeVarNPermutations = nperm; 
	fTreeVarNBTags=0;
	fTreeVarBestPermutation = new int[nperm]; 
	fTreeVarLogLikelihood = new double[nperm]; 
	fTreeVarIntegral = new double[nperm]; 
	fTreeVarEventProbability = new double[nperm]; 
	fTreeVarMinuitStatus = new double[nperm]; 
	for (int i = 0; i < nperm; ++i)
		{
			fTreeVarBestPermutation[i] = 0; 
			fTreeVarLogLikelihood[i] = 1e99; 
			fTreeVarIntegral[i] = -1; 
			fTreeVarEventProbability[i] = 0.; 
			fTreeVarMinuitStatus[i] = 0; 
		}

	// set branches for event variables
	fTreeModel -> Branch("EventNumber", &fTreeVarEventNumber, "EventNumber/I"); 
	fTreeModel -> Branch("N_permutations", &fTreeVarNPermutations, "N_permutations/I"); 
	fTreeModel -> Branch("N_btags", &fTreeVarNBTags, "N_btags/I"); 
	fTreeModel -> Branch("best_permutation", fTreeVarBestPermutation, "best_permutation[N_permutations]/I"); 
	fTreeModel -> Branch("LogLikelihood", fTreeVarLogLikelihood, "LogLikelihood[N_permutations]/D"); 
	fTreeModel -> Branch("Integral", fTreeVarIntegral, "Integral[N_permutations]/D"); 
	fTreeModel -> Branch("EventProbability", fTreeVarEventProbability, "EventProbability[N_permutations]/D"); 
	fTreeModel -> Branch("MinuitStatus", fTreeVarMinuitStatus, "MinuitStatus[N_permutations]/D"); 

	// loop over all parameters 
	for (int i = 0; i < fFitter -> Likelihood() -> NParameters(); ++i)
		{
			double * par = new double[nperm]; 
			double * parerr = new double[nperm]; 
			fTreeVarParameters -> push_back(par); 
			fTreeVarParameterErrors -> push_back(parerr); 
			fTreeModel -> Branch( this -> ModifyString( "par_" + fFitter -> Likelihood() -> GetParameter(i) -> GetName() ).data(), par, this -> ModifyString( "par_" + fFitter -> Likelihood() -> GetParameter(i) -> GetName() + "[N_permutations]/D").data()); 
			fTreeModel -> Branch( this -> ModifyString( "parerr_" + fFitter -> Likelihood() -> GetParameter(i) -> GetName() ).data(), parerr, this -> ModifyString( "parerr_" + fFitter -> Likelihood() -> GetParameter(i) -> GetName() + "[N_permutations]/D").data()); 
		}

 	// loop over all particle type 
	for (KLFitter::Particles::ParticleType itype = KLFitter::Particles::kParton; itype <= KLFitter::Particles::kPhoton; ++itype)
		{
			// get particle container
			std::vector <TLorentzVector *> * momcontainer = (*fParticlesModel) -> ParticleContainer(itype); 
			std::vector <std::string> * namecontainer = (*fParticlesModel) -> ParticleNameContainer(itype); 

			// get number of particles in container 
			int n = int(momcontainer -> size()); 			
			
			// loop over particles 
			for (int i = 0; i < n; ++i)
				{
					// get name 
					std::string name(namecontainer -> at(i)); 

					// create new pointer to variables 
					double * E = new double[nperm]; 
					double * px = new double[nperm]; 
					double * py = new double[nperm]; 
					double * pz = new double[nperm]; 
					double * m = new double[nperm]; 
					double * pt = new double[nperm]; 
					double * eta = new double[nperm]; 
					double * phi = new double[nperm]; 
					double * btag = new double[nperm]; 
					double * index = new double[nperm]; 

					// initialize variables 
					for (int j = 0; j < nperm; ++j)
						{
							E[j] = 0.0; 
							px[j] = 0.0; 
							py[j] = 0.0; 
							pz[j] = 0.0; 
							m[j] = 0.0; 
							pt[j] = 0.0; 
							eta[j] = 0.0; 
							phi[j] = 0.0; 
							btag[j] = 0.0; 
							index[j] = -1.0; 
						}

					// add variables to vector 
					fTreeVarModel -> push_back(E); 
					fTreeVarModel -> push_back(px); 
					fTreeVarModel -> push_back(py); 
					fTreeVarModel -> push_back(pz); 
					fTreeVarModel -> push_back(m); 
					fTreeVarModel -> push_back(pt); 
					fTreeVarModel -> push_back(eta); 
					fTreeVarModel -> push_back(phi); 
					fTreeVarModel -> push_back(btag); 
					fTreeVarModel -> push_back(index); 

					// create new branches 					
					fTreeModel -> Branch(this -> ModifyString((name+"_E")).data(), E, this -> ModifyString((name+"_E[N_permutations]/D")).data()); 
					fTreeModel -> Branch(this -> ModifyString((name+"_px")).data(), px, this -> ModifyString((name+"_px[N_permutations]/D")).data());
					fTreeModel -> Branch(this -> ModifyString((name+"_py")).data(), py, this -> ModifyString((name+"_py[N_permutations]/D")).data());
					fTreeModel -> Branch(this -> ModifyString((name+"_pz")).data(), pz, this -> ModifyString((name+"_pz[N_permutations]/D")).data());					
					fTreeModel -> Branch(this -> ModifyString((name+"_m")).data(), m, this -> ModifyString((name+"_m[N_permutations]/D")).data());					
					fTreeModel -> Branch(this -> ModifyString((name+"_pt")).data(), pt, this -> ModifyString((name+"_pt[N_permutations]/D")).data());					
					fTreeModel -> Branch(this -> ModifyString((name+"_eta")).data(), eta, this -> ModifyString((name+"_eta[N_permutations]/D")).data());					
					fTreeModel -> Branch(this -> ModifyString((name+"_phi")).data(), phi, this -> ModifyString((name+"_phi[N_permutations]/D")).data());					
					if (itype == KLFitter::Particles::kParton && (*fParticlesModel)->JetIndex(i)>=0) {
						fTreeModel -> Branch(this -> ModifyString((name+"_btag")).data(), btag, this -> ModifyString((name+"_btag[N_permutations]/D")).data()); 
						fTreeModel -> Branch(this -> ModifyString((name+"_index")).data(), index, this -> ModifyString((name+"_index[N_permutations]/D")).data()); 
					}					
					if (itype == KLFitter::Particles::kElectron && (*fParticlesModel)->ElectronIndex(i)>=0) 
						fTreeModel -> Branch(this -> ModifyString((name+"_index")).data(), index, this -> ModifyString((name+"_index[N_permutations]/D")).data()); 
					if (itype == KLFitter::Particles::kMuon && (*fParticlesModel)->MuonIndex(i)>=0) 
						fTreeModel -> Branch(this -> ModifyString((name+"_index")).data(), index, this -> ModifyString((name+"_index[N_permutations]/D")).data()); 
					if (itype == KLFitter::Particles::kPhoton && (*fParticlesModel)->PhotonIndex(i)>=0) 
						fTreeModel -> Branch(this -> ModifyString((name+"_index")).data(), index, this -> ModifyString((name+"_index[N_permutations]/D")).data()); 
				}
		}

	// no error 
	return 1; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput::CreateTreeMeasured()
{
	// check if particles exist
	if (!fParticlesMeasured)
		{
			std::cout << "KLFitter::InterfaceOutput::CreateTreeMeasured(). Particles do not exist." << std::endl; 
			return 0; 
		}

	// delete old tree if necessary 
	if (fTreeMeasured)
		delete fTreeMeasured; 

	// create new tree 
	fTreeMeasured = new TTree("TreeMeasured", "TreeMeasured"); 
	
	// event weight branch
	fTreeMeasured -> Branch("Weight", &fEventWeight, "Weight/D"); 
	
	fTreeMeasured -> Branch("N_jets", &fTreeVarNPartonsMeasured, "N_jets/I"); 
	fTreeMeasured -> Branch("N_electrons", &fTreeVarNElectronsMeasured, "N_electrons/I"); 
	fTreeMeasured -> Branch("N_muons", &fTreeVarNMuonsMeasured, "N_muons/I"); 
	fTreeMeasured -> Branch("N_photons", &fTreeVarNPhotonsMeasured, "N_photons/I"); 

	// loop over all particle type 
	for (KLFitter::Particles::ParticleType itype = KLFitter::Particles::kParton; itype <= KLFitter::Particles::kPhoton; ++itype)
		{
			// get the name of the branch
			std::string name = TreeMeasuredBranchName(itype);
			if (name == "")
				continue;

			// create new pointer to variables 
			std::vector<double> * vec_E   = new std::vector<double>(0); 
			std::vector<double> * vec_px  = new std::vector<double>(0); 
			std::vector<double> * vec_py  = new std::vector<double>(0); 
			std::vector<double> * vec_pz  = new std::vector<double>(0); 
			std::vector<double> * vec_m   = new std::vector<double>(0); 
			std::vector<double> * vec_pt  = new std::vector<double>(0); 
			std::vector<double> * vec_eta = new std::vector<double>(0); 
			std::vector<double> * vec_phi = new std::vector<double>(0); 

			// add variables to vector 
			fTreeVarMeasured -> push_back(vec_E); 
			fTreeVarMeasured -> push_back(vec_px); 
			fTreeVarMeasured -> push_back(vec_py); 
			fTreeVarMeasured -> push_back(vec_pz); 
			fTreeVarMeasured -> push_back(vec_m); 
			fTreeVarMeasured -> push_back(vec_pt); 
			fTreeVarMeasured -> push_back(vec_eta); 
			fTreeVarMeasured -> push_back(vec_phi); 

			// create new branches 					
			fTreeMeasured -> Branch((name+"_E").data(),   vec_E  ); 
			fTreeMeasured -> Branch((name+"_px").data(),  vec_px );
			fTreeMeasured -> Branch((name+"_py").data(),  vec_py );
			fTreeMeasured -> Branch((name+"_pz").data(),  vec_pz );
			fTreeMeasured -> Branch((name+"_m").data(),   vec_m  ); 
			fTreeMeasured -> Branch((name+"_pt").data(),  vec_pt );
			fTreeMeasured -> Branch((name+"_eta").data(), vec_eta);
			fTreeMeasured -> Branch((name+"_phi").data(), vec_phi);
		}

	// no error 
	return 1; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput::CreateTreeSelected()
{
	// check if particles exist
	if (!fParticlesSelected)
		{
			std::cout << "KLFitter::InterfaceOutput::CreateTreeSelected(). Particles do not exist." << std::endl; 
			return 0; 
		}

	// delete old tree if necessary 
	if (fTreeSelected)
		delete fTreeSelected; 

	// create new tree 
	fTreeSelected = new TTree("TreeSelected", "TreeSelected"); 

	fTreeSelected -> Branch("N_jets", &fTreeVarNPartonsSelected, "N_jets/I"); 
	fTreeSelected -> Branch("N_electrons", &fTreeVarNElectronsSelected, "N_electrons/I"); 
	fTreeSelected -> Branch("N_muons", &fTreeVarNMuonsSelected, "N_muons/I"); 
	fTreeSelected -> Branch("N_photons", &fTreeVarNPhotonsSelected, "N_photons/I"); 

	// loop over all particle type 
	for (KLFitter::Particles::ParticleType itype = KLFitter::Particles::kParton; itype <= KLFitter::Particles::kPhoton; ++itype)
		{
			// get the name of the branch
			std::string name = TreeMeasuredBranchName(itype);
			if (name == "")
				continue;

			// create new pointer to variables 
			std::vector<double> * vec_E    = new std::vector<double>(0); 
			std::vector<double> * vec_px   = new std::vector<double>(0); 
			std::vector<double> * vec_py   = new std::vector<double>(0); 
			std::vector<double> * vec_pz   = new std::vector<double>(0); 
			std::vector<double> * vec_m    = new std::vector<double>(0); 
			std::vector<double> * vec_pt   = new std::vector<double>(0); 
			std::vector<double> * vec_eta  = new std::vector<double>(0); 
			std::vector<double> * vec_phi  = new std::vector<double>(0); 
			std::vector<double> * vec_btag = new std::vector<double>(0); 

			// add variables to vector 
			fTreeVarSelected -> push_back(vec_E); 
			fTreeVarSelected -> push_back(vec_px); 
			fTreeVarSelected -> push_back(vec_py); 
			fTreeVarSelected -> push_back(vec_pz); 
			fTreeVarSelected -> push_back(vec_m); 
			fTreeVarSelected -> push_back(vec_pt); 
			fTreeVarSelected -> push_back(vec_eta); 
			fTreeVarSelected -> push_back(vec_phi); 
			fTreeVarSelected -> push_back(vec_btag); 

			// create new branches 					
			fTreeSelected -> Branch((name+"_E").data(),   vec_E  ); 
			fTreeSelected -> Branch((name+"_px").data(),  vec_px );
			fTreeSelected -> Branch((name+"_py").data(),  vec_py );
			fTreeSelected -> Branch((name+"_pz").data(),  vec_pz );
			fTreeSelected -> Branch((name+"_m").data(),   vec_m  ); 
			fTreeSelected -> Branch((name+"_pt").data(),  vec_pt );
			fTreeSelected -> Branch((name+"_eta").data(), vec_eta);
			fTreeSelected -> Branch((name+"_phi").data(), vec_phi);
			if (itype == KLFitter::Particles::kParton)
				fTreeSelected -> Branch((name+"_btag").data(), vec_btag);
		}

	// no error 
	return 1; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput::CreateTreeTruth()
{
	// check if particles exist
	if (!fParticlesTruth)
		{
			std::cout << "KLFitter::InterfaceOutput::CreateTreeTruth(). Particles do not exist." << std::endl; 
			return 0; 
		}

	// delete old tree if necessary 
	if (fTreeTruth)
		delete fTreeTruth; 

	// create new tree 
	fTreeTruth = new TTree("TreeTruth", "TreeTruth"); 

	// loop over all particle type 
	for (KLFitter::Particles::ParticleType itype = KLFitter::Particles::kParton; itype <= KLFitter::Particles::kPhoton; ++itype)
		{
			// get particle container
			std::vector <TLorentzVector *> * momcontainer = (*fParticlesTruth) -> ParticleContainer(itype); 
			std::vector <std::string> * namecontainer = (*fParticlesTruth) -> ParticleNameContainer(itype); 

			// get number of particles in container 
			int n = int(momcontainer -> size()); 			

			// loop over particles 
			for (int i = 0; i < n; ++i)
				{
					// get name 
					std::string name(namecontainer -> at(i)); 

					// create new pointer to variables 
					double * E = new double(0.0); 
					double * px = new double(0.0); 
					double * py = new double(0.0); 
					double * pz = new double(0.0); 
					double * m = new double(0.0); 
					double * pt = new double(0.0); 
					double * eta = new double(0.0); 
					double * phi = new double(0.0); 

					// add variables to vector 
					fTreeVarTruth -> push_back(E); 
					fTreeVarTruth -> push_back(px); 
					fTreeVarTruth -> push_back(py); 
					fTreeVarTruth -> push_back(pz); 
					fTreeVarTruth -> push_back(m); 
					fTreeVarTruth -> push_back(pt); 
					fTreeVarTruth -> push_back(eta); 
					fTreeVarTruth -> push_back(phi); 

					// create new branches 					
					fTreeTruth -> Branch(this -> ModifyString((name+"_E")).data(), E, this -> ModifyString((name+"_E/D")).data()); 
					fTreeTruth -> Branch(this -> ModifyString((name+"_px")).data(), px, this -> ModifyString((name+"_px/D")).data());
					fTreeTruth -> Branch(this -> ModifyString((name+"_py")).data(), py, this -> ModifyString((name+"_py/D")).data());
					fTreeTruth -> Branch(this -> ModifyString((name+"_pz")).data(), pz, this -> ModifyString((name+"_pz/D")).data());					
					fTreeTruth -> Branch(this -> ModifyString((name+"_m")).data(), m, this -> ModifyString((name+"_m/D")).data());					
					fTreeTruth -> Branch(this -> ModifyString((name+"_pt")).data(), pt, this -> ModifyString((name+"_pt/D")).data());					
					fTreeTruth -> Branch(this -> ModifyString((name+"_eta")).data(), eta, this -> ModifyString((name+"_eta/D")).data());					
					fTreeTruth -> Branch(this -> ModifyString((name+"_phi")).data(), phi, this -> ModifyString((name+"_phi/D")).data());					
				}
		}

	fTreeTruth->Branch("IsNotClassified",  &fIsNotClassified,  "IsNotClassified/B");
	fTreeTruth->Branch("IsRadTopProd",     &fIsRadTopProd,     "IsRadTopProd/B");
	fTreeTruth->Branch("IsHadTopRadDecay", &fIsHadTopRadDecay, "IsHadTopRadDecay/B");
	fTreeTruth->Branch("IsLepTopRadDecay", &fIsLepTopRadDecay, "IsLepTopRadDecay/B");
	fTreeTruth->Branch("IsHadWRadDecay",   &fIsHadWRadDecay,   "IsHadWRadDecay/B");
	fTreeTruth->Branch("IsLepWRadDecay",   &fIsLepWRadDecay,   "IsLepWRadDecay/B");

	// no error 
	return 1; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput::CreateTreeMatching()
{
	// check if particles exist
	if (!fMatchingTool)
		{
			std::cout << "KLFitter::InterfaceOutput::CreateTreeMatching(). Matching tool does not exist." << std::endl; 
			return 0; 
		}

	// delete old tree if necessary 
	if (fTreeMatching)
		delete fTreeMatching; 

	// create new tree 
	fTreeMatching = new TTree("TreeMatching", "TreeMatching"); 

	fTreeMatching -> Branch("N_jets", &fTreeVarNPartonsSelected, "N_jets/I"); 
	fTreeMatching -> Branch("N_electrons", &fTreeVarNElectronsSelected, "N_electrons/I"); 
	fTreeMatching -> Branch("N_muons", &fTreeVarNMuonsSelected, "N_muons/I"); 
	fTreeMatching -> Branch("N_photons", &fTreeVarNPhotonsSelected, "N_photons/I"); 

	// get parton container
	std::vector <std::string> * namecontainer = (*fParticlesTruth) -> ParticleNameContainer(KLFitter::Particles::kParton); 

	// get number of particles in container 
	int ntruth = int(namecontainer -> size()); 			

	// loop over particles 
	for (int i = 0; i < ntruth; ++i)
		{
			// get name 
			std::string name(namecontainer -> at(i)); 
			
			int * pi = new int(-1); 
			fTreeVarNMatchedPartons -> push_back(pi); 

			int nreco_partons_max = 10; 
			int * indices = new int[nreco_partons_max]; 
			for (int i = 0; i < nreco_partons_max; ++i)
				indices[i] = -1; 
			fTreeVarMatchedPartons -> push_back(indices); 

			// create new branches 
			fTreeMatching -> Branch( this -> ModifyString("Nmatches_"+name).data(), 
															 pi,
															 this -> ModifyString("Nmatches_"+name+"/I").data());
			fTreeMatching -> Branch( this -> ModifyString("Matches_"+name).data(), 
															 indices,
															 this -> ModifyString("Matches_"+name+"[N_jets]/I").data());
	}
	
	// get electron container
	namecontainer = (*fParticlesTruth) -> ParticleNameContainer(KLFitter::Particles::kElectron); 

	// get number of particles in container 
	int n = int(namecontainer -> size()); 			

	// loop over particles 
	for (int i = 0; i < n; ++i)
		{
			// get name 
			std::string name(namecontainer -> at(i)); 
			
			int * pi = new int(-1); 
			fTreeVarNMatchedElectrons -> push_back(pi); 

			int nreco_electrons_max = 10; 
			int * indices = new int[nreco_electrons_max]; 
			for (int i = 0; i < nreco_electrons_max; ++i)
				indices[i] = -1; 
			fTreeVarMatchedElectrons -> push_back(indices); 

			// create new branches 
			fTreeMatching -> Branch( this -> ModifyString("Nmatches_"+name).data(), 
															 pi,
															 this -> ModifyString("Nmatches_"+name+"/I").data());
			fTreeMatching -> Branch( this -> ModifyString("Matches_"+name).data(), 
															 indices,
															 this -> ModifyString("Matches_"+name+"[N_electrons]/I").data());
		}
	
	// get muon container
	namecontainer = (*fParticlesTruth) -> ParticleNameContainer(KLFitter::Particles::kMuon); 

	// get number of particles in container 
	n = int(namecontainer -> size()); 			

	// loop over particles 
	for (int i = 0; i < n; ++i)
		{
			// get name 
			std::string name(namecontainer -> at(i)); 
			
			int * pi = new int(-1); 
			fTreeVarNMatchedMuons -> push_back(pi); 

			int nreco_muons_max = 10; 
			int * indices = new int[nreco_muons_max]; 
			for (int i = 0; i < nreco_muons_max; ++i)
				indices[i] = -1; 
			fTreeVarMatchedMuons -> push_back(indices); 

			// create new branches 
			fTreeMatching -> Branch( this -> ModifyString("Nmatches_"+name).data(), 
															 pi,
															 this -> ModifyString("Nmatches_"+name+"/I").data());
			fTreeMatching -> Branch( this -> ModifyString("Matches_"+name).data(), 
															 indices,
															 this -> ModifyString("Matches_"+name+"[N_muons]/I").data());
		}	

	// get photon container
	namecontainer = (*fParticlesTruth) -> ParticleNameContainer(KLFitter::Particles::kPhoton); 

	// get number of particles in container 
	n = int(namecontainer -> size()); 			

	// loop over particles 
	for (int i = 0; i < n; ++i)
		{
			// get name 
			std::string name(namecontainer -> at(i)); 
			
			int * pi = new int(-1); 
			fTreeVarNMatchedPhotons -> push_back(pi); 

			int nreco_photons_max = 10; 
			int * indices = new int[nreco_photons_max]; 
			for (int i = 0; i < nreco_photons_max; ++i)
				indices[i] = -1; 
			fTreeVarMatchedPhotons -> push_back(indices); 

			// create new branches 
			fTreeMatching -> Branch( this -> ModifyString("Nmatches_"+name).data(), 
															 pi,
															 this -> ModifyString("Nmatches_"+name+"/I").data());
			fTreeMatching -> Branch( this -> ModifyString("Matches_"+name).data(), 
															 indices,
															 this -> ModifyString("Matches_"+name+"[N_photons]/I").data());
		}	

	// no error 
	return 1; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput::CreateTreeMap()
{
	// check if particles exist
	if (!fSelectionTool)
		{
			std::cout << "KLFitter::InterfaceOutput::CreateTreeMap(). Selection tool does not exist." << std::endl; 
			return 0; 
		}

	// delete old tree if necessary 
	if (fTreeMap)
		delete fTreeMap; 

	// reset variables
	fTreeVarMapJets = new int[fTreeVarNPartonsSelected]; 
	fTreeVarMapElectrons = new int[fTreeVarNElectronsSelected]; 
	fTreeVarMapMuons = new int[fTreeVarNMuonsSelected]; 
	fTreeVarMapPhotons = new int[fTreeVarNPhotonsSelected]; 

	for (int i = 0; i < fTreeVarNPartonsSelected; ++i)
		fTreeVarMapJets[i] = -1; 
	for (int i = 0; i < fTreeVarNElectronsSelected; ++i)
		fTreeVarMapElectrons[i] = -1; 
	for (int i = 0; i < fTreeVarNMuonsSelected; ++i)
		fTreeVarMapMuons[i] = -1; 
	for (int i = 0; i < fTreeVarNPhotonsSelected; ++i)
		fTreeVarMapPhotons[i] = -1; 

	// create new tree 
	fTreeMap = new TTree("TreeMap", "TreeMap"); 

	fTreeMap -> Branch("N_jets", &fTreeVarNPartonsSelected, "N_jets/I"); 
	fTreeMap -> Branch("N_electrons", &fTreeVarNElectronsSelected, "N_electrons/I"); 
	fTreeMap -> Branch("N_muons", &fTreeVarNMuonsSelected, "N_muons/I"); 
	fTreeMap -> Branch("N_photons", &fTreeVarNPhotonsSelected, "N_photons/I"); 

	// set branches for event variables
	fTreeMap -> Branch("EventNumber", &fTreeVarEventNumber, "EventNumber/I"); 
	fTreeMap -> Branch("Index_jet", fTreeVarMapJets, "Index_jet[N_jets]/I"); 
	fTreeMap -> Branch("Index_electron", fTreeVarMapElectrons, "Index_electron[N_electrons]/I");
	fTreeMap -> Branch("Index_muon", fTreeVarMapMuons, "Index_muon[N_muons]/I");
	fTreeMap -> Branch("Index_photon", fTreeVarMapPhotons, "Index_photon[N_photons]/I");

	// no error 
	return 1; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput::CreateTrees()
{
	// error code 
	int err = 1; 
	
	// create tree for truth particles 
	if (fParticlesTruth)
		err *= this -> CreateTreeTruth(); 

	// create tree for measured particles 
	err *= this -> CreateTreeMeasured(); 

	// create tree for selected particles 
	err *= this -> CreateTreeSelected(); 

	// create tree for model particles 
	err *= this -> CreateTreeModel(); 

	// create tree for matching informations
	if (fMatchingTool)
		err *= this -> CreateTreeMatching(); 

	if (fSelectionTool)
		err *= this -> CreateTreeMap(); 

	// return error code 
	return err; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput::FillTreeModelPermutation()
{
	// check tree 
	if (!fTreeModel)
		this -> CreateTreeModel(); 

	// initialize counter 
	int counter = 0; 

	// create new permutation table for writing out the index
	int npartons = (*fParticlesSelected) -> NPartons(); 
	std::vector < std::vector<int> * > * table_partons = new std::vector < std::vector<int> * >(0); 
	fFitter->Permutations()-> CreateSubTable(npartons, table_partons); 

	// get permutation index
	int pindex =  fFitter -> Permutations() -> PermutationIndex(); 
	
	// fill event variables
	fTreeVarNPermutations = fFitter -> Permutations() -> NPermutations(); 
	fTreeVarNBTags = (*fParticlesModel)->NBTags(); 
	fTreeVarLogLikelihood[pindex] = fFitter -> Likelihood() -> LogLikelihood( fFitter -> Likelihood() -> GetBestFitParameters() ); 
	fTreeVarMinuitStatus[pindex] = fFitter -> MinuitStatus(); 
	fTreeVarIntegral[pindex] = fFitter -> Likelihood() -> GetNormalization(); 
	fTreeVarEventProbability[pindex] = exp( fFitter -> Likelihood() -> LogEventProbability() ); 

	// check event probability for NaN
	if (fTreeVarLogLikelihood[pindex] != fTreeVarLogLikelihood[pindex])
		fTreeVarLogLikelihood[pindex] = double(pindex) * (-1e10); 

	if (fTreeVarEventProbability[pindex] != fTreeVarEventProbability[pindex])
		fTreeVarEventProbability[pindex] = 0.; 

	// normalize event probability 
	double sum = 0; 
	bool flagall = true;
	for (int i = 0; i < fTreeVarNPermutations; ++i)
		{
			//			if (fTreeVarIntegral[i] >= 0.0)
			if (fTreeVarLogLikelihood[i] < 1e99)
				sum += fTreeVarEventProbability[i];
			else
				flagall = false; 
		}

	if (flagall)
		{
			for (int i = 0; i < fTreeVarNPermutations; ++i)
				{
					fTreeVarEventProbability[i] = fTreeVarEventProbability[i] / sum; 
				}
		}
	
	// sort for best permutation
	for (int i = 0; i < fTreeVarNPermutations; ++i)
		{
			int counter = 0; 
			for (int j = 0; j < fTreeVarNPermutations; ++j)
				//				if (fTreeVarLogLikelihood[i] < fTreeVarLogLikelihood[j])
				if (fTreeVarEventProbability[i] < fTreeVarEventProbability[j])
					counter++;
			fTreeVarBestPermutation[counter] = i; 
		}

	// loop over all parameters 
	for (int i = 0; i < fFitter -> Likelihood() -> NParameters(); ++i)
		{
			double * par    = fTreeVarParameters -> at(i); 
			double * parerr = fTreeVarParameterErrors -> at(i); 
			par[pindex]     = fFitter -> Likelihood() -> GetBestFitParameter(i); 
			parerr[pindex]  = fFitter -> Likelihood() -> GetBestFitParameterError(i); 
		}

	// loop over all particle type 
	for (KLFitter::Particles::ParticleType itype = KLFitter::Particles::kParton; itype <= KLFitter::Particles::kPhoton; ++itype)
		{
			// get particle container
			std::vector <TLorentzVector *> * momcontainer = (*fParticlesModel) -> ParticleContainer(itype); 

			// get number of particles in container 
			int n = int(momcontainer -> size()); 			
			
			// loop over particles 
			for (int i = 0; i < n; ++i)
				{
					// get variables 
					double * E  = fTreeVarModel -> at(counter); 
					double * px = fTreeVarModel -> at(++counter); 
					double * py = fTreeVarModel -> at(++counter); 
					double * pz = fTreeVarModel -> at(++counter); 
					double * m = fTreeVarModel -> at(++counter); 
					double * pt = fTreeVarModel -> at(++counter); 
					double * eta = fTreeVarModel -> at(++counter); 
					double * phi = fTreeVarModel -> at(++counter); 
					double * btag = fTreeVarModel -> at(++counter); 
					double * index = fTreeVarModel -> at(++counter); 
					
					// get four vector 
					TLorentzVector * lv = momcontainer -> at(i); 
					
					// fill variables
					E[pindex] = lv -> E(); 
					px[pindex] = lv -> Px(); 
					py[pindex] = lv -> Py(); 
					pz[pindex] = lv -> Pz(); 
					m[pindex] = lv -> M(); 
					pt[pindex] = lv -> Pt(); 
					eta[pindex] = lv -> Eta(); 
					phi[pindex] = lv -> Phi(); 
					if (itype == KLFitter::Particles::kParton &&
							(*fParticlesModel)->JetIndex(i)>=0) {
						btag[pindex] = (*fParticlesModel) -> FlavorTag(i); 
						index[pindex] = (fFitter->Permutations()->PermutationTable()) -> at(pindex) -> at(i); //(*fParticlesModel) -> JetIndex(i);
					}
					if (itype == KLFitter::Particles::kElectron) {
						index[pindex] = (*fParticlesModel) -> ElectronIndex(i);
					}
					if (itype == KLFitter::Particles::kMuon) {
						index[pindex] = (*fParticlesModel) -> MuonIndex(i);
					}
					if (itype == KLFitter::Particles::kPhoton) {
						index[pindex] = (*fParticlesModel) -> PhotonIndex(i);
					}


					// increase counter
					counter++; 
				}
		}

	// no error 
	return 1; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput::FillTreeMeasured()
{
	// check tree
	if (!fTreeMeasured)
		this -> CreateTreeMeasured(); 

	// initialize counter 
	int counter = 0; 

	// fill number of reconstructed objects
	fTreeVarNPartonsMeasured = (*fParticlesMeasured) -> NPartons(); 
	fTreeVarNElectronsMeasured = (*fParticlesMeasured) -> NElectrons(); 
	fTreeVarNMuonsMeasured = (*fParticlesMeasured) -> NMuons(); 
	fTreeVarNPhotonsMeasured = (*fParticlesMeasured) -> NPhotons(); 

	// loop over all particle type 
	for (KLFitter::Particles::ParticleType itype = KLFitter::Particles::kParton; itype <= KLFitter::Particles::kPhoton; ++itype)
		{
			// check if the branch should exist
			std::string name = TreeMeasuredBranchName(itype);
			if (name == "")
				continue;

			// get particle container
			std::vector <TLorentzVector *> * momcontainer = (*fParticlesMeasured) -> ParticleContainer(itype); 

			// get number of particles in container 
			int n = int(momcontainer -> size()); 			

			// get variables 
			std::vector<double> * vec_E   = fTreeVarMeasured -> at(counter); 
			std::vector<double> * vec_px  = fTreeVarMeasured -> at(++counter); 
			std::vector<double> * vec_py  = fTreeVarMeasured -> at(++counter); 
			std::vector<double> * vec_pz  = fTreeVarMeasured -> at(++counter); 
			std::vector<double> * vec_m   = fTreeVarMeasured -> at(++counter); 
			std::vector<double> * vec_pt  = fTreeVarMeasured -> at(++counter); 
			std::vector<double> * vec_eta = fTreeVarMeasured -> at(++counter); 
			std::vector<double> * vec_phi = fTreeVarMeasured -> at(++counter); 

			// clear the particle vectors
			vec_E   -> clear();
			vec_px  -> clear();
			vec_py  -> clear();
			vec_pz  -> clear();
			vec_m   -> clear();
			vec_pt  -> clear();
			vec_eta -> clear();
			vec_phi -> clear();

			// increase counter
			counter++; 

			// loop over particles 
			for (int i = 0; i < n; ++i)
				{
					// get four vector 
					TLorentzVector * lv = momcontainer -> at(i); 

					// fill variables
					vec_E   -> push_back( lv -> E()   ); 
					vec_px  -> push_back( lv -> Px()  ); 
					vec_py  -> push_back( lv -> Py()  ); 
					vec_pz  -> push_back( lv -> Pz()  ); 
					vec_m   -> push_back( lv -> M()   ); 
					vec_pt  -> push_back( lv -> Pt()  ); 
					vec_eta -> push_back( lv -> Eta() ); 
					vec_phi -> push_back( lv -> Phi() ); 

				}
		}

	// no error 
	return 1; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput::FillTreeSelected()
{
	// check tree
	if (!fTreeSelected)
		this -> CreateTreeSelected(); 

	// initialize counter 
	int counter = 0; 

	// fill number of reconstructed objects
	fTreeVarNPartonsSelected = (*fParticlesSelected) -> NPartons(); 
	fTreeVarNElectronsSelected = (*fParticlesSelected) -> NElectrons(); 
	fTreeVarNMuonsSelected = (*fParticlesSelected) -> NMuons(); 
	fTreeVarNPhotonsSelected = (*fParticlesSelected) -> NPhotons(); 

	// loop over all particle type 
	for (KLFitter::Particles::ParticleType itype = KLFitter::Particles::kParton; itype <= KLFitter::Particles::kPhoton; ++itype)
		{
			// check if the branch should exist
			std::string name = TreeMeasuredBranchName(itype);
			if (name == "")
				continue;

			// get particle container
			std::vector <TLorentzVector *> * momcontainer = (*fParticlesSelected) -> ParticleContainer(itype); 

			// get number of particles in container 
			int n = int(momcontainer -> size()); 			

			// get variables 
			std::vector<double> * vec_E    = fTreeVarSelected -> at(counter); 
			std::vector<double> * vec_px   = fTreeVarSelected -> at(++counter); 
			std::vector<double> * vec_py   = fTreeVarSelected -> at(++counter); 
			std::vector<double> * vec_pz   = fTreeVarSelected -> at(++counter); 
			std::vector<double> * vec_m    = fTreeVarSelected -> at(++counter); 
			std::vector<double> * vec_pt   = fTreeVarSelected -> at(++counter); 
			std::vector<double> * vec_eta  = fTreeVarSelected -> at(++counter); 
			std::vector<double> * vec_phi  = fTreeVarSelected -> at(++counter); 
			std::vector<double> * vec_btag = fTreeVarSelected -> at(++counter); 

			// clear the particle vectors
			vec_E    -> clear();
			vec_px   -> clear();
			vec_py   -> clear();
			vec_pz   -> clear();
			vec_m    -> clear();
			vec_pt   -> clear();
			vec_eta  -> clear();
			vec_phi  -> clear();
			vec_btag -> clear();

			// increase counter
			counter++; 

			// loop over particles 
			for (int i = 0; i < n; ++i)
				{
					// get four vector 
					TLorentzVector * lv = momcontainer -> at(i); 

					// fill variables
					vec_E   -> push_back( lv -> E()   ); 
					vec_px  -> push_back( lv -> Px()  ); 
					vec_py  -> push_back( lv -> Py()  ); 
					vec_pz  -> push_back( lv -> Pz()  ); 
					vec_m   -> push_back( lv -> M()   ); 
					vec_pt  -> push_back( lv -> Pt()  ); 
					vec_eta -> push_back( lv -> Eta() ); 
					vec_phi -> push_back( lv -> Phi() ); 
					vec_btag -> push_back( (*fParticlesSelected) -> FlavorTag(i)  ); 
				}
		}

	// no error 
	return 1; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput::FillTreeTruth()
{
	// check tree 
	if (!fTreeTruth)
		{
			// error code 
			int err = 1; 

			// create tree 
			if (fParticlesTruth)
				err = this -> CreateTreeTruth(); 
			
			else 
				return 0; 
		}

	// initialize counter 
	int counter = 0; 

	// loop over all particle type 
	for (KLFitter::Particles::ParticleType itype = KLFitter::Particles::kParton; itype <= KLFitter::Particles::kPhoton; ++itype)
		{
			// get particle container
			std::vector <TLorentzVector *> * momcontainer = (*fParticlesTruth) -> ParticleContainer(itype); 

			// get number of particles in container 
			int n = int(momcontainer -> size()); 			
			
			// loop over particles 
			for (int i = 0; i < n; ++i)
				{
					// get variables 
					double * E  = fTreeVarTruth -> at(counter); 
					double * px = fTreeVarTruth -> at(++counter); 
					double * py = fTreeVarTruth -> at(++counter); 
					double * pz = fTreeVarTruth -> at(++counter); 
					double * m = fTreeVarTruth -> at(++counter); 
					double * pt = fTreeVarTruth -> at(++counter); 
					double * eta = fTreeVarTruth -> at(++counter); 
					double * phi = fTreeVarTruth -> at(++counter); 

					// get four vector 
					TLorentzVector * lv = momcontainer -> at(i); 
					
					// fill variables
					*E = lv -> E(); 
					*px = lv -> Px(); 
					*py = lv -> Py(); 
					*pz = lv -> Pz(); 
					*m = lv -> M(); 
					*pt = lv -> Pt(); 
					*eta = lv -> Eta(); 
					*phi = lv -> Phi(); 
					// increase counter
					counter++; 
				}
		}

	// no error 
	return 1; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput::FillTreeMatching()
{
	// check tree 
	if (!fTreeMatching)
		{
			// error code 
			int err = 1; 
			
			// create tree 
			if (fMatchingTool)
				err = this -> CreateTreeMatching(); 
			
			else 
				return 0; 
		}
	
	// fill number of reconstructed objects
	fTreeVarNPartonsSelected = (*fParticlesSelected) -> NPartons(); 
	fTreeVarNElectronsSelected = (*fParticlesSelected) -> NElectrons(); 
	fTreeVarNMuonsSelected = (*fParticlesSelected) -> NMuons(); 
	fTreeVarNPhotonsSelected = (*fParticlesSelected) -> NPhotons(); 

	// loop over partons
	for (int i = 0; i < (*fParticlesTruth) -> NPartons(); ++i)
		{
			// get number of matched partons
			*fTreeVarNMatchedPartons -> at(i) = fMatchingTool -> NMatchedTruth(i, KLFitter::Particles::kParton); 

			// loop over reconstructed partons 
			for (int j = 0; j < fTreeVarNPartonsSelected; ++j)
				{
					(fTreeVarMatchedPartons -> at(i))[j] = (fMatchingTool -> ListMatchedTruth(i, KLFitter::Particles::kParton)).at(j);  
				}
		}

	// loop over electrons
	for (int i = 0; i < (*fParticlesTruth) -> NElectrons(); ++i)
		{
			// get number of matched electrons
			*fTreeVarNMatchedElectrons -> at(i) = fMatchingTool -> NMatchedTruth(i, KLFitter::Particles::kElectron); 

			// loop over reconstructed electrons 
			for (int j = 0; j < fTreeVarNElectronsSelected; ++j)
				{
					(fTreeVarMatchedElectrons -> at(i))[j] = (fMatchingTool -> ListMatchedTruth(i, KLFitter::Particles::kElectron)).at(j);  
				}
		}

	// loop over muons
	for (int i = 0; i < (*fParticlesTruth) -> NMuons(); ++i)
		{
			// get number of matched muons
			*fTreeVarNMatchedMuons -> at(i) = fMatchingTool -> NMatchedTruth(i, KLFitter::Particles::kMuon); 

			// loop over reconstructed muons 
			for (int j = 0; j < fTreeVarNMuonsSelected; ++j)
				{
					(fTreeVarMatchedMuons -> at(i))[j] = (fMatchingTool -> ListMatchedTruth(i, KLFitter::Particles::kMuon)).at(j);  
				}
		}

	// loop over photons
	for (int i = 0; i < (*fParticlesTruth) -> NPhotons(); ++i)
		{

			// get number of matched photons
			*fTreeVarNMatchedPhotons -> at(i) = fMatchingTool -> NMatchedTruth(i, KLFitter::Particles::kPhoton); 

			// loop over reconstructed photons 
			for (int j = 0; j < fTreeVarNPhotonsSelected; ++j)
				{
					(fTreeVarMatchedPhotons -> at(i))[j] = (fMatchingTool -> ListMatchedTruth(i, KLFitter::Particles::kPhoton)).at(j);  
				}
		}

	// no error 
	return 1; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput::FillTreeMap()
{
	// check tree 
	if (!fTreeMap)
		{
			// error code 
			int err = 1; 
			
			// create tree 
			if (fSelectionTool)
				err = this -> CreateTreeMap(); 
			
			else 
				return 0; 
		}
	
	// fill number of reconstructed objects
	fTreeVarNPartonsSelected = (*fParticlesSelected) -> NPartons(); 
	fTreeVarNElectronsSelected = (*fParticlesSelected) -> NElectrons(); 
	fTreeVarNMuonsSelected = (*fParticlesSelected) -> NMuons(); 
	fTreeVarNPhotonsSelected = (*fParticlesSelected) -> NPhotons(); 

	// get maps 
	for (int i = 0; i < fTreeVarNPartonsSelected; ++i)
		fTreeVarMapJets[i] = (fSelectionTool -> MapJets()).at(i); 
	for (int i = 0; i < fTreeVarNElectronsSelected; ++i)
		fTreeVarMapElectrons[i] = (fSelectionTool -> MapElectrons()).at(i); 
	for (int i = 0; i < fTreeVarNMuonsSelected; ++i)
		fTreeVarMapMuons[i] = (fSelectionTool -> MapMuons()).at(i); 
	for (int i = 0; i < fTreeVarNPhotonsSelected; ++i)
		fTreeVarMapPhotons[i] = (fSelectionTool -> MapPhotons()).at(i); 

	// no error 
	return 1; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput::FillTrees()
{
	// fill tree with truth particles 
	if (fParticlesTruth && fTreeTruth)
		fTreeTruth -> Fill(); 

	// fill tree with measured particles 
	if (fParticlesMeasured && fTreeMeasured)
		fTreeMeasured -> Fill(); 

	// fill tree with selected particles 
	if (fParticlesSelected && fTreeSelected)
		fTreeSelected -> Fill(); 

	// fill tree with model particles 
	if (fParticlesModel && fTreeModel)
		fTreeModel -> Fill(); 

	if (fMatchingTool)
		fTreeMatching -> Fill(); 

	if (fSelectionTool)
		fTreeMap -> Fill(); 

	fTreeVarEventNumber++; 

	// reset tree variables 
	for (int i = 0; i < fTreeVarNPermutations; ++i)
		{
			fTreeVarLogLikelihood[i] = 1e99; 
			fTreeVarIntegral[i] = -1; 
			fTreeVarEventProbability[i] = 0.; 
			fTreeVarMinuitStatus[i] = 0; 
		}

	// no error 
	return 1; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput::SetEventWeight(double weight)
{
	fEventWeight = weight;

	// no error
	return 0;
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput::SetPhotonType(bool isNotClassified, bool isRadTopProd, bool isHadTopRadDecay, bool isLepTopRadDecay, bool isHadWRadDecay, bool isLepWRadDecay)
{
	fIsNotClassified  = isNotClassified;
	fIsRadTopProd     = isRadTopProd;
	fIsHadTopRadDecay = isHadTopRadDecay;
	fIsLepTopRadDecay = isLepTopRadDecay;
	fIsHadWRadDecay   = isHadWRadDecay;
	fIsLepWRadDecay   = isLepWRadDecay;

	// no error
	return 0;
}

// --------------------------------------------------------- 
std::string KLFitter::InterfaceOutput::ModifyString(std::string str)
{
	int idx;

	while( (idx=str.find_first_of(' ')) >= 0 ) 
		str.replace(idx, 1, "_" );

	while( (idx=str.find_first_of('-')) >= 0 ) 
		str.replace(idx, 1, "_" );

	return str; 
}

// --------------------------------------------------------- 

std::string KLFitter::InterfaceOutput::TreeMeasuredBranchName(KLFitter::Particles::ParticleType pType)
{
	std::string name = "";
	switch(pType)
		{
		case KLFitter::Particles::kParton: name = "jet"; break;
		case KLFitter::Particles::kElectron: name = "electron"; break;
		case KLFitter::Particles::kMuon: name = "muon"; break;
		case KLFitter::Particles::kPhoton: name = "photon"; break;
		case KLFitter::Particles::kTau: name = ""; break; // there is no measured tau
		case KLFitter::Particles::kNeutrino: name = ""; break; // there is no measured neutrino
		case KLFitter::Particles::kBoson: name = ""; break; // there is no measured (W) boson
		default: name = ""; break;
		}

	return name;
}

// --------------------------------------------------------- 
