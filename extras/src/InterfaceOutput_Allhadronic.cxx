#include "InterfaceOutput_Allhadronic.h" 
#include "LikelihoodBase.h"
#include "Fitter.h" 
#include "MatchingTool.h" 
#include "SelectionTool.h" 
#include "Permutations.h" 

#include <TFile.h>
#include <TTree.h> 

#include <iostream> 
#include <map> 

// --------------------------------------------------------- 
KLFitter::InterfaceOutput_Allhadronic::InterfaceOutput_Allhadronic()
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
  fTreeVarLogLikelihoodComp_TF_bhad1 = 0; 
  fTreeVarLogLikelihoodComp_TF_bhad2 = 0; 
  fTreeVarLogLikelihoodComp_TF_lq1 = 0; 
  fTreeVarLogLikelihoodComp_TF_lq2 = 0; 
  fTreeVarLogLikelihoodComp_TF_lq3 = 0;
  fTreeVarLogLikelihoodComp_TF_lq4 = 0;
  fTreeVarLogLikelihoodComp_BW_Whad1 = 0;
  fTreeVarLogLikelihoodComp_BW_Whad2 = 0;
  fTreeVarLogLikelihoodComp_BW_Thad1 = 0;
  fTreeVarLogLikelihoodComp_BW_Thad2 = 0;
//  fTreeVarLogLikelihoodComponents = 0;
  fTreeVarIntegral = 0;
  fTreeVarEventProbability = 0; 
  fTreeVarMinuitStatus = 0; 
  fTreeVarConvergenceStatus = 0; 
  fTreeVarEventNumber = 0; 
  fTreeVarParameters = new std::vector<std::vector<double> *>(0); 
  fTreeVarParameterErrors = new std::vector<std::vector<double> *>(0); 
  fTreeVarModel = new std::vector<std::vector<double> *>(0);
  fTreeIntVarModel = new std::vector<std::vector<int> *>(0);  
  fTreeVarTruth = new std::vector<std::vector<double> *>(0); 
  fTreeVarMeasured = new std::vector<std::vector<double> *>(0); 
  fTreeVarSelected = new std::vector<std::vector<double> *>(0); 
  fTreeVarNMatchedPartons = new std::vector<int *>(0); 
  fTreeVarNMatchedElectrons = new std::vector<int *>(0); 
  fTreeVarNMatchedMuons = new std::vector<int *>(0); 
  fTreeVarNMatchedPhotons = new std::vector<int *>(0); 
  fTreeVarMatchedPartons = new std::vector <std::vector<int> *>(0); 
  fTreeVarMatchedElectrons = new std::vector <std::vector<int> *>(0); 
  fTreeVarMatchedMuons = new std::vector <std::vector<int> *>(0); 
  fTreeVarMatchedPhotons = new std::vector <std::vector<int> *>(0); 
  fTreeVarMapJets = 0; 
  fTreeVarMapElectrons =  0;
  fTreeVarMapMuons = 0; 
  fTreeVarMapPhotons = 0; 
  fEventWeight = 0.;
  fPileupWeight = 0.;  
  fIsNotClassified = false;
  fIsRadTopProd = false;
  fIsHadTopRadDecay = false;
  fIsLepTopRadDecay = false;
  fIsHadWRadDecay = false;
  fIsLepWRadDecay = false;
}

// --------------------------------------------------------- 
KLFitter::InterfaceOutput_Allhadronic::~InterfaceOutput_Allhadronic()
{
  if (fTreeVarBestPermutation)
    delete fTreeVarBestPermutation; 

  if (fTreeVarLogLikelihood)
    delete fTreeVarLogLikelihood; 

  if(fTreeVarLogLikelihoodComp_TF_bhad1) 
    delete fTreeVarLogLikelihoodComp_TF_bhad1;

  if(fTreeVarLogLikelihoodComp_TF_bhad2) 
    delete fTreeVarLogLikelihoodComp_TF_bhad2;

  if(fTreeVarLogLikelihoodComp_TF_lq1) 
    delete fTreeVarLogLikelihoodComp_TF_lq1;

  if(fTreeVarLogLikelihoodComp_TF_lq2) 
    delete fTreeVarLogLikelihoodComp_TF_lq2;

  if(fTreeVarLogLikelihoodComp_TF_lq3)
    delete fTreeVarLogLikelihoodComp_TF_lq3;

  if(fTreeVarLogLikelihoodComp_TF_lq4)
    delete fTreeVarLogLikelihoodComp_TF_lq4;

  if(fTreeVarLogLikelihoodComp_BW_Whad1)
    delete fTreeVarLogLikelihoodComp_BW_Whad1;

  if(fTreeVarLogLikelihoodComp_BW_Whad2)
    delete fTreeVarLogLikelihoodComp_BW_Whad2;

  if(fTreeVarLogLikelihoodComp_BW_Thad1)
    delete fTreeVarLogLikelihoodComp_BW_Thad1;

  if(fTreeVarLogLikelihoodComp_BW_Thad2)
    delete fTreeVarLogLikelihoodComp_BW_Thad2;

//  if (fTreeVarLogLikelihoodComponents)
//    delete fTreeVarLogLikelihoodComponents; 

  if (fTreeVarIntegral)
    delete fTreeVarIntegral; 

  if (fTreeVarEventProbability)
    delete fTreeVarEventProbability; 

  if (fTreeVarMinuitStatus)
    delete fTreeVarMinuitStatus; 

  if (fTreeVarConvergenceStatus)
    delete fTreeVarConvergenceStatus; 

  while (!fTreeVarParameters->empty())
    {
      std::vector<double>* d = fTreeVarParameters->front(); 
      fTreeVarParameters->erase(fTreeVarParameters->begin()); 
      delete d; 
    }
  delete fTreeVarParameters; 

  while (!fTreeVarParameterErrors->empty())
    {
      std::vector<double>* d = fTreeVarParameterErrors->front(); 
      fTreeVarParameterErrors->erase(fTreeVarParameterErrors->begin()); 
      delete d; 
    }
  delete fTreeVarParameterErrors; 
        
  while (!fTreeVarModel->empty())
    {
      std::vector<double>* d = fTreeVarModel->front(); 
      fTreeVarModel->erase(fTreeVarModel->begin()); 
      delete d; 
    }
  delete fTreeVarModel;
  
  while (!fTreeIntVarModel->empty())
    {
      std::vector<int>* d = fTreeIntVarModel->front(); 
      fTreeIntVarModel->erase(fTreeIntVarModel->begin()); 
      delete d; 
    }
  delete fTreeIntVarModel;  
        
  while (!fTreeVarTruth->empty())
    {
      std::vector<double>* d = fTreeVarTruth->front(); 
      fTreeVarTruth->erase(fTreeVarTruth->begin()); 
      delete d; 
    }
  delete fTreeVarTruth; 
        
  while (!fTreeVarMeasured->empty())
    {
      std::vector<double>* d = fTreeVarMeasured-> front(); 
      fTreeVarMeasured->erase(fTreeVarMeasured->begin()); 
      delete d; 
    }
  delete fTreeVarMeasured;

  while (!fTreeVarSelected->empty())
    {
      std::vector<double>* d = fTreeVarSelected->front(); 
      fTreeVarSelected->erase(fTreeVarSelected->begin()); 
      delete d; 
    }
  delete fTreeVarSelected;
        
  while (!fTreeVarNMatchedPartons->empty())
    {
      int* i = fTreeVarNMatchedPartons->front(); 
      fTreeVarNMatchedPartons->erase(fTreeVarNMatchedPartons->begin()); 
      delete i; 
    }
  delete fTreeVarNMatchedPartons;

  while (!fTreeVarNMatchedElectrons->empty())
    {
      int* i = fTreeVarNMatchedElectrons->front(); 
      fTreeVarNMatchedElectrons->erase(fTreeVarNMatchedElectrons->begin()); 
      delete i; 
    }
  delete fTreeVarNMatchedElectrons;

  while (!fTreeVarNMatchedMuons->empty())
    {
      int* i = fTreeVarNMatchedMuons->front(); 
      fTreeVarNMatchedMuons->erase(fTreeVarNMatchedMuons->begin()); 
      delete i; 
    }
  delete fTreeVarNMatchedMuons;

  while (!fTreeVarNMatchedPhotons->empty())
    {
      int* i = fTreeVarNMatchedPhotons->front(); 
      fTreeVarNMatchedPhotons->erase(fTreeVarNMatchedPhotons->begin()); 
      delete i; 
    }
  delete fTreeVarNMatchedPhotons;

  while (!fTreeVarMatchedPartons->empty())
    {
      std::vector<int>* i = fTreeVarMatchedPartons->front(); 
      fTreeVarMatchedPartons->erase(fTreeVarMatchedPartons->begin()); 
      delete i; 
    }
  delete fTreeVarMatchedPartons;

  while (!fTreeVarMatchedElectrons->empty())
    {
      std::vector<int>* i = fTreeVarMatchedElectrons->front(); 
      fTreeVarMatchedElectrons->erase(fTreeVarMatchedElectrons->begin()); 
      delete i; 
    }
  delete fTreeVarMatchedElectrons;

  while (!fTreeVarMatchedMuons->empty())
    {
      std::vector<int>* i = fTreeVarMatchedMuons->front(); 
      fTreeVarMatchedMuons->erase(fTreeVarMatchedMuons->begin()); 
      delete i; 
    }
  delete fTreeVarMatchedMuons;

  while (!fTreeVarMatchedPhotons->empty())
    {
      std::vector<int>* i = fTreeVarMatchedPhotons->front(); 
      fTreeVarMatchedPhotons->erase(fTreeVarMatchedPhotons->begin()); 
      delete i; 
    }
  delete fTreeVarMatchedPhotons;

  if (fTreeVarMapJets)
    delete fTreeVarMapJets; 

  if (fTreeVarMapElectrons)
    delete fTreeVarMapElectrons; 

  if (fTreeVarMapMuons)
    delete fTreeVarMapMuons; 

  if (fTreeVarMapPhotons)
    delete fTreeVarMapPhotons; 

  if (fTreeTruth)
    delete fTreeTruth; 
        
  if (fTreeModel)
    delete fTreeModel; 

  if (fTreeMeasured)
    delete fTreeMeasured; 

  if (fTreeSelected)
    delete fTreeSelected; 
        
  if (fTreeMatching)
    delete fTreeMatching; 

  if (fTreeMap)
    delete fTreeMap; 

}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput_Allhadronic::SetFitter(KLFitter::Fitter* fitter)
{
  // check if fitter exists 
  if (!fitter)
    {
      std::cout << "KLFitter::InterfaceOutput::SetFitter(). Fitter does not exist." << std::endl; 
      return 0; 
    }

  // set pointers to pointer 
  fParticlesModel = fitter->Likelihood()->PParticlesModel(); 
  fParticlesSelected = fitter->PParticles(); 

  // set fitter
  fFitter = fitter;

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput_Allhadronic::SetMatchingTool(KLFitter::MatchingTool* matchingtool)
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
int KLFitter::InterfaceOutput_Allhadronic::SetSelectionTool(KLFitter::SelectionTool* selectiontool)
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
int KLFitter::InterfaceOutput_Allhadronic::OpenRootFile(const char* filename, Option_t* opt)
{
  // define error code 
  int err = 1; 

  // open file 
  err *= KLFitter::InterfaceRoot::OpenRootFile(filename, opt); 

  // return error code 
  return err; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput_Allhadronic::CloseRootFile()
{
  // define error code 
  int err = 1; 

  // check if file exists 
  if (!fRootFile)
    return 0; 
        
  // check if file is open 
  if (!fRootFile->IsOpen())
    return 0; 

  if (fTreeTruth)
    fRootFile->WriteTObject(fTreeTruth); 
        
  if (fTreeModel)
    fRootFile->WriteTObject(fTreeModel); 

  if (fTreeMeasured)
    fRootFile->WriteTObject(fTreeMeasured); 

  if (fTreeSelected)
    fRootFile->WriteTObject(fTreeSelected); 

  if (fTreeMatching)
    fRootFile->WriteTObject(fTreeMatching); 

  if (fTreeMap)
    fRootFile->WriteTObject(fTreeMap); 

  // close file 
  KLFitter::InterfaceRoot::CloseRootFile(); 

  // return error code 
  return err; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput_Allhadronic::CreateTreeModel()
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
  int nperm = fFitter->Permutations()->NPermutations(); 

  // create new tree 
  fTreeModel = new TTree("TreeModel", "TreeModel"); 

  // reset variables
  fTreeVarNPermutations = nperm; 
  fTreeVarNBTags=0;
  fTreeVarBestPermutation = new std::vector<int>(0);
  fTreeVarLogLikelihood = new std::vector<double>(0);
  fTreeVarLogLikelihoodComp_TF_bhad1 = new std::vector<double>(0); 
  fTreeVarLogLikelihoodComp_TF_bhad2 = new std::vector<double>(0); 
  fTreeVarLogLikelihoodComp_TF_lq1 = new std::vector<double>(0); 
  fTreeVarLogLikelihoodComp_TF_lq2 = new std::vector<double>(0); 
  fTreeVarLogLikelihoodComp_TF_lq3 = new std::vector<double>(0);
  fTreeVarLogLikelihoodComp_TF_lq4 = new std::vector<double>(0);
  fTreeVarLogLikelihoodComp_BW_Whad1 = new std::vector<double>(0);
  fTreeVarLogLikelihoodComp_BW_Whad2 = new std::vector<double>(0);
  fTreeVarLogLikelihoodComp_BW_Thad1 = new std::vector<double>(0);
  fTreeVarLogLikelihoodComp_BW_Thad2 = new std::vector<double>(0);
//  fTreeVarLogLikelihoodComponents = new std::vector<std::vector<double> >(0);
  fTreeVarIntegral = new std::vector<double>(0);
  fTreeVarEventProbability = new std::vector<double>(0);
  fTreeVarMinuitStatus = new std::vector<double>(0);
  fTreeVarConvergenceStatus = new std::vector<unsigned int>(0);

  // set branches for event variables
  fTreeModel->Branch("EventNumber", &fTreeVarEventNumber, "EventNumber/I"); 
  fTreeModel->Branch("N_permutations", &fTreeVarNPermutations, "N_permutations/I"); 
  fTreeModel->Branch("N_btags", &fTreeVarNBTags, "N_btags/I"); 
  fTreeModel->Branch("best_permutation", fTreeVarBestPermutation);
  fTreeModel->Branch("LogLikelihood", fTreeVarLogLikelihood);
  fTreeModel->Branch("LogLikelihoodComp_TF_bhad1",  fTreeVarLogLikelihoodComp_TF_bhad1);
  fTreeModel->Branch("LogLikelihoodComp_TF_bhad2",  fTreeVarLogLikelihoodComp_TF_bhad2);
  fTreeModel->Branch("LogLikelihoodComp_TF_lq1",  fTreeVarLogLikelihoodComp_TF_lq1);
  fTreeModel->Branch("LogLikelihoodComp_TF_lq2",  fTreeVarLogLikelihoodComp_TF_lq2);
  fTreeModel->Branch("LogLikelihoodComp_TF_lq3",  fTreeVarLogLikelihoodComp_TF_lq3);
  fTreeModel->Branch("LogLikelihoodComp_TF_lq4",  fTreeVarLogLikelihoodComp_TF_lq4);
  fTreeModel->Branch("LogLikelihoodComp_BW_Whad1",  fTreeVarLogLikelihoodComp_BW_Whad1);
  fTreeModel->Branch("LogLikelihoodComp_BW_Whad2",  fTreeVarLogLikelihoodComp_BW_Whad2);
  fTreeModel->Branch("LogLikelihoodComp_BW_Thad1",  fTreeVarLogLikelihoodComp_BW_Thad1);
  fTreeModel->Branch("LogLikelihoodComp_BW_Thad2",  fTreeVarLogLikelihoodComp_BW_Thad2);
//  fTreeModel->Branch("LogLikelihoodComponents", fTreeVarLogLikelihoodComponents);
  fTreeModel->Branch("Integral", fTreeVarIntegral);
  fTreeModel->Branch("EventProbability", fTreeVarEventProbability);
  fTreeModel->Branch("MinuitStatus", fTreeVarMinuitStatus);
  fTreeModel->Branch("ConvergenceStatusBit", fTreeVarConvergenceStatus);

  // loop over all parameters 
  for (int i = 0; i < fFitter->Likelihood()->NParameters(); ++i)
    {
      std::vector<double>* par = new std::vector<double>(0); 
      std::vector<double>* parerr = new std::vector<double>(0); 
      fTreeVarParameters->push_back(par); 
      fTreeVarParameterErrors->push_back(parerr); 

      fTreeModel->Branch( this->ModifyString( "par_" + fFitter->Likelihood()->GetParameter(i)->GetName() ).data(), par);
      fTreeModel->Branch( this->ModifyString( "parerr_" + fFitter->Likelihood()->GetParameter(i)->GetName() ).data(), parerr);
    }


  // loop over all particle type 
  for (KLFitter::Particles::ParticleType itype = KLFitter::Particles::kParton; itype <= KLFitter::Particles::kPhoton; ++itype)
    {
      // get particle container
      std::vector <TLorentzVector *>* momcontainer = (*fParticlesModel)->ParticleContainer(itype); 
      std::vector <std::string>* namecontainer = (*fParticlesModel)->ParticleNameContainer(itype); 

      // get number of particles in container 
      int n = int(momcontainer->size());                        
                        
      // loop over particles 
      for (int i = 0; i < n; ++i)
        {
          // get name 
          std::string name(namecontainer->at(i)); 

          // create new pointer to variables 
          std::vector<double>* E = new std::vector<double>(0); 
          std::vector<double>* px = new std::vector<double>(0); 
          std::vector<double>* py = new std::vector<double>(0); 
          std::vector<double>* pz = new std::vector<double>(0); 
          std::vector<double>* m = new std::vector<double>(0); 
          std::vector<double>* pt = new std::vector<double>(0); 
          std::vector<double>* eta = new std::vector<double>(0); 
          std::vector<double>* phi = new std::vector<double>(0); 
          std::vector<double>* btag = new std::vector<double>(0); 
          std::vector<int>* index = new std::vector<int>(0);

          // add variables to vector 
          fTreeVarModel->push_back(E); 
          fTreeVarModel->push_back(px); 
          fTreeVarModel->push_back(py); 
          fTreeVarModel->push_back(pz); 
          fTreeVarModel->push_back(m); 
          fTreeVarModel->push_back(pt); 
          fTreeVarModel->push_back(eta); 
          fTreeVarModel->push_back(phi); 
          fTreeVarModel->push_back(btag); 
          fTreeIntVarModel->push_back(index); 

          // create new branches                                        
          fTreeModel->Branch(this->ModifyString((name+"_E")).data(), E);
          fTreeModel->Branch(this->ModifyString((name+"_px")).data(), px);
          fTreeModel->Branch(this->ModifyString((name+"_py")).data(), py);
          fTreeModel->Branch(this->ModifyString((name+"_pz")).data(), pz);
          fTreeModel->Branch(this->ModifyString((name+"_m")).data(), m);
          fTreeModel->Branch(this->ModifyString((name+"_pt")).data(), pt);
          fTreeModel->Branch(this->ModifyString((name+"_eta")).data(), eta);
          fTreeModel->Branch(this->ModifyString((name+"_phi")).data(), phi);
          if (itype == KLFitter::Particles::kParton && (*fParticlesModel)->JetIndex(i)>=0) {
            fTreeModel->Branch(this->ModifyString((name+"_btag")).data(), btag);
            fTreeModel->Branch(this->ModifyString((name+"_index")).data(), index);
          }                                     
          if (itype == KLFitter::Particles::kElectron && (*fParticlesModel)->ElectronIndex(i)>=0) 
            fTreeModel->Branch(this->ModifyString((name+"_index")).data(), index);
          if (itype == KLFitter::Particles::kMuon && (*fParticlesModel)->MuonIndex(i)>=0) 
            fTreeModel->Branch(this->ModifyString((name+"_index")).data(), index);
          if (itype == KLFitter::Particles::kPhoton && (*fParticlesModel)->PhotonIndex(i)>=0) 
            fTreeModel->Branch(this->ModifyString((name+"_index")).data(), index);
        }
    }

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput_Allhadronic::CreateTreeMeasured()
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
  fTreeMeasured->Branch("Weight", &fEventWeight, "Weight/D"); 
  fTreeMeasured->Branch("PileupWeight", &fPileupWeight, "PileupWeight/D");          
  fTreeMeasured->Branch("N_jets", &fTreeVarNPartonsMeasured, "N_jets/I"); 
  fTreeMeasured->Branch("N_electrons", &fTreeVarNElectronsMeasured, "N_electrons/I"); 
  fTreeMeasured->Branch("N_muons", &fTreeVarNMuonsMeasured, "N_muons/I"); 
  fTreeMeasured->Branch("N_photons", &fTreeVarNPhotonsMeasured, "N_photons/I"); 

  // loop over all particle type 
  for (KLFitter::Particles::ParticleType itype = KLFitter::Particles::kParton; itype <= KLFitter::Particles::kPhoton; ++itype)
    {
      // get the name of the branch
      std::string name = TreeMeasuredBranchName(itype);
      if (name == "")
        continue;

      // create new pointer to variables 
      std::vector<double>* vec_E   = new std::vector<double>(0); 
      std::vector<double>* vec_px  = new std::vector<double>(0); 
      std::vector<double>* vec_py  = new std::vector<double>(0); 
      std::vector<double>* vec_pz  = new std::vector<double>(0); 
      std::vector<double>* vec_m   = new std::vector<double>(0); 
      std::vector<double>* vec_pt  = new std::vector<double>(0); 
      std::vector<double>* vec_eta = new std::vector<double>(0); 
      std::vector<double>* vec_phi = new std::vector<double>(0); 

      // add variables to vector 
      fTreeVarMeasured->push_back(vec_E); 
      fTreeVarMeasured->push_back(vec_px); 
      fTreeVarMeasured->push_back(vec_py); 
      fTreeVarMeasured->push_back(vec_pz); 
      fTreeVarMeasured->push_back(vec_m); 
      fTreeVarMeasured->push_back(vec_pt); 
      fTreeVarMeasured->push_back(vec_eta); 
      fTreeVarMeasured->push_back(vec_phi); 

      // create new branches                                    
      fTreeMeasured->Branch((name+"_E").data(),   vec_E  ); 
      fTreeMeasured->Branch((name+"_px").data(),  vec_px );
      fTreeMeasured->Branch((name+"_py").data(),  vec_py );
      fTreeMeasured->Branch((name+"_pz").data(),  vec_pz );
      fTreeMeasured->Branch((name+"_m").data(),   vec_m  ); 
      fTreeMeasured->Branch((name+"_pt").data(),  vec_pt );
      fTreeMeasured->Branch((name+"_eta").data(), vec_eta);
      fTreeMeasured->Branch((name+"_phi").data(), vec_phi);
    }

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput_Allhadronic::CreateTreeSelected()
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

  fTreeSelected->Branch("N_jets", &fTreeVarNPartonsSelected, "N_jets/I"); 
  fTreeSelected->Branch("N_electrons", &fTreeVarNElectronsSelected, "N_electrons/I"); 
  fTreeSelected->Branch("N_muons", &fTreeVarNMuonsSelected, "N_muons/I"); 
  fTreeSelected->Branch("N_photons", &fTreeVarNPhotonsSelected, "N_photons/I"); 

  // loop over all particle type 
  for (KLFitter::Particles::ParticleType itype = KLFitter::Particles::kParton; itype <= KLFitter::Particles::kPhoton; ++itype)
    {
      // get the name of the branch
      std::string name = TreeMeasuredBranchName(itype);
      if (name == "")
        continue;

      // create new pointer to variables 
      std::vector<double>* vec_E    = new std::vector<double>(0); 
      std::vector<double>* vec_px   = new std::vector<double>(0); 
      std::vector<double>* vec_py   = new std::vector<double>(0); 
      std::vector<double>* vec_pz   = new std::vector<double>(0); 
      std::vector<double>* vec_m    = new std::vector<double>(0); 
      std::vector<double>* vec_pt   = new std::vector<double>(0); 
      std::vector<double>* vec_eta  = new std::vector<double>(0); 
      std::vector<double>* vec_phi  = new std::vector<double>(0); 
      std::vector<double>* vec_btag = new std::vector<double>(0); 

      // add variables to vector 
      fTreeVarSelected->push_back(vec_E); 
      fTreeVarSelected->push_back(vec_px); 
      fTreeVarSelected->push_back(vec_py); 
      fTreeVarSelected->push_back(vec_pz); 
      fTreeVarSelected->push_back(vec_m); 
      fTreeVarSelected->push_back(vec_pt); 
      fTreeVarSelected->push_back(vec_eta); 
      fTreeVarSelected->push_back(vec_phi); 
      fTreeVarSelected->push_back(vec_btag); 

      // create new branches                                    
      fTreeSelected->Branch((name+"_E").data(),   vec_E  ); 
      fTreeSelected->Branch((name+"_px").data(),  vec_px );
      fTreeSelected->Branch((name+"_py").data(),  vec_py );
      fTreeSelected->Branch((name+"_pz").data(),  vec_pz );
      fTreeSelected->Branch((name+"_m").data(),   vec_m  ); 
      fTreeSelected->Branch((name+"_pt").data(),  vec_pt );
      fTreeSelected->Branch((name+"_eta").data(), vec_eta);
      fTreeSelected->Branch((name+"_phi").data(), vec_phi);
      if (itype == KLFitter::Particles::kParton)
        fTreeSelected->Branch((name+"_btag").data(), vec_btag);
    }

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput_Allhadronic::CreateTreeTruth()
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
      std::vector <TLorentzVector *>* momcontainer = (*fParticlesTruth)->ParticleContainer(itype); 
      std::vector <std::string>* namecontainer = (*fParticlesTruth)->ParticleNameContainer(itype); 

      // get number of particles in container 
      int n = int(momcontainer->size());                        

      // loop over particles 
      for (int i = 0; i < n; ++i)
        {
          // get name 
          std::string name(namecontainer->at(i)); 

          // create new pointer to variables 
          std::vector<double>* E = new std::vector<double>(0); 
          std::vector<double>* px = new std::vector<double>(0); 
          std::vector<double>* py = new std::vector<double>(0); 
          std::vector<double>* pz = new std::vector<double>(0); 
          std::vector<double>* m = new std::vector<double>(0); 
          std::vector<double>* pt = new std::vector<double>(0); 
          std::vector<double>* eta = new std::vector<double>(0); 
          std::vector<double>* phi = new std::vector<double>(0); 

          // add variables to vector 
          fTreeVarTruth->push_back(E); 
          fTreeVarTruth->push_back(px); 
          fTreeVarTruth->push_back(py); 
          fTreeVarTruth->push_back(pz); 
          fTreeVarTruth->push_back(m); 
          fTreeVarTruth->push_back(pt); 
          fTreeVarTruth->push_back(eta); 
          fTreeVarTruth->push_back(phi); 

          // create new branches                                        
          fTreeTruth->Branch(this->ModifyString((name+"_E")).data(), E);
          fTreeTruth->Branch(this->ModifyString((name+"_px")).data(), px);
          fTreeTruth->Branch(this->ModifyString((name+"_py")).data(), py);
          fTreeTruth->Branch(this->ModifyString((name+"_pz")).data(), pz);
          fTreeTruth->Branch(this->ModifyString((name+"_m")).data(), m);
          fTreeTruth->Branch(this->ModifyString((name+"_pt")).data(), pt);
          fTreeTruth->Branch(this->ModifyString((name+"_eta")).data(), eta);
          fTreeTruth->Branch(this->ModifyString((name+"_phi")).data(), phi);
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
int KLFitter::InterfaceOutput_Allhadronic::CreateTreeMatching()
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

  fTreeMatching->Branch("N_jets", &fTreeVarNPartonsSelected, "N_jets/I"); 
  fTreeMatching->Branch("N_electrons", &fTreeVarNElectronsSelected, "N_electrons/I"); 
  fTreeMatching->Branch("N_muons", &fTreeVarNMuonsSelected, "N_muons/I"); 
  fTreeMatching->Branch("N_photons", &fTreeVarNPhotonsSelected, "N_photons/I"); 

  // get parton container
  std::vector <std::string>* namecontainer = (*fParticlesTruth)->ParticleNameContainer(KLFitter::Particles::kParton); 

  // get number of particles in container 
  int ntruth = int(namecontainer->size());                      

  // loop over particles 
  for (int i = 0; i < ntruth; ++i)
    {
      // get name 
      std::string name(namecontainer->at(i)); 

      int* pi = new int(-1); 
      fTreeVarNMatchedPartons->push_back(pi); 

      // if fParticlesSelected are not yet defined, set the maximum value to 10
      int nreco_partons_max = fParticlesSelected ? (*fParticlesSelected)->NPartons() : 10;
      std::vector<int>* indices = new std::vector<int>(nreco_partons_max, -1); 
      fTreeVarMatchedPartons->push_back(indices); 

      // create new branches 
      fTreeMatching->Branch( this->ModifyString("Nmatches_"+name).data(), pi);
      fTreeMatching->Branch( this->ModifyString("Matches_"+name).data(), indices);
    }
        
  // get electron container
  namecontainer = (*fParticlesTruth)->ParticleNameContainer(KLFitter::Particles::kElectron); 

  // get number of particles in container 
  int n = int(namecontainer->size());                   

  // loop over particles 
  for (int i = 0; i < n; ++i)
    {
      // get name 
      std::string name(namecontainer->at(i)); 
                        
      int* pi = new int(-1); 
      fTreeVarNMatchedElectrons->push_back(pi); 

      // if fParticlesSelected are not yet defined, set the maximum value to 10
      int nreco_electrons_max = fParticlesSelected ? (*fParticlesSelected)->NElectrons() : 10;
      std::vector<int>* indices = new std::vector<int>(nreco_electrons_max, -1); 
      fTreeVarMatchedElectrons->push_back(indices); 

      // create new branches 
      fTreeMatching->Branch( this->ModifyString("Nmatches_"+name).data(), pi);
      fTreeMatching->Branch( this->ModifyString("Matches_"+name).data(), indices);
    }
        
  // get muon container
  namecontainer = (*fParticlesTruth)->ParticleNameContainer(KLFitter::Particles::kMuon); 

  // get number of particles in container 
  n = int(namecontainer->size());                       

  // loop over particles 
  for (int i = 0; i < n; ++i)
    {
      // get name 
      std::string name(namecontainer->at(i)); 
                        
      int* pi = new int(-1); 
      fTreeVarNMatchedMuons->push_back(pi); 

      // if fParticlesSelected are not yet defined, set the maximum value to 10
      int nreco_muons_max = fParticlesSelected ? (*fParticlesSelected)->NMuons() : 10;
      std::vector<int>* indices = new std::vector<int>(nreco_muons_max, -1); 
      fTreeVarMatchedMuons->push_back(indices); 

      // create new branches 
      fTreeMatching->Branch( this->ModifyString("Nmatches_"+name).data(), pi);
      fTreeMatching->Branch( this->ModifyString("Matches_"+name).data(), indices);
    }   

  // get photon container
  namecontainer = (*fParticlesTruth)->ParticleNameContainer(KLFitter::Particles::kPhoton); 

  // get number of particles in container 
  n = int(namecontainer->size());                       

  // loop over particles 
  for (int i = 0; i < n; ++i)
    {
      // get name 
      std::string name(namecontainer->at(i)); 
                        
      int* pi = new int(-1); 
      fTreeVarNMatchedPhotons->push_back(pi); 

      // if fParticlesSelected are not yet defined, set the maximum value to 10
      int nreco_photons_max = fParticlesSelected ? (*fParticlesSelected)->NPhotons() : 10;
      std::vector<int>* indices = new std::vector<int>(nreco_photons_max, -1); 
      fTreeVarMatchedPhotons->push_back(indices); 

      // create new branches 
      fTreeMatching->Branch( this->ModifyString("Nmatches_"+name).data(), pi);
      fTreeMatching->Branch( this->ModifyString("Matches_"+name).data(), indices);
    }   

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput_Allhadronic::CreateTreeMap()
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
  fTreeVarMapJets = new std::vector<int>(0); 
  fTreeVarMapElectrons = new std::vector<int>(0); 
  fTreeVarMapMuons = new std::vector<int>(0);
  fTreeVarMapPhotons = new std::vector<int>(0); 

  // create new tree 
  fTreeMap = new TTree("TreeMap", "TreeMap"); 

  fTreeMap->Branch("N_jets", &fTreeVarNPartonsSelected, "N_jets/I"); 
  fTreeMap->Branch("N_electrons", &fTreeVarNElectronsSelected, "N_electrons/I"); 
  fTreeMap->Branch("N_muons", &fTreeVarNMuonsSelected, "N_muons/I"); 
  fTreeMap->Branch("N_photons", &fTreeVarNPhotonsSelected, "N_photons/I"); 

  // set branches for event variables
  fTreeMap->Branch("EventNumber", &fTreeVarEventNumber, "EventNumber/I"); 
  fTreeMap->Branch("Index_jet", fTreeVarMapJets);
  fTreeMap->Branch("Index_electron", fTreeVarMapElectrons);
  fTreeMap->Branch("Index_muon", fTreeVarMapMuons);
  fTreeMap->Branch("Index_photon", fTreeVarMapPhotons);

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput_Allhadronic::CreateTrees()
{
  // error code 
  int err = 1; 
        
  // create tree for truth particles 
  if (fParticlesTruth)
    err *= this->CreateTreeTruth(); 

  // create tree for measured particles 
  err *= this->CreateTreeMeasured(); 

  // create tree for selected particles 
  err *= this->CreateTreeSelected(); 

  // create tree for model particles 
  err *= this->CreateTreeModel(); 

  // create tree for matching informations
  if (fMatchingTool && fParticlesTruth)
    err *= this->CreateTreeMatching(); 

  if (fSelectionTool)
    err *= this->CreateTreeMap(); 

  // return error code 
  return err; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput_Allhadronic::FillTreeModelPermutation()
{
  // check tree 
  if (!fTreeModel) 
    this->CreateTreeModel(); 

  // make sure that the model particles are being built in the likelihood
  fParticlesModel = fFitter->Likelihood()->PParticlesModel(); 

  // initialize counter 
  int counter = 0; 

  // create new permutation table for writing out the index
  //    int npartons = (*fParticlesSelected)->NPartons(); 
  //    std::vector < std::vector<int>* >* table_partons = new std::vector < std::vector<int>* >(0); 
  //    fFitter->Permutations()-> CreateSubTable(npartons, table_partons); 

  // get permutation index
  int pindex =  fFitter->Permutations()->PermutationIndex(); 
        
  // fill event variables
  fTreeVarNBTags = (*fParticlesModel)->NBTags(); 

  // resize and reset branches only in first permutation
  if (pindex == 0)
    {
      // - get number of permutations
      fTreeVarNPermutations = fFitter->Permutations()->NPermutations(); 

      // reset tree variables 
      fTreeVarLogLikelihood->clear();
      fTreeVarLogLikelihood->assign(fTreeVarNPermutations, 1.e99);

      fTreeVarLogLikelihoodComp_TF_bhad1->clear();
      fTreeVarLogLikelihoodComp_TF_bhad1->assign(fTreeVarNPermutations, 1.e99);

      fTreeVarLogLikelihoodComp_TF_bhad2->clear();
      fTreeVarLogLikelihoodComp_TF_bhad2->assign(fTreeVarNPermutations, 1.e99);

      fTreeVarLogLikelihoodComp_TF_lq1->clear();
      fTreeVarLogLikelihoodComp_TF_lq1->assign(fTreeVarNPermutations, 1.e99);

      fTreeVarLogLikelihoodComp_TF_lq2->clear();
      fTreeVarLogLikelihoodComp_TF_lq2->assign(fTreeVarNPermutations, 1.e99);

      fTreeVarLogLikelihoodComp_TF_lq3->clear();
      fTreeVarLogLikelihoodComp_TF_lq3->assign(fTreeVarNPermutations, 1.e99);

      fTreeVarLogLikelihoodComp_TF_lq4->clear();
      fTreeVarLogLikelihoodComp_TF_lq4->assign(fTreeVarNPermutations, 1.e99);

      fTreeVarLogLikelihoodComp_BW_Whad1->clear();
      fTreeVarLogLikelihoodComp_BW_Whad1->assign(fTreeVarNPermutations, 1.e99);

      fTreeVarLogLikelihoodComp_BW_Whad2->clear();
      fTreeVarLogLikelihoodComp_BW_Whad2->assign(fTreeVarNPermutations, 1.e99);

      fTreeVarLogLikelihoodComp_BW_Thad1->clear();
      fTreeVarLogLikelihoodComp_BW_Thad1->assign(fTreeVarNPermutations, 1.e99);

      fTreeVarLogLikelihoodComp_BW_Thad2->clear();
      fTreeVarLogLikelihoodComp_BW_Thad2->assign(fTreeVarNPermutations, 1.e99);

//      fTreeVarLogLikelihoodComponents->clear();
//      std::vector<double> tempvec(3, 0.0);
//      fTreeVarLogLikelihoodComponents->assign(fTreeVarNPermutations, tempvec);

      fTreeVarIntegral->clear();
      fTreeVarIntegral->assign(fTreeVarNPermutations, -1.);

      fTreeVarEventProbability->clear();
      fTreeVarEventProbability->assign(fTreeVarNPermutations, 0.);

      fTreeVarMinuitStatus->clear();
      fTreeVarMinuitStatus->assign(fTreeVarNPermutations, 0);

      fTreeVarConvergenceStatus->clear();
      fTreeVarConvergenceStatus->assign(fTreeVarNPermutations, 0);

      fTreeVarBestPermutation->clear(); 
      fTreeVarBestPermutation->assign(fTreeVarNPermutations, -1);

      unsigned int n = fTreeVarParameters->size();
      for (unsigned int i = 0; i < n; ++i) {
        std::vector<double>* d = fTreeVarParameters->at(i); 
        d->clear();
        d->assign(fTreeVarNPermutations, -1.);
        
        d = fTreeVarParameterErrors->at(i); 
        d->clear();
        d->assign(fTreeVarNPermutations, -1.);
      }

      n = fTreeVarModel->size();
      for (unsigned int i = 0; i < n; ++i) {
        std::vector<double>* d = fTreeVarModel->at(i); 
        d->clear();
        d->assign(fTreeVarNPermutations, -1.);    
      }
      
      n = fTreeIntVarModel->size();
      for (unsigned int i = 0; i < n; ++i) {
        std::vector<int>* d = fTreeIntVarModel->at(i); 
        d->clear();
        d->assign(fTreeVarNPermutations, -1);    
      }
      
    }
        
  (*fTreeVarLogLikelihood)[pindex] = fFitter->Likelihood()->LogLikelihood( fFitter->Likelihood()->GetBestFitParameters() ); 

//  (*fTreeVarLogLikelihoodComponents)[pindex] = fFitter->Likelihood()->LogLikelihoodComponents( fFitter->Likelihood()->GetBestFitParameters() ); 

  (*fTreeVarLogLikelihoodComp_TF_bhad1)[pindex] = fFitter->Likelihood()->LogLikelihoodComponents( fFitter->Likelihood()->GetBestFitParameters() ).at(0);  
  (*fTreeVarLogLikelihoodComp_TF_bhad2)[pindex] = fFitter->Likelihood()->LogLikelihoodComponents( fFitter->Likelihood()->GetBestFitParameters() ).at(1);  
  (*fTreeVarLogLikelihoodComp_TF_lq1)[pindex] = fFitter->Likelihood()->LogLikelihoodComponents( fFitter->Likelihood()->GetBestFitParameters() ).at(2);  
  (*fTreeVarLogLikelihoodComp_TF_lq2)[pindex] = fFitter->Likelihood()->LogLikelihoodComponents( fFitter->Likelihood()->GetBestFitParameters() ).at(3);  
  (*fTreeVarLogLikelihoodComp_TF_lq3)[pindex] = fFitter->Likelihood()->LogLikelihoodComponents( fFitter->Likelihood()->GetBestFitParameters() ).at(4); 
  (*fTreeVarLogLikelihoodComp_TF_lq4)[pindex] = fFitter->Likelihood()->LogLikelihoodComponents( fFitter->Likelihood()->GetBestFitParameters() ).at(5); 
  (*fTreeVarLogLikelihoodComp_BW_Whad1)[pindex] = fFitter->Likelihood()->LogLikelihoodComponents( fFitter->Likelihood()->GetBestFitParameters() ).at(6); 
  (*fTreeVarLogLikelihoodComp_BW_Whad2)[pindex] = fFitter->Likelihood()->LogLikelihoodComponents( fFitter->Likelihood()->GetBestFitParameters() ).at(7); 
  (*fTreeVarLogLikelihoodComp_BW_Thad1)[pindex] = fFitter->Likelihood()->LogLikelihoodComponents( fFitter->Likelihood()->GetBestFitParameters() ).at(8); 
  (*fTreeVarLogLikelihoodComp_BW_Thad2)[pindex] = fFitter->Likelihood()->LogLikelihoodComponents( fFitter->Likelihood()->GetBestFitParameters() ).at(9); 


  (*fTreeVarMinuitStatus)[pindex] = fFitter->MinuitStatus(); 
  (*fTreeVarConvergenceStatus)[pindex] = fFitter->ConvergenceStatus(); 
  (*fTreeVarIntegral)[pindex] = fFitter->Likelihood()->GetNormalization(); 
  (*fTreeVarEventProbability)[pindex] = exp( fFitter->Likelihood()->LogEventProbability() ); 

  // check event probability for NaN
  if ((*fTreeVarLogLikelihood)[pindex] != (*fTreeVarLogLikelihood)[pindex])
    (*fTreeVarLogLikelihood)[pindex] = double(pindex)* (-1e10); 

  if ((*fTreeVarEventProbability)[pindex] != (*fTreeVarEventProbability)[pindex])
    (*fTreeVarEventProbability)[pindex] = 0.; 

  // normalize event probability 
  double sum = 0; 
  bool flagall = true;
  for (int i = 0; i < fTreeVarNPermutations; ++i)
    {
      if ((*fTreeVarLogLikelihood)[i] < 1e99)
        sum += (*fTreeVarEventProbability)[i];
      else
        flagall = false; 
    }

  if (flagall)
    {
      for (int i = 0; i < fTreeVarNPermutations; ++i)
        {
          (*fTreeVarEventProbability)[i] = (*fTreeVarEventProbability)[i] / sum; 
        }
    }

  // sort for best permutation
  for (int i = 0; i < fTreeVarNPermutations; ++i)
    {
      int counter = 0; 
      for (int j = 0; j < fTreeVarNPermutations; ++j)
        if ((*fTreeVarEventProbability)[i] < (*fTreeVarEventProbability)[j])
          counter++;
      (*fTreeVarBestPermutation)[counter] = i; 
    }

  // loop over all parameters 
  for (int i = 0; i < fFitter->Likelihood()->NParameters(); ++i)
    {
      (*fTreeVarParameters->at(i))[pindex] = fFitter->Likelihood()->GetBestFitParameter(i); 
      (*fTreeVarParameterErrors->at(i))[pindex] = fFitter->Likelihood()->GetBestFitParameterError(i); 
    }
	
	int IntVarcounter = 0;	

  // loop over all particle type 
  for (KLFitter::Particles::ParticleType itype = KLFitter::Particles::kParton; itype <= KLFitter::Particles::kPhoton; ++itype)
    {
      // get particle container
      std::vector <TLorentzVector *>* momcontainer = (*fParticlesModel)->ParticleContainer(itype); 

      // get number of particles in container 
      int n = int(momcontainer->size());                        
                                
      // loop over particles 
      for (int i = 0; i < n; ++i)
        {
          // get variables
          std::vector<int>* index = fTreeIntVarModel->at(IntVarcounter); 
          std::vector<double>* E  = fTreeVarModel->at(counter);
          std::vector<double>* px = fTreeVarModel->at(++counter); 
          std::vector<double>* py = fTreeVarModel->at(++counter); 
          std::vector<double>* pz = fTreeVarModel->at(++counter); 
          std::vector<double>* m = fTreeVarModel->at(++counter); 
          std::vector<double>* pt = fTreeVarModel->at(++counter); 
          std::vector<double>* eta = fTreeVarModel->at(++counter); 
          std::vector<double>* phi = fTreeVarModel->at(++counter); 
          std::vector<double>* btag = fTreeVarModel->at(++counter); 
                                        
          // get four vector 
          TLorentzVector* lv = momcontainer->at(i); 
                                        
          // fill variables
          (*E)[pindex] = lv->E(); 
          (*px)[pindex] = lv->Px(); 
          (*py)[pindex] = lv->Py(); 
          (*pz)[pindex] = lv->Pz(); 
          (*m)[pindex] = lv->M(); 
          (*pt)[pindex] = lv->Pt(); 
          (*eta)[pindex] = lv->Eta(); 
          (*phi)[pindex] = lv->Phi(); 
          if (itype == KLFitter::Particles::kParton && (*fParticlesModel)->JetIndex(i)>=0) {
            (*btag)[pindex] = (*fParticlesModel)->IsBTagged(i);
//            //std::cout << pindex << " : " << (fFitter->Permutations()->PermutationTable())->at(pindex)->at(i) << std::endl;	
//
//            //std::cout <<fParticles->NameParticle(0, KLFitter::Particles::kParton)	<< std::endl;
            (*index)[pindex] = (fFitter->Permutations()->PermutationTable())->at(pindex)->at(i); //(*fParticlesModel)->JetIndex(i);

          	IntVarcounter++;  
          }

          if (itype == KLFitter::Particles::kElectron) {
            (*index)[pindex] = (*fParticlesModel)->ElectronIndex(i);
            IntVarcounter++;
          }
          if (itype == KLFitter::Particles::kMuon) {
            (*index)[pindex] = (*fParticlesModel)->MuonIndex(i);
            IntVarcounter++;
          }
          if (itype == KLFitter::Particles::kPhoton) {
            (*index)[pindex] = (*fParticlesModel)->PhotonIndex(i);
            IntVarcounter++;
          }

          // increase counter
          counter++; 
        }
    }

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput_Allhadronic::TreeModelDeleteAllButNBestPermutations(unsigned int n)
{
  //Sanity check:
  if(n>unsigned(fTreeVarNPermutations)){
    std::cout<< "KLFitter::InterfaceOutput::TreeModelDeleteAllButNBestPermutations(). Number of permutations to be kept larger than total number of permutations!" << std::endl;
    return 0;
  }
  
  //Define tmp vectors for int, double
  std::vector<int> vi_tmp;
  std::vector<double> vd_tmp;
  
  //Copy every single vector to tmp and refill the vector in the order of the best permutation from temp and resize to n
  // fTreeVarEventNumber and fTreeVarBTags remain untouched
  //LogLikelihood
  vd_tmp.resize(fTreeVarLogLikelihood->size());
  copy(fTreeVarLogLikelihood->begin(), fTreeVarLogLikelihood->end(), vd_tmp.begin()); 
  for (unsigned int i=0; i<n; i++)
    (*fTreeVarLogLikelihood)[i] = vd_tmp[(*fTreeVarBestPermutation)[i]];
  fTreeVarLogLikelihood->resize(n);
  //Integral
  vd_tmp.resize(fTreeVarIntegral->size());
  copy(fTreeVarIntegral->begin(), fTreeVarIntegral->end(), vd_tmp.begin()); 
  for (unsigned int i=0; i<n; i++)
    (*fTreeVarIntegral)[i] = vd_tmp[(*fTreeVarBestPermutation)[i]];
  fTreeVarIntegral->resize(n);
  //EventProbability
  vd_tmp.resize(fTreeVarEventProbability->size());
  copy(fTreeVarEventProbability->begin(), fTreeVarEventProbability->end(), vd_tmp.begin()); 
  for (unsigned int i=0; i<n; i++)
    (*fTreeVarEventProbability)[i] = vd_tmp[(*fTreeVarBestPermutation)[i]];
  fTreeVarEventProbability->resize(n);
  //MinuitStatus
  vd_tmp.resize(fTreeVarMinuitStatus->size());
  copy(fTreeVarMinuitStatus->begin(), fTreeVarMinuitStatus->end(), vd_tmp.begin()); 
  for (unsigned int i=0; i<n; i++)
    (*fTreeVarMinuitStatus)[i] = vd_tmp[(*fTreeVarBestPermutation)[i]];
  fTreeVarMinuitStatus->resize(n);
  //ConvergenceStatus
  vi_tmp.resize(fTreeVarConvergenceStatus->size());
  copy(fTreeVarConvergenceStatus->begin(), fTreeVarConvergenceStatus->end(), vi_tmp.begin()); 
  for (unsigned int i=0; i<n; i++)
    (*fTreeVarConvergenceStatus)[i] = vi_tmp[(*fTreeVarBestPermutation)[i]];
  fTreeVarConvergenceStatus->resize(n);

  //Parameters Vector
  for (unsigned int iPar(0), nPar(fTreeVarParameters->size()); iPar < nPar; ++iPar) {
    vd_tmp.resize((*fTreeVarParameters)[iPar]->size());
    copy((*fTreeVarParameters)[iPar]->begin(), (*fTreeVarParameters)[iPar]->end(), vd_tmp.begin()); 
    for (unsigned int i=0; i<n; i++)
      (*(*fTreeVarParameters)[iPar])[i] = vd_tmp[(*fTreeVarBestPermutation)[i]];
    (*fTreeVarParameters)[iPar]->resize(n);
  }
  //ParameterErrors Vector
  for (unsigned int iPar(0), nPar(fTreeVarParameterErrors->size()); iPar < nPar; ++iPar) {
    vd_tmp.resize((*fTreeVarParameterErrors)[iPar]->size());
    copy((*fTreeVarParameterErrors)[iPar]->begin(), (*fTreeVarParameterErrors)[iPar]->end(), vd_tmp.begin()); 
    for (unsigned int i=0; i<n; i++)
      (*(*fTreeVarParameterErrors)[iPar])[i] = vd_tmp[(*fTreeVarBestPermutation)[i]];
    (*fTreeVarParameterErrors)[iPar]->resize(n);
  }
  //VarModel Vector
  for (unsigned int iPar(0), nPar(fTreeVarModel->size()); iPar < nPar; ++iPar) {
    vd_tmp.resize((*fTreeVarModel)[iPar]->size());
    copy((*fTreeVarModel)[iPar]->begin(), (*fTreeVarModel)[iPar]->end(), vd_tmp.begin()); 
    for (unsigned int i=0; i<n; i++)
      (*(*fTreeVarModel)[iPar])[i] = vd_tmp[(*fTreeVarBestPermutation)[i]];
    (*fTreeVarModel)[iPar]->resize(n);
  }
  //IntVarModel Vector
  for (unsigned int iPar(0), nPar(fTreeIntVarModel->size()); iPar < nPar; ++iPar) {
    vi_tmp.resize((*fTreeIntVarModel)[iPar]->size());
    copy((*fTreeIntVarModel)[iPar]->begin(), (*fTreeIntVarModel)[iPar]->end(), vi_tmp.begin()); 
    for (unsigned int i=0; i<n; i++)
      (*(*fTreeIntVarModel)[iPar])[i] = vi_tmp[(*fTreeVarBestPermutation)[i]];
    (*fTreeIntVarModel)[iPar]->resize(n);
  }

  //BestPermutation -- NEEDS TO BE DONE LAST!
  for (unsigned int i=0; i<n; i++)
    (*fTreeVarBestPermutation)[i] = i;
  fTreeVarBestPermutation->resize(n);
  //NPermutations
  fTreeVarNPermutations = n;

  // no error
  return 1;
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput_Allhadronic::FillTreeMeasured()
{
  // check tree
  if (!fTreeMeasured)
    this->CreateTreeMeasured(); 

  // initialize counter 
  int counter = 0; 

  // fill number of reconstructed objects
  fTreeVarNPartonsMeasured = (*fParticlesMeasured)->NPartons(); 
  fTreeVarNElectronsMeasured = (*fParticlesMeasured)->NElectrons(); 
  fTreeVarNMuonsMeasured = (*fParticlesMeasured)->NMuons(); 
  fTreeVarNPhotonsMeasured = (*fParticlesMeasured)->NPhotons(); 

  // loop over all particle type 
  for (KLFitter::Particles::ParticleType itype = KLFitter::Particles::kParton; itype <= KLFitter::Particles::kPhoton; ++itype)
    {
      // check if the branch should exist
      std::string name = TreeMeasuredBranchName(itype);
      if (name == "")
        continue;

      // get particle container
      std::vector <TLorentzVector *>* momcontainer = (*fParticlesMeasured)->ParticleContainer(itype); 

      // get number of particles in container 
      int n = int(momcontainer->size());                        

      // get variables 
      std::vector<double>* vec_E   = fTreeVarMeasured->at(counter); 
      std::vector<double>* vec_px  = fTreeVarMeasured->at(++counter); 
      std::vector<double>* vec_py  = fTreeVarMeasured->at(++counter); 
      std::vector<double>* vec_pz  = fTreeVarMeasured->at(++counter); 
      std::vector<double>* vec_m   = fTreeVarMeasured->at(++counter); 
      std::vector<double>* vec_pt  = fTreeVarMeasured->at(++counter); 
      std::vector<double>* vec_eta = fTreeVarMeasured->at(++counter); 
      std::vector<double>* vec_phi = fTreeVarMeasured->at(++counter); 

      // clear the particle vectors
      vec_E  ->clear();
      vec_px ->clear();
      vec_py ->clear();
      vec_pz ->clear();
      vec_m  ->clear();
      vec_pt ->clear();
      vec_eta->clear();
      vec_phi->clear();

      // increase counter
      counter++; 

      // loop over particles 
      for (int i = 0; i < n; ++i)
        {
          // get four vector 
          TLorentzVector* lv = momcontainer->at(i); 

          // fill variables
          vec_E  ->push_back( lv->E()   ); 
          vec_px ->push_back( lv->Px()  ); 
          vec_py ->push_back( lv->Py()  ); 
          vec_pz ->push_back( lv->Pz()  ); 
          vec_m  ->push_back( lv->M()   ); 
          vec_pt ->push_back( lv->Pt()  ); 
          vec_eta->push_back( lv->Eta() ); 
          vec_phi->push_back( lv->Phi() ); 

        }
    }

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput_Allhadronic::FillTreeSelected()
{
  // check tree
  if (!fTreeSelected)
    this->CreateTreeSelected(); 

  // initialize counter 
  int counter = 0; 

  // fill number of reconstructed objects
  fTreeVarNPartonsSelected = (*fParticlesSelected)->NPartons(); 
  fTreeVarNElectronsSelected = (*fParticlesSelected)->NElectrons(); 
  fTreeVarNMuonsSelected = (*fParticlesSelected)->NMuons(); 
  fTreeVarNPhotonsSelected = (*fParticlesSelected)->NPhotons(); 

  // loop over all particle type 
  for (KLFitter::Particles::ParticleType itype = KLFitter::Particles::kParton; itype <= KLFitter::Particles::kPhoton; ++itype)
    {
      // check if the branch should exist
      std::string name = TreeMeasuredBranchName(itype);
      if (name == "")
        continue;

      // get particle container
      std::vector <TLorentzVector *>* momcontainer = (*fParticlesSelected)->ParticleContainer(itype); 

      // get number of particles in container 
      int n = int(momcontainer->size());                        

      // get variables 
      std::vector<double>* vec_E    = fTreeVarSelected->at(counter); 
      std::vector<double>* vec_px   = fTreeVarSelected->at(++counter); 
      std::vector<double>* vec_py   = fTreeVarSelected->at(++counter); 
      std::vector<double>* vec_pz   = fTreeVarSelected->at(++counter); 
      std::vector<double>* vec_m    = fTreeVarSelected->at(++counter); 
      std::vector<double>* vec_pt   = fTreeVarSelected->at(++counter); 
      std::vector<double>* vec_eta  = fTreeVarSelected->at(++counter); 
      std::vector<double>* vec_phi  = fTreeVarSelected->at(++counter); 
      std::vector<double>* vec_btag = fTreeVarSelected->at(++counter); 

      // clear the particle vectors
      vec_E   ->clear();
      vec_px  ->clear();
      vec_py  ->clear();
      vec_pz  ->clear();
      vec_m   ->clear();
      vec_pt  ->clear();
      vec_eta ->clear();
      vec_phi ->clear();
      vec_btag->clear();

      // increase counter
      counter++; 

      // loop over particles 
      for (int i = 0; i < n; ++i)
        {
          // get four vector 
          TLorentzVector* lv = momcontainer->at(i); 

          // fill variables
          vec_E  ->push_back( lv->E()   ); 
          vec_px ->push_back( lv->Px()  ); 
          vec_py ->push_back( lv->Py()  ); 
          vec_pz ->push_back( lv->Pz()  ); 
          vec_m  ->push_back( lv->M()   ); 
          vec_pt ->push_back( lv->Pt()  ); 
          vec_eta->push_back( lv->Eta() ); 
          vec_phi->push_back( lv->Phi() ); 
          vec_btag->push_back( (*fParticlesSelected)->IsBTagged(i)  ); 
        }
    }

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput_Allhadronic::FillTreeTruth()
{
  // check tree 
  if (!fTreeTruth)
    {
      // error code 
      int err = 1; 

      // create tree 
      if (fParticlesTruth)
        err = this->CreateTreeTruth(); 
                        
      else 
        return 0; 
    }

  // initialize counter 
  int counter = 0; 

  // loop over all particle type 
  for (KLFitter::Particles::ParticleType itype = KLFitter::Particles::kParton; itype <= KLFitter::Particles::kPhoton; ++itype)
    {
      // get particle container
      std::vector <TLorentzVector *>* momcontainer = (*fParticlesTruth)->ParticleContainer(itype); 

      // get number of particles in container 
      int n = int(momcontainer->size());                        
                        
      // loop over particles 
      for (int i = 0; i < n; ++i)
        {
          // get variables 
          std::vector<double>* vec_E  = fTreeVarTruth->at(counter); 
          std::vector<double>* vec_px = fTreeVarTruth->at(++counter); 
          std::vector<double>* vec_py = fTreeVarTruth->at(++counter); 
          std::vector<double>* vec_pz = fTreeVarTruth->at(++counter); 
          std::vector<double>* vec_m = fTreeVarTruth->at(++counter); 
          std::vector<double>* vec_pt = fTreeVarTruth->at(++counter); 
          std::vector<double>* vec_eta = fTreeVarTruth->at(++counter); 
          std::vector<double>* vec_phi = fTreeVarTruth->at(++counter); 

          // clear the particle vectors
          vec_E  ->clear();
          vec_px ->clear();
          vec_py ->clear();
          vec_pz ->clear();
          vec_m  ->clear();
          vec_pt ->clear();
          vec_eta->clear();
          vec_phi->clear();

          // get four vector 
          TLorentzVector* lv = momcontainer->at(i); 
                                        
          // fill variables
          vec_E->push_back(lv->E()); 
          vec_px->push_back(lv->Px()); 
          vec_py->push_back(lv->Py()); 
          vec_pz->push_back(lv->Pz()); 
          vec_m->push_back(lv->M()); 
          vec_pt->push_back(lv->Pt()); 
          vec_eta->push_back(lv->Eta()); 
          vec_phi->push_back(lv->Phi()); 
          // increase counter
          counter++; 
        }
    }

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput_Allhadronic::FillTreeMatching()
{
  // check tree 
  if (!fTreeMatching)
    {
      // error code 
      int err = 1; 
                        
      // create tree 
      if (fMatchingTool)
        err = this->CreateTreeMatching(); 
                        
      else 
        return 0; 
    }
        
  // fill number of reconstructed objects
  fTreeVarNPartonsSelected = (*fParticlesSelected)->NPartons(); 
  fTreeVarNElectronsSelected = (*fParticlesSelected)->NElectrons(); 
  fTreeVarNMuonsSelected = (*fParticlesSelected)->NMuons(); 
  fTreeVarNPhotonsSelected = (*fParticlesSelected)->NPhotons(); 

  // set matching vector for all partons to [-1,-1,-1,-1]	
  for (unsigned int k = 0; k < fTreeVarMatchedPartons->size(); ++k) {
    std::vector<int>* d = fTreeVarMatchedPartons->at(k); 
    d->clear();
    d->assign(fTreeVarNPartonsSelected, -1);
    // set number of matched partons to -1	
    *fTreeVarNMatchedPartons->at(k) = -1; 	
  }
	
  // loop over partons
  for (int i = 0; i < (*fParticlesTruth)->NPartons(); ++i)
    {
      // get number of matched partons	
      *(fTreeVarNMatchedPartons->at(i)) = fMatchingTool->NMatchedTruth(i, KLFitter::Particles::kParton); 
	
      // loop over reconstructed partons 
      for (int j = 0; j < fTreeVarNPartonsSelected; ++j)
        {
          (*(fTreeVarMatchedPartons->at(i)))[j] = (fMatchingTool->ListMatchedTruth(i, KLFitter::Particles::kParton)).at(j);  
        }
    }
  
  // set matching vector for all electrons to [-1,-1,-1,-1]
  for (unsigned int k = 0; k < fTreeVarMatchedElectrons->size(); ++k) {
    std::vector<int>* d = fTreeVarMatchedElectrons->at(k); 
    d->clear();
    d->assign(fTreeVarNElectronsSelected, -1);
    // set number of matched electrons to -1	
    *fTreeVarNMatchedElectrons->at(k) = -1; 	
  }
  // loop over electrons
  for (int i = 0; i < (*fParticlesTruth)->NElectrons(); ++i)
    {
      // get number of matched electrons
      *(fTreeVarNMatchedElectrons->at(i)) = fMatchingTool->NMatchedTruth(i, KLFitter::Particles::kElectron); 

      // loop over reconstructed electrons 
      for (int j = 0; j < fTreeVarNElectronsSelected; ++j)
        {
          (*(fTreeVarMatchedElectrons->at(i)))[j] = (fMatchingTool->ListMatchedTruth(i, KLFitter::Particles::kElectron)).at(j);  
        }
    }

  // set matching vector for all muons to [-1,-1,-1,-1]
  for (unsigned int k = 0; k < fTreeVarMatchedMuons->size(); ++k) {
    std::vector<int>* d = fTreeVarMatchedMuons->at(k); 
    d->clear();
    d->assign(fTreeVarNMuonsSelected, -1);
    // set number of matched muons to -1	
    *fTreeVarNMatchedMuons->at(k) = -1; 
  }
  // loop over muons
  for (int i = 0; i < (*fParticlesTruth)->NMuons(); ++i)
    {
      // get number of matched muons
      *(fTreeVarNMatchedMuons->at(i)) = fMatchingTool->NMatchedTruth(i, KLFitter::Particles::kMuon); 

      // loop over reconstructed muons 
      for (int j = 0; j < fTreeVarNMuonsSelected; ++j)
        {
          (*(fTreeVarMatchedMuons->at(i)))[j] = (fMatchingTool->ListMatchedTruth(i, KLFitter::Particles::kMuon)).at(j);  
        }
    }

  // set matching vector for all photons to [-1,-1,-1,-1]
  for (unsigned int k = 0; k < fTreeVarMatchedPhotons->size(); ++k) {
    std::vector<int>* d = fTreeVarMatchedPhotons->at(k); 
    d->clear();
    d->assign(fTreeVarNPhotonsSelected, -1);
    // set number of matched photons to -1	
    *fTreeVarNMatchedPhotons->at(k) = -1; 
  } 
  // loop over photons
  for (int i = 0; i < (*fParticlesTruth)->NPhotons(); ++i)
    {

      // get number of matched photons
      *(fTreeVarNMatchedPhotons->at(i)) = fMatchingTool->NMatchedTruth(i, KLFitter::Particles::kPhoton); 

      // loop over reconstructed photons 
      for (int j = 0; j < fTreeVarNPhotonsSelected; ++j)
        {
          (*(fTreeVarMatchedPhotons->at(i)))[j] = (fMatchingTool->ListMatchedTruth(i, KLFitter::Particles::kPhoton)).at(j);  
        }
    }

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput_Allhadronic::FillTreeMap()
{
  // check tree 
  if (!fTreeMap)
    {
      // error code 
      int err = 1; 
                        
      // create tree 
      if (fSelectionTool)
        err = this->CreateTreeMap(); 
                        
      else 
        return 0; 
    }
        
  // fill number of reconstructed objects
  fTreeVarNPartonsSelected = (*fParticlesSelected)->NPartons(); 
  fTreeVarNElectronsSelected = (*fParticlesSelected)->NElectrons(); 
  fTreeVarNMuonsSelected = (*fParticlesSelected)->NMuons(); 
  fTreeVarNPhotonsSelected = (*fParticlesSelected)->NPhotons(); 

  fTreeVarMapJets->clear();
  fTreeVarMapJets->assign(fTreeVarNPartonsSelected, -1);

  fTreeVarMapJets->clear();
  fTreeVarMapJets->assign(fTreeVarNPartonsSelected, -1);

  fTreeVarMapElectrons->clear();
  fTreeVarMapElectrons->assign(fTreeVarNElectronsSelected, -1);

  fTreeVarMapMuons->clear();
  fTreeVarMapMuons->assign(fTreeVarNMuonsSelected, -1);

  fTreeVarMapPhotons->clear();
  fTreeVarMapPhotons->assign(fTreeVarNPhotonsSelected, -1);

  // get maps 
  for (int i = 0; i < fTreeVarNPartonsSelected; ++i)
    (*fTreeVarMapJets)[i] = (fSelectionTool->MapJets()).at(i); 
  for (int i = 0; i < fTreeVarNElectronsSelected; ++i)
    (*fTreeVarMapElectrons)[i] = (fSelectionTool->MapElectrons()).at(i); 
  for (int i = 0; i < fTreeVarNMuonsSelected; ++i)
    (*fTreeVarMapMuons)[i] = (fSelectionTool->MapMuons()).at(i); 
  for (int i = 0; i < fTreeVarNPhotonsSelected; ++i)
    (*fTreeVarMapPhotons)[i] = (fSelectionTool->MapPhotons()).at(i); 

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput_Allhadronic::FillTrees()
{
  // fill tree with truth particles 
  if (fParticlesTruth && fTreeTruth)
    fTreeTruth->Fill(); 

  // fill tree with measured particles 
  if (fParticlesMeasured && fTreeMeasured)
    fTreeMeasured->Fill(); 

  // fill tree with selected particles 
  if (fParticlesSelected && fTreeSelected)
    fTreeSelected->Fill(); 

  // fill tree with model particles 
  if (fParticlesModel && fTreeModel)
    fTreeModel->Fill(); 

  if (fMatchingTool)
    fTreeMatching->Fill(); 

  if (fSelectionTool)
    fTreeMap->Fill(); 

  fTreeVarEventNumber++; 

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceOutput_Allhadronic::SetEventWeight(double weight)
{
  fEventWeight = weight;

  // no error
  return 0;
}
// --------------------------------------------------------- 
int KLFitter::InterfaceOutput_Allhadronic::SetPileupWeight(double weight)
{
  fPileupWeight = weight;

  // no error
  return 0;
}
// --------------------------------------------------------- 

int KLFitter::InterfaceOutput_Allhadronic::SetPhotonType(bool isNotClassified, bool isRadTopProd, bool isHadTopRadDecay, bool isLepTopRadDecay, bool isHadWRadDecay, bool isLepWRadDecay)
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
std::string KLFitter::InterfaceOutput_Allhadronic::ModifyString(std::string str)
{
  int idx;

  while( (idx=str.find_first_of(' ')) >= 0 ) 
    str.replace(idx, 1, "_" );

  while( (idx=str.find_first_of('-')) >= 0 ) 
    str.replace(idx, 1, "_" );

  return str; 
}

// --------------------------------------------------------- 

std::string KLFitter::InterfaceOutput_Allhadronic::TreeMeasuredBranchName(KLFitter::Particles::ParticleType pType)
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

template<class type> void KLFitter::InterfaceOutput_Allhadronic::Resize(std::vector<std::vector<type>* >* v, unsigned int length)
{
  for (unsigned int i = 0; i < v->size(); i++)
    {
      v->at(i)->resize(length);
    }
}
// --------------------------------------------------------- 
