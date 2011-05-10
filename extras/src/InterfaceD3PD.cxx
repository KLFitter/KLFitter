#include "InterfaceD3PD.h" 
#include "Particles.h"

#include <TLorentzVector.h>

#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>

#include <iostream> 
#include <cmath>

// --------------------------------------------------------- 
KLFitter::InterfaceD3PD::InterfaceD3PD()
{
  fTree = 0; 

  EventNumber = 0;
	if(mcevt_weight)
		mcevt_weight = 0;

	topMu_n = 0;
	topMu_index = 0;
	topMu_use = 0;
	topMu_inTrigger = 0;    
  mu_E = 0;  
  mu_px = 0;  
  mu_py = 0;  
  mu_pz = 0;  
  mu_pt = 0;  
  mu_eta = 0;
  mu_phi = 0;  

  topEl_n = 0;
  topEl_index = 0;
	topEl_use = 0;
	topEl_inTrigger = 0;    
  el_E = 0;  
  el_eta = 0;
  el_deteta = 0;  
  el_phi = 0;  

  topJet_n = 0;
  topJet_index = 0;
	topJet_use = 0;
	topJet_inTrigger = 0;   
  jet_E = 0;  
  jet_pt = 0;  
  jet_eta = 0;
  jet_deteta = 0;  
  jet_phi = 0;  
  jet_flavor_weight_SV0 = 0;  

  topMET_et = 0; 
  topMET_phi = 0; 
  topMET_etx = 0; 
  topMET_ety = 0;

  mc_eta = 0;
  mc_phi = 0;
  mc_pt = 0;
  mc_pdgId = 0;
  mcevt_weight = 0;
  mc_m = 0;
  mc_barcode = 0;
  mc_status = 0;
  mc_parents = 0;
  mc_parent_index = 0;
  mc_child_index = 0;
  mc_children = 0;
      
  // needed due to incomplete dictionaries in ROOT (reading vector<bool> from TTree won't work without)
  gROOT->ProcessLine("#include <vector>");
  gROOT->ProcessLine(".L loader.C+");
}

// --------------------------------------------------------- 
KLFitter::InterfaceD3PD::~InterfaceD3PD()
{
}

// --------------------------------------------------------- 
int KLFitter::InterfaceD3PD::NEvents()
{
  if(fTree)
  	return fTree->GetEntries();
  if(fChain)
  	return fChain->GetEntries();	        
  else
    return 0; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceD3PD::OpenRootFile(const char* filename, Option_t* opt)
{
  // define error code 
  int err = 1; 

  // open file 
  err *= KLFitter::InterfaceRoot::OpenRootFile(filename, opt); 

  // connect Root tree 
  err *= this ->ConnectTree("physics"); 

  // return error code 
  return err; 
}
// --------------------------------------------------------- 
int KLFitter::InterfaceD3PD::OpenRootFiles(std::vector<std::string> filenames, Option_t* opt)
{
  // define error code 
  int err = 1; 
	
	fChain =  new TChain("physics");
  // open files
  for(unsigned int i=0; i<filenames.size(); i++){ 
  	err *= KLFitter::InterfaceRoot::OpenRootFile(filenames.at(i).c_str(), opt);
  	fChain->Add(filenames.at(i).c_str());
  }	 
		
  // connect Root tree 
  err *= this ->ConnectChain(fChain); 

  // return error code 
  return err; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceD3PD::ConnectTree(const char* treename)
{
  // check if file exists 
  if (!fRootFile)
    {
      std::cout << "KLFitter::InterfaceD3PD::ConnectTree(). No Root file defined." << std::endl; 
      return 0; 
    } 

  // check if file is open 
  if (!fRootFile->IsOpen())
    { 
      std::cout << "KLFitter::InterfaceD3PD::ConnectTree(). Root file not open."<< std::endl; 
      return 0; 
    }

  // get tree from file 
  fTree = (TTree *) fRootFile->Get(treename); 

  if (!fTree)
    {
      std::cout << "KLFitter::InterfaceD3PD::ConnectTree(). Tree not found." << std::endl; 
      return 0; 
    }

  return ConnectTree(fTree);
}
// --------------------------------------------------------- 
int KLFitter::InterfaceD3PD::ConnectTree(TTree * fTree)
{
  if (!this->fTree) this->fTree = fTree;
  // set branch adresses
  fTree->SetBranchAddress("EventNumber",  &EventNumber); 
  //fTree->SetBranchAddress("mcevt_weight", &fWeight);
  fTree->SetBranchAddress("mcevt_weight", &mcevt_weight);
		
  fTree->SetBranchAddress("topMu_n",  &topMu_n);
  fTree->SetBranchAddress("topMu_index",  &topMu_index);
  fTree->SetBranchAddress("topMu_use",  &topMu_use);
  fTree->SetBranchAddress("topMu_inTrigger",  &topMu_inTrigger);  
  fTree->SetBranchAddress("mu_E",  &mu_E); 
  fTree->SetBranchAddress("mu_px", &mu_px); 
  fTree->SetBranchAddress("mu_py", &mu_py); 
  fTree->SetBranchAddress("mu_pz", &mu_pz); 
  fTree->SetBranchAddress("mu_pt", &mu_pt); 
  fTree->SetBranchAddress("mu_eta", &mu_eta);
  fTree->SetBranchAddress("mu_phi", &mu_phi); 
 
  fTree->SetBranchAddress("topEl_n",  &topEl_n); 
  fTree->SetBranchAddress("topEl_index",  &topEl_index);
  fTree->SetBranchAddress("topEl_use",  &topEl_use);
  fTree->SetBranchAddress("topEl_inTrigger",  &topEl_inTrigger);
  fTree->SetBranchAddress("el_cl_E",  &el_E); 
  fTree->SetBranchAddress("el_tracketa", &el_eta);
  fTree->SetBranchAddress("el_cl_eta", &el_deteta); 
  fTree->SetBranchAddress("el_trackphi", &el_phi); 
  
  fTree->SetBranchAddress("topJet_n",   &topJet_n);
  fTree->SetBranchAddress("topJet_index",   &topJet_index);
  fTree->SetBranchAddress("topJet_use",  &topJet_use);
  fTree->SetBranchAddress("topJet_inTrigger",  &topJet_inTrigger);  
  fTree->SetBranchAddress("jet_E",   &jet_E); 
  fTree->SetBranchAddress("jet_pt",  &jet_pt); 
  fTree->SetBranchAddress("jet_eta", &jet_eta);
  fTree->SetBranchAddress("jet_emscale_eta", &jet_deteta);  
  fTree->SetBranchAddress("jet_phi", &jet_phi); 
  fTree->SetBranchAddress("jet_flavor_weight_SV0", &jet_flavor_weight_SV0); 

  fTree->SetBranchAddress("topMET_et",  &topMET_et); 
  fTree->SetBranchAddress("topMET_phi", &topMET_phi); 
  fTree->SetBranchAddress("topMET_etx", &topMET_etx); 
  fTree->SetBranchAddress("topMET_ety", &topMET_ety); 

  //Truth Variables
  fTree->SetBranchAddress("mc_eta", &mc_eta );
  fTree->SetBranchAddress("mc_phi", &mc_phi );
  fTree->SetBranchAddress("mc_pt",  &mc_pt );
  fTree->SetBranchAddress("mc_pdgId", &mc_pdgId );
  fTree->SetBranchAddress("mcevt_weight", &mcevt_weight );
  fTree->SetBranchAddress("mc_m", &mc_m );
  fTree->SetBranchAddress("mc_barcode", &mc_barcode );
  fTree->SetBranchAddress("mc_status", &mc_status );
  fTree->SetBranchAddress("mc_parents", &mc_parents );
  fTree->SetBranchAddress("mc_parent_index", &mc_parent_index );
  fTree->SetBranchAddress("mc_child_index", &mc_child_index );
  fTree->SetBranchAddress("mc_children", &mc_children );

  // no error     
  return 1;
}
// --------------------------------------------------------- 
int KLFitter::InterfaceD3PD::ConnectChain(TChain * fChain)
{
  if (!this->fChain) this->fChain = fChain;
  // set branch adresses
  fChain->SetBranchAddress("EventNumber",  &EventNumber); 
  //fChain->SetBranchAddress("mcevt_weight", &fWeight);
  fChain->SetBranchAddress("mcevt_weight", &mcevt_weight);
		
  fChain->SetBranchAddress("topMu_n",  &topMu_n);
  fChain->SetBranchAddress("topMu_index",  &topMu_index);
  fChain->SetBranchAddress("topMu_use",  &topMu_use);
  fChain->SetBranchAddress("topMu_inTrigger",  &topMu_inTrigger);  
  fChain->SetBranchAddress("mu_E",  &mu_E); 
  fChain->SetBranchAddress("mu_px", &mu_px); 
  fChain->SetBranchAddress("mu_py", &mu_py); 
  fChain->SetBranchAddress("mu_pz", &mu_pz); 
  fChain->SetBranchAddress("mu_pt", &mu_pt); 
  fChain->SetBranchAddress("mu_eta", &mu_eta);
  fChain->SetBranchAddress("mu_phi", &mu_phi); 
 
  fChain->SetBranchAddress("topEl_n",  &topEl_n); 
  fChain->SetBranchAddress("topEl_index",  &topEl_index);
  fChain->SetBranchAddress("topEl_use",  &topEl_use);
  fChain->SetBranchAddress("topEl_inTrigger",  &topEl_inTrigger);
  fChain->SetBranchAddress("el_cl_E",  &el_E); 
  fChain->SetBranchAddress("el_tracketa", &el_eta);  
  fChain->SetBranchAddress("el_cl_eta", &el_deteta); 
  fChain->SetBranchAddress("el_trackphi", &el_phi); 
  
  fChain->SetBranchAddress("topJet_n",   &topJet_n);
  fChain->SetBranchAddress("topJet_index",   &topJet_index);
  fChain->SetBranchAddress("topJet_use",  &topJet_use);
  fChain->SetBranchAddress("topJet_inTrigger",  &topJet_inTrigger);  
  fChain->SetBranchAddress("jet_E",   &jet_E); 
  fChain->SetBranchAddress("jet_pt",  &jet_pt); 
  fChain->SetBranchAddress("jet_eta", &jet_eta);
  fChain->SetBranchAddress("jet_emscale_eta", &jet_deteta); 
  fChain->SetBranchAddress("jet_phi", &jet_phi); 
  fChain->SetBranchAddress("jet_flavor_weight_SV0", &jet_flavor_weight_SV0); 

  fChain->SetBranchAddress("topMET_et",  &topMET_et); 
  fChain->SetBranchAddress("topMET_phi", &topMET_phi); 
  fChain->SetBranchAddress("topMET_etx", &topMET_etx); 
  fChain->SetBranchAddress("topMET_ety", &topMET_ety); 

  //Truth Variables
  fChain->SetBranchAddress("mc_eta", &mc_eta );
  fChain->SetBranchAddress("mc_phi", &mc_phi );
  fChain->SetBranchAddress("mc_pt",  &mc_pt );
  fChain->SetBranchAddress("mc_pdgId", &mc_pdgId );
  fChain->SetBranchAddress("mcevt_weight", &mcevt_weight );
  fChain->SetBranchAddress("mc_m", &mc_m );
  fChain->SetBranchAddress("mc_barcode", &mc_barcode );
  fChain->SetBranchAddress("mc_status", &mc_status );
  fChain->SetBranchAddress("mc_parents", &mc_parents );
  fChain->SetBranchAddress("mc_parent_index", &mc_parent_index );
  fChain->SetBranchAddress("mc_child_index", &mc_child_index );
  fChain->SetBranchAddress("mc_children", &mc_children );

  // no error     
  return 1;
}

// --------------------------------------------------------- 
int KLFitter::InterfaceD3PD::Event(int index)
{
		
  // check tree 
  if (!fTree && !fChain)
    {
      std::cout << "KLFitter::InterfaceD3PD::GetEvent(). Tree or Chain not defined." << std::endl; 
      return 0; 
    } 

  if(fTree){
  	// check event number
  	if (index < 0 || index >= fTree->GetEntries())
    	{
      	std::cout << "KLFitter::InterfaceD3PD::GetEvent(). Event number negative or too large." << std::endl; 
      	return 0; 
    	} 
   	// get event 
  	fTree->GetEntry(index);
  }
  
  if(fChain){
  	// check event number
  	if (index < 0 || index >= fChain->GetEntries())
    	{
      	std::cout << "KLFitter::InterfaceD3PD::GetEvent(). Event number negative or too large." << std::endl; 
      	return 0; 
    	} 
    // get event 
  	fChain->GetEntry(index);
  } 
  
  // fill particles 
  if (!this->FillParticles())
    return 0; 

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceD3PD::FillParticles()
{
  // delete old particle container
  if (fParticles) 
    delete fParticles; 

  // delete old truth particles container
  if (fParticlesTruth) 
    delete fParticlesTruth; 

  // create new particle container
  fParticles = new KLFitter::Particles(); 
  
  //set weight
    if (mcevt_weight){
      fWeight = mcevt_weight->at(0);
    } else {
    fWeight=1.0;
    }
  
	// fill jets
	for (int i = 0; i < topJet_n; ++i){
		if ( topJet_use->at(i) ){
                  if (jet_E->at(topJet_index->at(i)) <= 0.)
                    continue;
  		TLorentzVector * tlv_tmp = new TLorentzVector(0,0,0,0);
  		tlv_tmp->SetPtEtaPhiE(jet_pt ->at(topJet_index->at(i)) / 1000. , 
  										jet_eta->at(topJet_index->at(i)), 
  										jet_phi->at(topJet_index->at(i)),
  										jet_E  ->at(topJet_index->at(i)) / 1000.);
					
	  	fParticles->AddParticle(tlv_tmp, jet_deteta->at(topJet_index->at(i)), KLFitter::Particles::kParton,"",jet_flavor_weight_SV0->at(topJet_index->at(i)));
      delete tlv_tmp;
		}
	}
	std::sort(fParticles->ParticleContainer(KLFitter::Particles::kParton)->begin(),  fParticles->ParticleContainer(KLFitter::Particles::kParton)->end() , KLFitter::Particles::PtOrder);

	//fill electrons
  for (int i = 0; i < topEl_n; ++i){
  	if ( topEl_use->at(i) ){
          if (el_E->at(topEl_index->at(i)) <= 0.)
            continue;
      TLorentzVector * tlv_tmp = new TLorentzVector(0,0,0,0);
      tlv_tmp->SetPtEtaPhiE((el_E ->at(topEl_index->at(i)) / 1000.) / cosh(el_eta->at(topEl_index->at(i))),
                        el_eta->at(topEl_index->at(i)),
                        el_phi->at(topEl_index->at(i)),
                        el_E  ->at(topEl_index->at(i)) / 1000.);
      fParticles->AddParticle(tlv_tmp, el_deteta->at(topEl_index->at(i)), KLFitter::Particles::kElectron);
      delete tlv_tmp;
		}   																					 
	}
	std::sort(fParticles->ParticleContainer(KLFitter::Particles::kElectron)->begin(),  fParticles->ParticleContainer(KLFitter::Particles::kElectron)->end() , KLFitter::Particles::PtOrder);

  // fill muons
  for (int i = 0; i < topMu_n; ++i){
  	if ( topMu_use->at(i)){ 
          if (mu_E->at(topMu_index->at(i)) <= 0.)
            continue;
      TLorentzVector * tlv_tmp = new TLorentzVector(0,0,0,0);
      tlv_tmp->SetPxPyPzE(mu_px->at(topMu_index->at(i)) / 1000., 
                      mu_py->at(topMu_index->at(i)) / 1000., 
                      mu_pz->at(topMu_index->at(i)) / 1000., 
                      mu_E->at(topMu_index->at(i)) / 1000.);
      fParticles->AddParticle(tlv_tmp, mu_eta->at(topMu_index->at(i)), KLFitter::Particles::kMuon);
      delete tlv_tmp;
    }
	}
	std::sort(fParticles->ParticleContainer(KLFitter::Particles::kMuon)->begin(),  fParticles->ParticleContainer(KLFitter::Particles::kMuon)->end() , KLFitter::Particles::PtOrder);
	
  // check if input is Signal MC
  if (!fFlagWriteSignalMCTruth)
    return 1; 

  // create truth particle container
  fParticlesTruth = new KLFitter::Particles(); 

  //Find the correct truth particle indices
  this->TruthMapper();

  
  // check if event is proper ttbar event
  if (!this->IsProperMCEvent())
    return 1;

  // do not fill mc information if di-leptonic
  if ( (Truth_WplusHad == false) && (Truth_WminusHad == false) )
    return 1; 

  int index_Whad = -1;
  int index_Wlep = -1;
  int index_tophad = -1; 
  int index_toplep = -1; 
  int index_bhad = -1; 
  int index_blep = -1; 
  int index_q1 = -1; 
  int index_q2 = -1; 
  int index_l = -1; 
  int index_nu = -1;
  
  if (Truth_WplusHad) 
    {
      index_Whad = TruthIdx_Wplus; 
      index_tophad = TruthIdx_t; 
      index_bhad = TruthIdx_b; 

      index_Wlep = TruthIdx_Wminus; 
      index_toplep = TruthIdx_tbar; 
      index_blep = TruthIdx_bbar; 
      index_l = TruthIdx_lminus;
      index_nu = TruthIdx_nbar; 
      index_q1 = TruthIdx_QfromWplus; 
      index_q2 = TruthIdx_QbarfromWplus; 

    }
  else
    {
      index_Wlep = TruthIdx_Wplus; 
      index_toplep = TruthIdx_t; 
      index_blep = TruthIdx_b; 
      index_l = TruthIdx_lplus;
      index_nu = TruthIdx_n; 
      index_q1 = TruthIdx_QfromWminus; 
      index_q2 = TruthIdx_QbarfromWminus; 

      index_Whad = TruthIdx_Wminus; 
      index_tophad = TruthIdx_tbar; 
      index_bhad = TruthIdx_bbar; 

    }
  //Create new temp TLorentzVector
  TLorentzVector * tlv_tmp = new TLorentzVector(0,0,0,0);
  tlv_tmp->SetPtEtaPhiM(mc_pt->at(index_bhad) / 1000.,
                      mc_eta->at(index_bhad), 
                      mc_phi->at(index_bhad), 
                      mc_m->at(index_bhad) / 1000.);
  fParticlesTruth->AddParticle(tlv_tmp, KLFitter::Particles::kParton, "hadronic b quark");
  tlv_tmp->SetPtEtaPhiM(mc_pt->at(index_blep) / 1000., 
                      mc_eta->at(index_blep), 
                      mc_phi->at(index_blep), 
                      mc_m->at(index_blep) / 1000.);
  fParticlesTruth->AddParticle(tlv_tmp, KLFitter::Particles::kParton, "leptonic b quark");
  tlv_tmp->SetPtEtaPhiM(mc_pt->at(index_q1) / 1000., 
                      mc_eta->at(index_q1), 
                      mc_phi->at(index_q1), 
                      mc_m->at(index_q1) / 1000.);
  fParticlesTruth->AddParticle(tlv_tmp, KLFitter::Particles::kParton, "light quark 1");
  tlv_tmp->SetPtEtaPhiM(mc_pt->at(index_q2) / 1000., 
                      mc_eta->at(index_q2), 
                      mc_phi->at(index_q2), 
                      mc_m->at(index_q2) / 1000.);
  fParticlesTruth->AddParticle(tlv_tmp, KLFitter::Particles::kParton, "light quark 2"); 

  if (index_l!=-1 && abs(mc_pdgId->at(index_l)) == 11){
    tlv_tmp->SetPtEtaPhiM(mc_pt->at(index_l) / 1000., 
                        mc_eta->at(index_l), 
                        mc_phi->at(index_l), 
                        mc_m->at(index_l) / 1000.);
    fParticlesTruth->AddParticle(tlv_tmp, KLFitter::Particles::kElectron, "electron");
  } 
  else if (index_l!=-1 && abs(mc_pdgId->at(index_l)) == 13){
    tlv_tmp->SetPtEtaPhiM(mc_pt->at(index_l) / 1000.,                                                                                           mc_eta->at(index_l), 
                        mc_phi->at(index_l), 
                        mc_m->at(index_l) / 1000.);
    fParticlesTruth->AddParticle(tlv_tmp, KLFitter::Particles::kMuon, "muon");
  } 
  else if (index_l!=-1 && abs(mc_pdgId->at(index_l)) == 15){
    tlv_tmp->SetPtEtaPhiM(mc_pt->at(index_l) / 1000., 
                        mc_eta->at(index_l), 
                        mc_phi->at(index_l), 
                        mc_m->at(index_l) / 1000.);
    fParticlesTruth->AddParticle(tlv_tmp, KLFitter::Particles::kTau, "tau");
  } 
  
  if (index_nu!=-1){
    tlv_tmp->SetPtEtaPhiM(mc_pt->at(index_nu) / 1000., 
                        mc_eta->at(index_nu), 
                        mc_phi->at(index_nu), 
                        mc_m->at(index_nu) / 1000.);
    fParticlesTruth->AddParticle(tlv_tmp, KLFitter::Particles::kNeutrino, "neutrino");
  } 
  if (index_tophad >= 0) {
    tlv_tmp->SetPtEtaPhiM(mc_pt->at(index_tophad) / 1000., 
                        mc_eta->at(index_tophad), 
                        mc_phi->at(index_tophad), 
                        mc_m->at(index_tophad) / 1000.);
    fParticlesTruth->AddParticle(tlv_tmp, KLFitter::Particles::kParton, "hadronic top quark"); 
  }

  if (index_toplep >= 0) {
    tlv_tmp->SetPtEtaPhiM(mc_pt->at(index_toplep) / 1000., 
                        mc_eta->at(index_toplep), 
                        mc_phi->at(index_toplep), 
                        mc_m->at(index_toplep) / 1000.);
    fParticlesTruth->AddParticle(tlv_tmp, KLFitter::Particles::kParton, "leptonic top quark"); 
  }
  //free memory
  delete tlv_tmp;

  // no error 
  return 1;
}
// --------------------------------------------------------- 

bool KLFitter::InterfaceD3PD::OriginatesFromPDG(int truthIdx,long pdg)
{
  // a helper function needed for the TruthMapper
  //returns true if particle with index truthIdx originates directly from a particle with pgdid 'pdg'
  //===========================================================================
  long ppdg=mc_pdgId->at(truthIdx);
  for (int i=0;i!=6;++i) { //check to most 6 generations
    int NParents=mc_parent_index->at(truthIdx).size();
    if (NParents>0) {
      int truthMotherIdx=mc_parent_index->at(truthIdx).at(0);
      long ParentPDG=mc_pdgId->at(truthMotherIdx);
      if (ParentPDG==pdg) return true; //found pdg
      else if (ParentPDG==ppdg) truthIdx=truthMotherIdx; //particle decayed from itself, keep looking upwards
      else return false; //other parent
    }
  else return false; //no parent
  }
  return false;
}

// --------------------------------------------------------- 

int KLFitter::InterfaceD3PD::TruthMapper(){

  //Init: No particles found yet
  TruthIdx_t = -1;
  TruthIdx_tbar = -1; 
  TruthIdx_b = -1; 
  TruthIdx_bbar = -1; 
  TruthIdx_Wplus = -1; 
  TruthIdx_Wminus = -1; 
  TruthIdx_QfromWplus = -1; 
  TruthIdx_QbarfromWplus = -1; 
  TruthIdx_QfromWminus = -1; 
  TruthIdx_QbarfromWminus = -1; 
  TruthIdx_lplus = -1; 
  TruthIdx_lminus = -1; 
  TruthIdx_n = -1; 
  TruthIdx_nbar = -1;
  Truth_WplusHad = true;
  Truth_WminusHad = true;

   //counters for cross checks
   int Nt      = 0;
   int Ntbar   = 0;
   int Nb      = 0;
   int Nbbar   = 0;
   int NWplus  = 0;
   int NWminus = 0;
   int NQfromWminus    = 0;
   int NQbarfromWminus = 0;
   int NQfromWplus     = 0;
   int NQbarfromWplus  = 0;
   int Nlplus  = 0;
   int Nlminus  = 0;
   int Nn  = 0;
   int Nnbar  = 0;

   bool TruthHERWIGFlag = fFlagIsHerwigMC;  // !!! Define this parameter in your config file/job option. Its crucial for a correct MC truth particle identification!

   for (unsigned int i=0;i!=mc_pt->size();++i) {
      int pdg = (*mc_pdgId)[i];
     
       //TruthIdx
       bool decays_in_itself=false;
       if (mc_child_index->at(i).size()==1) {
           decays_in_itself = (pdg==mc_pdgId->at(mc_child_index->at(i).at(0)));
       }

       // according to Un-ki Yang <ukyang@hep.manchester.ac.uk> the particles with status code 123 or 124 are the ones to be used

       if ((TruthHERWIGFlag && (mc_status->at(i)==123 || mc_status->at(i)==124)) || (!TruthHERWIGFlag && !decays_in_itself)) {

         //-----------
         // top branch
         //-----------
         if (pdg==6  && Nt==0) { //top (ok)
             TruthIdx_t=i;
             ++Nt;
         }
         else if (pdg==5 && Nb==0 && TruthIdx_t!=-1  && OriginatesFromPDG(i,6)) { //bottom (ok)
             TruthIdx_b=i;
             ++Nb;
         }
         else if (pdg==24 && NWplus==0 && TruthIdx_t!=-1 && OriginatesFromPDG(i,6)) { //W+
             TruthIdx_Wplus=i;
             ++NWplus;
             if (!TruthHERWIGFlag) {
           for (unsigned int c=0;c < mc_child_index->at(i).size();++c) {
               int WpChildIdx=mc_child_index->at(i).at(c);
              long Wc_pdg=mc_pdgId->at(WpChildIdx);
             if (abs(Wc_pdg)>10 && abs(Wc_pdg)<17) { //W+ decays in leptons
                 Truth_WplusHad=false;
                 break;
             }
           }
             }
         }         
        else if (!TruthHERWIGFlag && !Truth_WplusHad && (pdg==-11 || pdg==-13 || pdg==-15) && Nlplus==0 &&  OriginatesFromPDG(i,24) ) { //lplus (!HERWIG)
             TruthIdx_lplus=i;
             ++Nlplus;
         }
         else if (TruthHERWIGFlag && (pdg==-11 || pdg==-13 || pdg==-15) && Nlplus==0) { //lplus (HERWIG)
             TruthIdx_lplus=i;
             Truth_WplusHad=false;
             ++Nlplus;
         }
         else if (!TruthHERWIGFlag && !Truth_WplusHad && (pdg==12 || pdg==14 || pdg==16) && Nn==0 &&  OriginatesFromPDG(i,24)) { //neutrino (!HERWIG)
             TruthIdx_n=i;
             ++Nn;
         }
         else if (TruthHERWIGFlag && (pdg==12 || pdg==14 || pdg==16) && Nn==0 &&  OriginatesFromPDG(i,24)) { //neutrino (HERWIG)
             TruthIdx_n=i;
             ++Nn;
         }
         else if (Truth_WplusHad && NWplus==1 && NQfromWplus==0 && (pdg==2 || pdg==4) &&  OriginatesFromPDG(i,24)) { // up/charm from W+
             TruthIdx_QfromWplus=i;
             ++NQfromWplus;
         }
         else if (Truth_WplusHad && NWplus==1 && NQbarfromWplus==0 && (pdg==-1 || pdg==-3) &&  OriginatesFromPDG(i,24)) { // antidown/antistrange from W+
             TruthIdx_QbarfromWplus=i;
             ++NQbarfromWplus;
         }
           //----------------
           // anti-top branch
           //----------------
         if (pdg==-6 && Ntbar==0 ) { //anti top (ok)
             TruthIdx_tbar=i;
             ++Ntbar;
         }
         else if (pdg==-5 && Nbbar==0 && TruthIdx_tbar!=-1  && OriginatesFromPDG(i,-6)) { //anti bottom (ok)
             TruthIdx_bbar=i;
             ++Nbbar;
         }
         else if (pdg==-24 && NWminus==0 && TruthIdx_tbar!=-1 && OriginatesFromPDG(i,-6)) { //W-
             TruthIdx_Wminus=i;
             ++NWminus;
             if (!TruthHERWIGFlag) {
           for (unsigned int c=0;c<mc_child_index->at(i).size();++c) {
               int WmChildIdx=mc_child_index->at(i).at(c);
             long Wc_pdg=mc_pdgId->at(WmChildIdx);
             if (abs(Wc_pdg)>10 && abs(Wc_pdg)<17) { //W- decays in leptons
                 Truth_WminusHad=false;
                 break;
             }
           }
             }
         }         
        else if (!TruthHERWIGFlag && !Truth_WminusHad && (pdg==11 || pdg==13 || pdg==15) && Nlminus==0 &&  OriginatesFromPDG(i,-24)) { //lminus (!HERWIG)
             TruthIdx_lminus=i;
             ++Nlminus;
         }
         else if (TruthHERWIGFlag && (pdg==11 || pdg==13 || pdg==15) && Nlminus==0) { //lminus (HERWIG)
             Truth_WminusHad=false;
             TruthIdx_lminus=i;
             ++Nlminus;
         }
         else if (!TruthHERWIGFlag && !Truth_WminusHad && (pdg==-12 || pdg==-14 || pdg==-16) && Nnbar==0 &&  OriginatesFromPDG(i,-24)) { //anti neutrino (!HERWIG)
             TruthIdx_nbar=i;
             ++Nnbar;
         }
         else if (TruthHERWIGFlag && (pdg==-12 || pdg==-14 || pdg==-16) && Nnbar==0) { //anti neutrino (HERWIG)
             TruthIdx_nbar=i;
             ++Nnbar;
         }
         else if (Truth_WminusHad && NWminus==1 && NQfromWminus==0 && (pdg==1 || pdg==3) &&  OriginatesFromPDG(i,-24)) { // down/strange from W-
             TruthIdx_QfromWminus=i;
             ++NQfromWminus;
         }
         else if (Truth_WminusHad && NWminus==1 && NQbarfromWminus==0 && (pdg==-2 || pdg==-4) &&  OriginatesFromPDG(i,-24)) { // antiup/anticharm from W-
             TruthIdx_QbarfromWminus=i;
             ++NQbarfromWminus;
         }
       }
   } //loop over all particles
  
  //no error
  return 1;
}
// ---------------------------------------------------------
 
bool KLFitter::InterfaceD3PD::IsProperMCEvent(){
//some sanity checks, most of them commented out:
/*
std::cout<< "=======================" << std::endl;
std::cout<< "TruthIdx_t " <<TruthIdx_t << std::endl;
std::cout<< "TruthIdx_tbar "<< TruthIdx_tbar << std::endl; 
std::cout<< "TruthIdx_b "<< TruthIdx_b<<std::endl; 
std::cout<< "TruthIdx_bbar " <<TruthIdx_bbar <<std::endl; 
std::cout<< "TruthIdx_Wplus " <<TruthIdx_Wplus <<std::endl; 
std::cout<< "TruthIdx_Wminus " <<TruthIdx_Wminus<< std::endl; 
std::cout<< "TruthIdx_QfromWplus " <<TruthIdx_QfromWplus <<std::endl; 
std::cout<< "TruthIdx_QbarfromWplus  " <<TruthIdx_QbarfromWplus <<std::endl; 
std::cout<< "TruthIdx_QfromWminus " << TruthIdx_QfromWminus<<std::endl; 
std::cout<< "TruthIdx_QbarfromWminus " <<TruthIdx_QbarfromWminus <<std::endl; 
std::cout<< "TruthIdx_lplus " << TruthIdx_lplus<<std::endl; 
std::cout<< "TruthIdx_lminus " <<TruthIdx_lminus<< std::endl; 
std::cout<< "TruthIdx_n " << TruthIdx_n<<std::endl; 
std::cout<< "TruthIdx_nbar " << TruthIdx_nbar<<std::endl;
std::cout<< "Truth_WplusHad " << Truth_WplusHad<<std::endl; 
std::cout<< "Truth_WminusHad " << Truth_WminusHad<<std::endl;
*/

bool sane = (TruthIdx_t!=-1 && TruthIdx_tbar!=-1 && TruthIdx_b!=-1 && TruthIdx_bbar!=-1 && TruthIdx_Wplus!=-1 && TruthIdx_Wminus!=-1 && //ttbar->W+W-bbbar
       ( (Truth_WplusHad && Truth_WminusHad && TruthIdx_QfromWplus!=-1 && TruthIdx_QbarfromWplus!=-1 && TruthIdx_QfromWminus!=-1 && TruthIdx_QbarfromWminus!=-1) || //alljets
   (Truth_WplusHad && !Truth_WminusHad && TruthIdx_lminus!=-1 && TruthIdx_nbar!=-1 && TruthIdx_QfromWplus!=-1 && TruthIdx_QbarfromWplus!=-1) || //(l^+)+jets 
   (Truth_WminusHad && !Truth_WplusHad && TruthIdx_lplus!=-1 && TruthIdx_n!=-1 && TruthIdx_QfromWminus!=-1 && TruthIdx_QbarfromWminus!=-1) || //(l^-)+jets
   (!Truth_WplusHad && !Truth_WminusHad && TruthIdx_lminus!=-1 && TruthIdx_nbar!=-1 && TruthIdx_lplus!=-1 && TruthIdx_n!=-1)));

//std::cout<<  "Sanity check: " << sane << std::endl;


/*
if (TruthIdx_t!=-1) {
  std::cout << "-------------"<< std::endl;
  std::cout << "Checking top: "<< std::endl;
  std::cout << "mass: "<< mc_m->at(TruthIdx_t)<<std::endl;
  std::cout << "PDGid: "<< mc_pdgId->at(TruthIdx_t)<<std::endl;
}

if (TruthIdx_tbar!=-1) {
  std::cout << "-------------"<< std::endl;
  std::cout << "Checking antitop: "<< std::endl;
  std::cout << "mass: "<< mc_m->at(TruthIdx_tbar)<<std::endl;
  std::cout << "PDGid: "<< mc_pdgId->at(TruthIdx_tbar)<<std::endl;
}

if (TruthIdx_b!=-1) {
  std::cout << "-------------"<< std::endl;
  std::cout << "Checking b: "<< std::endl;
  std::cout << "mass: "<< mc_m->at(TruthIdx_b)<<std::endl;
  std::cout << "PDGid: "<< mc_pdgId->at(TruthIdx_b)<<std::endl;
}

if (TruthIdx_bbar!=-1) {
  std::cout << "-------------"<< std::endl;
  std::cout << "Checking bbar: "<< std::endl;
  std::cout << "mass: "<< mc_m->at(TruthIdx_bbar)<<std::endl;
  std::cout << "PDGid: "<< mc_pdgId->at(TruthIdx_bbar)<<std::endl;
}

if (TruthIdx_Wplus!=-1) {
  std::cout << "-------------"<< std::endl;
  std::cout << "Checking Wplus: "<< std::endl;
  std::cout << "mass: "<< mc_m->at(TruthIdx_Wplus)<<std::endl;
  std::cout << "PDGid: "<< mc_pdgId->at(TruthIdx_Wplus)<<std::endl;
}

if (TruthIdx_Wminus!=-1) {
  std::cout << "-------------"<< std::endl;
  std::cout << "Checking Wminus: "<< std::endl;
  std::cout << "mass: "<< mc_m->at(TruthIdx_Wminus)<<std::endl;
  std::cout << "PDGid: "<< mc_pdgId->at(TruthIdx_Wminus)<<std::endl;
}
*/
/*
if (TruthIdx_QfromWplus!=-1 && TruthIdx_QbarfromWplus!=-1) {
  TLorentzVector truffiQ, truffiQbar, truffiW;
  truffiQ.SetPtEtaPhiM((*mc_pt)[TruthIdx_QfromWplus]/1000, (*mc_eta)[TruthIdx_QfromWplus], (*mc_phi)[TruthIdx_QfromWplus], (*mc_m)[TruthIdx_QfromWplus]/1000);

  truffiQbar.SetPtEtaPhiM((*mc_pt)[TruthIdx_QbarfromWplus]/1000, (*mc_eta)[TruthIdx_QbarfromWplus], (*mc_phi)[TruthIdx_QbarfromWplus], (*mc_m)[TruthIdx_QbarfromWplus]/1000);

  truffiW = truffiQ + truffiQbar;

  std::cout << "-------------"<< std::endl;
  std::cout << "Wplus reco possible! Mass: "<< truffiW.M() <<std::endl;
}

if (TruthIdx_QfromWminus!=-1 && TruthIdx_QbarfromWminus!=-1) {
  TLorentzVector truffiQ, truffiQbar, truffiW;
  truffiQ.SetPtEtaPhiM((*mc_pt)[TruthIdx_QfromWminus]/1000, (*mc_eta)[TruthIdx_QfromWminus], (*mc_phi)[TruthIdx_QfromWminus], (*mc_m)[TruthIdx_QfromWminus]/1000);

  truffiQbar.SetPtEtaPhiM((*mc_pt)[TruthIdx_QbarfromWminus]/1000, (*mc_eta)[TruthIdx_QbarfromWminus], (*mc_phi)[TruthIdx_QbarfromWminus], (*mc_m)[TruthIdx_QbarfromWminus]/1000);

  truffiW = truffiQ + truffiQbar;

  std::cout << "-------------"<< std::endl;
  std::cout << "Wminus reco possible! Mass: "<< truffiW.M() <<std::endl;
}
*/
return sane;
}
// --------------------------------------------------------- 
