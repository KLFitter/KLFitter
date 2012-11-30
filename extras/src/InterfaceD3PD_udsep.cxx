#include "InterfaceD3PD_udsep.h" 
#include "Particles.h"

#include <TLorentzVector.h>

#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>

#include <iostream> 
#include <cmath>

// --------------------------------------------------------- 
KLFitter::InterfaceD3PD_udsep::InterfaceD3PD_udsep()
{
  fFlagWriteSignalMCTruth = false;
  fSignalMCGen = KLFitter::InterfaceRoot::kHerwig;

  fTree = 0; 

  eventNumber = 0;
  if(eventWeight)
    eventWeight = 0;

  fBtagCut = 0.601713;
  fBtagEff = 0.696;
  fBtagRej = 134;

  mu_E = 0;  
  mu_pt = 0;  
  mu_eta = 0;  
  mu_phi = 0;

  el_E = 0;  
  el_eta = 0;
  el_deteta = 0;  
  el_phi = 0;  

  jet_n = 0;
  jet_E = 0;  
  jet_pt = 0;  
  jet_eta = 0;
  jet_det_eta = 0;  
  jet_phi = 0;  
  jet_MV1 = 0;  

  met_et = 0; 
  met_x = 0; 
  met_y = 0;
  met_sumet = 0;

  mc_eta = 0;
  mc_phi = 0;
  mc_pt = 0;
  mc_pdgId = 0;
  mc_m = 0;
  mc_status = 0;
  mc_parent_index = 0;
  mc_child_index = 0;
      
  // needed due to incomplete dictionaries in ROOT (reading vector<bool> from TTree won't work without)
  gROOT->ProcessLine("#include <vector>");
  gROOT->ProcessLine(".L loader.C+");
}

// --------------------------------------------------------- 
KLFitter::InterfaceD3PD_udsep::~InterfaceD3PD_udsep()
{
}

// --------------------------------------------------------- 
int KLFitter::InterfaceD3PD_udsep::NEvents()
{
  if(fTree)
  	return fTree->GetEntries();
  if(fChain)
  	return fChain->GetEntries();	        
  else
    return 0; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceD3PD_udsep::OpenRootFile(const char* filename, Option_t* opt)
{
  // define error code 
  int err = 1; 

  // open file 
  err *= KLFitter::InterfaceRoot::OpenRootFile(filename, opt); 

  // connect Root tree 
  err *= this ->ConnectTree("GoeTree"); 

  // return error code 
  return err; 
}
// --------------------------------------------------------- 
int KLFitter::InterfaceD3PD_udsep::OpenRootFiles(std::vector<std::string> filenames, Option_t* opt)
{
  // define error code 
  int err = 1; 
	
	fChain =  new TChain("GoeTree");
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
int KLFitter::InterfaceD3PD_udsep::ConnectTree(const char* treename)
{
  // check if file exists 
  if (!fRootFile)
    {
      std::cout << "KLFitter::InterfaceD3PD_udsep::ConnectTree(). No Root file defined." << std::endl; 
      return 0; 
    } 

  // check if file is open 
  if (!fRootFile->IsOpen())
    { 
      std::cout << "KLFitter::InterfaceD3PD_udsep::ConnectTree(). Root file not open."<< std::endl; 
      return 0; 
    }

  // get tree from file 
  fTree = (TTree *) fRootFile->Get(treename); 

  if (!fTree)
    {
      std::cout << "KLFitter::InterfaceD3PD_udsep::ConnectTree(). Tree not found." << std::endl; 
      return 0; 
    }

  return ConnectTree(fTree);
}
// --------------------------------------------------------- 
int KLFitter::InterfaceD3PD_udsep::ConnectTree(TTree * fTree)
{
  if (!this->fTree) this->fTree = fTree;
  // set branch adresses
  fTree->SetBranchAddress("eventNumber",  &eventNumber); 
  fTree->SetBranchAddress("eventWeight", &eventWeight);
		
  fTree->SetBranchAddress("mu_E",  &mu_E); 
  fTree->SetBranchAddress("mu_pt", &mu_pt); 
  fTree->SetBranchAddress("mu_eta", &mu_eta); 
  fTree->SetBranchAddress("mu_phi", &mu_phi);
 
  fTree->SetBranchAddress("el_cl_E",  &el_E); 
  fTree->SetBranchAddress("el_eta", &el_eta);
  fTree->SetBranchAddress("el_cl_eta", &el_deteta); 
  fTree->SetBranchAddress("el_phi", &el_phi); 
  
  fTree->SetBranchAddress("jet_n",   &jet_n);
  fTree->SetBranchAddress("jet_E",   &jet_E); 
  fTree->SetBranchAddress("jet_pt",  &jet_pt); 
  fTree->SetBranchAddress("jet_eta", &jet_eta);
  fTree->SetBranchAddress("jet_det_eta", &jet_det_eta);  
  fTree->SetBranchAddress("jet_phi", &jet_phi); 
  fTree->SetBranchAddress("jet_MV1", &jet_MV1); 

  fTree->SetBranchAddress("met_et",  &met_et); 
  fTree->SetBranchAddress("met_x", &met_x); 
  fTree->SetBranchAddress("met_y", &met_y); 
  fTree->SetBranchAddress("met_sumet", &met_sumet); 

  //Truth Variables
  fTree->SetBranchAddress("mc_eta", &mc_eta );
  fTree->SetBranchAddress("mc_phi", &mc_phi );
  fTree->SetBranchAddress("mc_pt",  &mc_pt );
  fTree->SetBranchAddress("mc_pdgId", &mc_pdgId );
  fTree->SetBranchAddress("eventWeight", &eventWeight );
  fTree->SetBranchAddress("mc_m", &mc_m );
  fTree->SetBranchAddress("mc_status", &mc_status );
  fTree->SetBranchAddress("mc_parent_index", &mc_parent_index );
  fTree->SetBranchAddress("mc_child_index", &mc_child_index );

  // no error     
  return 1;
}
// --------------------------------------------------------- 
int KLFitter::InterfaceD3PD_udsep::ConnectChain(TChain * fChain)
{
  if (!this->fChain) this->fChain = fChain;
   fChain->SetBranchStatus("*", 0);

		const char* branches[] =
		{"eventNumber", "eventWeight", "mu_E","mu_pt","mu_eta","mu_phi","el_cl_E","el_eta","el_cl_eta", "el_phi", "jet_n", "jet_E", "jet_pt","jet_eta", "jet_det_eta", "jet_phi", "jet_MV1", "met_et", "met_x", "met_y", "met_sumet", "mc_eta", "mc_phi", "mc_pt", "mc_pdgId","eventWeight", "mc_m", "mc_status", "mc_parent_index","mc_child_index"};

   for (unsigned int b = 0; b < sizeof(branches) / sizeof(const char*); b++)
   fChain->SetBranchStatus(branches[b], 1);

  // set branch adresses
  fChain->SetBranchAddress("eventNumber",  &eventNumber); 
  fChain->SetBranchAddress("eventWeight", &eventWeight);
		
  fChain->SetBranchAddress("mu_E",  &mu_E); 
  fChain->SetBranchAddress("mu_pt", &mu_pt); 
  fChain->SetBranchAddress("mu_phi", &mu_phi); 
  fChain->SetBranchAddress("mu_eta", &mu_eta);  
 
  fChain->SetBranchAddress("el_cl_E",  &el_E); 
  fChain->SetBranchAddress("el_eta", &el_eta);  
  fChain->SetBranchAddress("el_cl_eta", &el_deteta); 
  fChain->SetBranchAddress("el_phi", &el_phi); 
  
  fChain->SetBranchAddress("jet_n",   &jet_n);
  fChain->SetBranchAddress("jet_E",   &jet_E); 
  fChain->SetBranchAddress("jet_pt",  &jet_pt); 
  fChain->SetBranchAddress("jet_eta", &jet_eta);
  fChain->SetBranchAddress("jet_det_eta", &jet_det_eta); 
  fChain->SetBranchAddress("jet_phi", &jet_phi); 
  fChain->SetBranchAddress("jet_MV1", &jet_MV1); 

  fChain->SetBranchAddress("met_et",  &met_et); 
  fChain->SetBranchAddress("met_x", &met_x); 
  fChain->SetBranchAddress("met_y", &met_y); 
  fChain->SetBranchAddress("met_sumet", &met_sumet); 

  //Truth Variables
  fChain->SetBranchAddress("mc_eta", &mc_eta );
  fChain->SetBranchAddress("mc_phi", &mc_phi );
  fChain->SetBranchAddress("mc_pt",  &mc_pt );
  fChain->SetBranchAddress("mc_pdgId", &mc_pdgId );
  fChain->SetBranchAddress("eventWeight", &eventWeight );
  fChain->SetBranchAddress("mc_m", &mc_m );
  fChain->SetBranchAddress("mc_status", &mc_status );
  fChain->SetBranchAddress("mc_parent_index", &mc_parent_index );
  fChain->SetBranchAddress("mc_child_index", &mc_child_index );

  // no error     
  return 1;
}

// --------------------------------------------------------- 
int KLFitter::InterfaceD3PD_udsep::Event(int index)
{
		
  // check tree 
  if (!fTree && !fChain)
    {
      std::cout << "KLFitter::InterfaceD3PD_udsep::GetEvent(). Tree or Chain not defined." << std::endl; 
      return 0; 
    } 

  if(fTree){
  	// check event number
  	if (index < 0 || index >= fTree->GetEntries())
    	{
      	std::cout << "KLFitter::InterfaceD3PD_udsep::GetEvent(). Event number negative or too large." << std::endl; 
      	return 0; 
    	} 
   	// get event 
  	fTree->GetEntry(index);
  }
  
  if(fChain){
  	// check event number
  	if (index < 0 || index >= fChain->GetEntries())
    	{
      	std::cout << "KLFitter::InterfaceD3PD_udsep::GetEvent(). Event number negative or too large." << std::endl; 
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
int KLFitter::InterfaceD3PD_udsep::FillParticles()
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
    if (eventWeight){
      fWeight = eventWeight;
    } else {
    fWeight=1.0;
    }
  
	// fill jets
	for (int i = 0; i < jet_n; ++i){
          if (jet_E->at(i) <= 0.)
            continue;

          bool isTagged = jet_MV1->at(i) > fBtagCut;

	  TLorentzVector * tlv_tmp = new TLorentzVector(0,0,0,0);
  	tlv_tmp->SetPtEtaPhiE(jet_pt->at(i), jet_eta->at(i), jet_phi->at(i), jet_E->at(i));
	//If mass is negative, manually correct it to 0.
	if (tlv_tmp->M() < 0){
	std::cout << "KLFitter::InterfaceD3PD_udsep::FillParticles(). Jet mass was negative and corrected to 0." << std::endl;
	tlv_tmp->SetPtEtaPhiM(tlv_tmp->Pt(), tlv_tmp->Eta(), tlv_tmp->Phi(), 0); 
	} 
        fParticles->AddParticle(tlv_tmp, jet_det_eta->at(i), KLFitter::Particles::kParton, "", i, isTagged, fBtagEff, fBtagRej, KLFitter::Particles::kNone, jet_MV1->at(i));
    delete tlv_tmp;
	}

	//fill electrons  
  for (int i = 0; i < el_E->size(); ++i){
    if (el_E->at(i) <= 0.)
      continue;

    TLorentzVector * tlv_tmp = new TLorentzVector(0,0,0,0);
    tlv_tmp->SetPtEtaPhiE((el_E->at(i)) / cosh(el_eta->at(i)), el_eta->at(i), el_phi->at(i), el_E->at(i));
	//If mass is negative, manually correct it to 0.
	if (tlv_tmp->M() < 0){
	std::cout << "KLFitter::InterfaceD3PD_udsep::FillParticles(). Electron mass was negative and corrected to 0." << std::endl;
	tlv_tmp->SetPtEtaPhiM(tlv_tmp->Pt(), tlv_tmp->Eta(), tlv_tmp->Phi(), 0); 
	} 
    fParticles->AddParticle(tlv_tmp, el_deteta->at(i), KLFitter::Particles::kElectron, "", i);
    delete tlv_tmp;
	}

  // fill muons
  for (int i = 0; i < mu_E->size(); ++i){
    if (mu_E->at(i) <= 0.)
      continue;

    TLorentzVector * tlv_tmp = new TLorentzVector(0,0,0,0);
    tlv_tmp->SetPtEtaPhiE(mu_pt->at(i), mu_eta->at(i), mu_phi->at(i), mu_E->at(i));
	//If mass is negative, manually correct it to 0.
	if (tlv_tmp->M() < 0){
	std::cout << "KLFitter::InterfaceD3PD_udsep::FillParticles(). Muon mass was negative and corrected to 0." << std::endl;
	tlv_tmp->SetPtEtaPhiM(tlv_tmp->Pt(), tlv_tmp->Eta(), tlv_tmp->Phi(), 0); 
	} 
    fParticles->AddParticle(tlv_tmp, mu_eta->at(i), KLFitter::Particles::kMuon, "", i);
    delete tlv_tmp;
  }
	
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
      index_q1 = TruthIdx_QfromWplus;		//q1: up type quark 
      index_q2 = TruthIdx_QbarfromWplus;	//q2: down type quark 

    }
  else
    {
      index_Wlep = TruthIdx_Wplus; 
      index_toplep = TruthIdx_t; 
      index_blep = TruthIdx_b; 
      index_l = TruthIdx_lplus;
      index_nu = TruthIdx_n; 
      index_q1 = TruthIdx_QbarfromWminus; 	//q1: up type quark 
      index_q2 = TruthIdx_QfromWminus; 		//q2: down type quark 

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
  fParticlesTruth->AddParticle(tlv_tmp, KLFitter::Particles::kParton, "light up type quark");
  tlv_tmp->SetPtEtaPhiM(mc_pt->at(index_q2) / 1000., 
                      mc_eta->at(index_q2), 
                      mc_phi->at(index_q2), 
                      mc_m->at(index_q2) / 1000.);
  fParticlesTruth->AddParticle(tlv_tmp, KLFitter::Particles::kParton, "light down type quark"); 

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

bool KLFitter::InterfaceD3PD_udsep::OriginatesFromPDG(int truthIdx,long pdg)
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

int KLFitter::InterfaceD3PD_udsep::TruthMapper(){

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

  // !!! Define this parameter in your config file/job option. Its crucial for a correct MC truth particle identification!
  bool TruthHERWIGFlag = true; 
  if (fSignalMCGen == KLFitter::InterfaceRoot::kAcer)
    TruthHERWIGFlag = false;  

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
 
bool KLFitter::InterfaceD3PD_udsep::IsProperMCEvent(){
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

