#include "InterfaceD3PD_Allhadronic.h" 
#include "Particles.h"

#include <TLorentzVector.h>

#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>

#include <iostream> 
#include <cmath>

// --------------------------------------------------------- 
KLFitter::InterfaceD3PD_Allhadronic::InterfaceD3PD_Allhadronic()
{
  fFlagWriteSignalMCTruth = false;
  fSignalMCGen = KLFitter::InterfaceRoot::kHerwig;

  fTree = 0; 

  EventNumber = 0;
  if(mcevt_weight)
    mcevt_weight = 0;

  fBtagCut = 2.20;
  fBtagEff = 0.564;
  fBtagRej = 624.;

  mu_n = 0;
  mu_E = 0;  
  mu_px = 0;  
  mu_py = 0;  
  mu_pz = 0;
  mu_eta = 0;  

  el_n = 0;  
  el_E = 0;  
  el_eta = 0;
  el_deteta = 0;  
  el_phi = 0;  

  jet_n = 0;
  jet_E = 0;  
  jet_pt = 0;  
  jet_eta = 0;
  jet_deteta = 0;  
  jet_phi = 0;  
  jet_flavor_weight_JetFitterCOMBNN = 0;  

  MET_RefFinal_em_tightpp_et = 0; 
  MET_RefFinal_em_tightpp_etx = 0; 
  MET_RefFinal_em_tightpp_ety = 0;
  MET_RefFinal_em_tightpp_sumet = 0;

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
KLFitter::InterfaceD3PD_Allhadronic::~InterfaceD3PD_Allhadronic()
{
}

// --------------------------------------------------------- 
int KLFitter::InterfaceD3PD_Allhadronic::NEvents()
{
  if(fTree)
  	return fTree->GetEntries();
  if(fChain)
  	return fChain->GetEntries();	        
  else
    return 0; 
}

// --------------------------------------------------------- 
int KLFitter::InterfaceD3PD_Allhadronic::OpenRootFile(const char* filename, Option_t* opt)
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
int KLFitter::InterfaceD3PD_Allhadronic::OpenRootFiles(std::vector<std::string> filenames, Option_t* opt)
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
int KLFitter::InterfaceD3PD_Allhadronic::ConnectTree(const char* treename)
{
  // check if file exists 
  if (!fRootFile)
    {
      std::cout << "KLFitter::InterfaceD3PD_Allhadronic::ConnectTree(). No Root file defined." << std::endl; 
      return 0; 
    } 

  // check if file is open 
  if (!fRootFile->IsOpen())
    { 
      std::cout << "KLFitter::InterfaceD3PD_Allhadronic::ConnectTree(). Root file not open."<< std::endl; 
      return 0; 
    }

  // get tree from file 
  fTree = (TTree *) fRootFile->Get(treename); 

  if (!fTree)
    {
      std::cout << "KLFitter::InterfaceD3PD_Allhadronic::ConnectTree(). Tree not found." << std::endl; 
      return 0; 
    }

  return ConnectTree(fTree);
}
// --------------------------------------------------------- 
int KLFitter::InterfaceD3PD_Allhadronic::ConnectTree(TTree * fTree)
{
  if (!this->fTree) this->fTree = fTree;
  // set branch adresses
  fTree->SetBranchAddress("EventNumber",  &EventNumber); 
  fTree->SetBranchAddress("mcevt_weight", &mcevt_weight);
		
  fTree->SetBranchAddress("mu_n",  &mu_n);
  fTree->SetBranchAddress("mu_E",  &mu_E); 
  fTree->SetBranchAddress("mu_px", &mu_px); 
  fTree->SetBranchAddress("mu_py", &mu_py); 
  fTree->SetBranchAddress("mu_pz", &mu_pz);
  fTree->SetBranchAddress("mu_eta", &mu_eta); 
 
  fTree->SetBranchAddress("el_n",  &el_n); 
  fTree->SetBranchAddress("el_cl_E",  &el_E); 
  fTree->SetBranchAddress("el_tracketa", &el_eta);
  fTree->SetBranchAddress("el_cl_eta", &el_deteta); 
  fTree->SetBranchAddress("el_trackphi", &el_phi); 
  
  fTree->SetBranchAddress("jet_n",   &jet_n);
  fTree->SetBranchAddress("jet_E",   &jet_E); 
  fTree->SetBranchAddress("jet_pt",  &jet_pt); 
  fTree->SetBranchAddress("jet_eta", &jet_eta);
  fTree->SetBranchAddress("jet_emscale_eta", &jet_deteta);  
  fTree->SetBranchAddress("jet_phi", &jet_phi); 
  fTree->SetBranchAddress("jet_flavor_weight_JetFitterCOMBNN", &jet_flavor_weight_JetFitterCOMBNN); 

  fTree->SetBranchAddress("MET_RefFinal_em_tightpp_et",  &MET_RefFinal_em_tightpp_et); 
  fTree->SetBranchAddress("MET_RefFinal_em_tightpp_etx", &MET_RefFinal_em_tightpp_etx); 
  fTree->SetBranchAddress("MET_RefFinal_em_tightpp_ety", &MET_RefFinal_em_tightpp_ety); 
  fTree->SetBranchAddress("MET_RefFinal_em_tightpp_sumet", &MET_RefFinal_em_tightpp_sumet); 


  //Truth Variables
  fTree->SetBranchAddress("mc_eta", &mc_eta );
  fTree->SetBranchAddress("mc_phi", &mc_phi );
  fTree->SetBranchAddress("mc_pt",  &mc_pt );
  fTree->SetBranchAddress("mc_pdgId", &mc_pdgId );
  fTree->SetBranchAddress("mcevt_weight", &mcevt_weight );
  fTree->SetBranchAddress("mc_m", &mc_m );
  fTree->SetBranchAddress("mc_status", &mc_status );
  fTree->SetBranchAddress("mc_parent_index", &mc_parent_index );
  fTree->SetBranchAddress("mc_child_index", &mc_child_index );

  // no error     
  return 1;
}
// --------------------------------------------------------- 
int KLFitter::InterfaceD3PD_Allhadronic::ConnectChain(TChain * fChain)
{
  if (!this->fChain) this->fChain = fChain;
   fChain->SetBranchStatus("*", 0);

		const char* branches[] =
		{"EventNumber", "mcevt_weight","mu_n", "mu_E","mu_px","mu_py","mu_pz","mu_eta","el_n","el_cl_E","el_tracketa","el_cl_eta", "el_trackphi", "jet_n", "jet_E", "jet_pt","jet_eta", "jet_emscale_eta", "jet_phi", "jet_flavor_weight_JetFitterCOMBNN", "MET_RefFinal_em_tightpp_et", "MET_RefFinal_em_tightpp_etx", "MET_RefFinal_em_tightpp_ety", "MET_RefFinal_em_tightpp_sumet", "mc_eta", "mc_phi", "mc_pt", "mc_pdgId","mcevt_weight", "mc_m", "mc_status", "mc_parent_index","mc_child_index"};

   for (unsigned int b = 0; b < sizeof(branches) / sizeof(const char*); b++)
   fChain->SetBranchStatus(branches[b], 1);

  // set branch adresses
  fChain->SetBranchAddress("EventNumber",  &EventNumber); 
  fChain->SetBranchAddress("mcevt_weight", &mcevt_weight);
		
  fChain->SetBranchAddress("mu_n",  &mu_n);
  fChain->SetBranchAddress("mu_E",  &mu_E); 
  fChain->SetBranchAddress("mu_px", &mu_px); 
  fChain->SetBranchAddress("mu_py", &mu_py); 
  fChain->SetBranchAddress("mu_pz", &mu_pz);
  fChain->SetBranchAddress("mu_eta", &mu_eta);  
 
  fChain->SetBranchAddress("el_n",  &el_n); 
  fChain->SetBranchAddress("el_cl_E",  &el_E); 
  fChain->SetBranchAddress("el_tracketa", &el_eta);  
  fChain->SetBranchAddress("el_cl_eta", &el_deteta); 
  fChain->SetBranchAddress("el_trackphi", &el_phi); 
  
  fChain->SetBranchAddress("jet_n",   &jet_n);
  fChain->SetBranchAddress("jet_E",   &jet_E); 
  fChain->SetBranchAddress("jet_pt",  &jet_pt); 
  fChain->SetBranchAddress("jet_eta", &jet_eta);
  fChain->SetBranchAddress("jet_emscale_eta", &jet_deteta); 
  fChain->SetBranchAddress("jet_phi", &jet_phi); 
  fChain->SetBranchAddress("jet_flavor_weight_JetFitterCOMBNN", &jet_flavor_weight_JetFitterCOMBNN); 

  fChain->SetBranchAddress("MET_RefFinal_em_tightpp_et",  &MET_RefFinal_em_tightpp_et); 
  fChain->SetBranchAddress("MET_RefFinal_em_tightpp_etx", &MET_RefFinal_em_tightpp_etx); 
  fChain->SetBranchAddress("MET_RefFinal_em_tightpp_ety", &MET_RefFinal_em_tightpp_ety); 
  fChain->SetBranchAddress("MET_RefFinal_em_tightpp_sumet", &MET_RefFinal_em_tightpp_sumet); 

  //Truth Variables
  fChain->SetBranchAddress("mc_eta", &mc_eta );
  fChain->SetBranchAddress("mc_phi", &mc_phi );
  fChain->SetBranchAddress("mc_pt",  &mc_pt );
  fChain->SetBranchAddress("mc_pdgId", &mc_pdgId );
  fChain->SetBranchAddress("mcevt_weight", &mcevt_weight );
  fChain->SetBranchAddress("mc_m", &mc_m );
  fChain->SetBranchAddress("mc_status", &mc_status );
  fChain->SetBranchAddress("mc_parent_index", &mc_parent_index );
  fChain->SetBranchAddress("mc_child_index", &mc_child_index );

  // no error     
  return 1;
}

// --------------------------------------------------------- 
int KLFitter::InterfaceD3PD_Allhadronic::Event(int index)
{
		
  // check tree 
  if (!fTree && !fChain)
    {
      std::cout << "KLFitter::InterfaceD3PD_Allhadronic::GetEvent(). Tree or Chain not defined." << std::endl; 
      return 0; 
    } 

  if(fTree){
  	// check event number
  	if (index < 0 || index >= fTree->GetEntries())
    	{
      	std::cout << "KLFitter::InterfaceD3PD_Allhadronic::GetEvent(). Event number negative or too large." << std::endl; 
      	return 0; 
    	} 
   	// get event 
  	fTree->GetEntry(index);
  }
  
  if(fChain){
  	// check event number
  	if (index < 0 || index >= fChain->GetEntries())
    	{
      	std::cout << "KLFitter::InterfaceD3PD_Allhadronic::GetEvent(). Event number negative or too large." << std::endl; 
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
int KLFitter::InterfaceD3PD_Allhadronic::FillParticles()
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
      fWeight = mcevt_weight->at(0).at(0);
    } else {
    fWeight=1.0;
    }

    // fill jets
    for (int i = 0; i < jet_n; ++i){
      if (jet_E->at(i) <= 0.)
        continue;

      bool isTagged = jet_flavor_weight_JetFitterCOMBNN->at(i) > fBtagCut;

      TLorentzVector * tlv_tmp = new TLorentzVector(0,0,0,0);
  	tlv_tmp->SetPtEtaPhiE(jet_pt->at(i)/1000., jet_eta->at(i), jet_phi->at(i), jet_E->at(i)/1000.);
	//If mass is negative, manually correct it to 0.
	if (tlv_tmp->M() < 0){
	std::cout << "KLFitter::InterfaceD3PD_Allhadronic::FillParticles(). Jet mass was negative and corrected to 0." << std::endl;
	tlv_tmp->SetPtEtaPhiM(tlv_tmp->Pt(), tlv_tmp->Eta(), tlv_tmp->Phi(), 0); 
	} 
        fParticles->AddParticle(tlv_tmp, jet_deteta->at(i), KLFitter::Particles::kParton, "", i, isTagged, fBtagEff, fBtagRej, KLFitter::Particles::kNone, jet_flavor_weight_JetFitterCOMBNN->at(i));
    delete tlv_tmp;
	}

	
     //fill electrons  
     for (int i = 0; i < el_n; ++i){
     if  (el_E->at(i) <= 0.)
       continue;

    TLorentzVector * tlv_tmp = new TLorentzVector(0,0,0,0);
    tlv_tmp->SetPtEtaPhiE((el_E->at(i)/1000.) / cosh(el_eta->at(i)), el_eta->at(i), el_phi->at(i), el_E->at(i)/1000.);
	//If mass is negative, manually correct it to 0.
	if (tlv_tmp->M() < 0){
	std::cout << "KLFitter::InterfaceD3PD_Allhadronic::FillParticles(). Electron mass was negative and corrected to 0." << std::endl;
	tlv_tmp->SetPtEtaPhiM(tlv_tmp->Pt(), tlv_tmp->Eta(), tlv_tmp->Phi(), 0); 
	} 
    fParticles->AddParticle(tlv_tmp, el_deteta->at(i), KLFitter::Particles::kElectron, "", i);
    delete tlv_tmp;
	}

  // fill muons
  for (int i = 0; i < mu_n; ++i){
    if (mu_E->at(i) <= 0.)
      continue;

    TLorentzVector * tlv_tmp = new TLorentzVector(0,0,0,0);
    tlv_tmp->SetPxPyPzE(mu_px->at(i)/1000., mu_py->at(i)/1000., mu_pz->at(i)/1000., mu_E->at(i)/1000.);
	//If mass is negative, manually correct it to 0.
	if (tlv_tmp->M() < 0){
	std::cout << "KLFitter::InterfaceD3PD_Allhadronic::FillParticles(). Muon mass was negative and corrected to 0." << std::endl;
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

  // do not fill mc information if not fully hadronic
  if ( (Truth_WplusHad == false) || (Truth_WminusHad == false) )
      return false; 
  

    int index_Whad1 = -1;
    int index_Whad2 = -1;
    int index_tophad1 = -1; 
    int index_tophad2 = -1; 
    int index_bhad1 = -1; 
    int index_bhad2 = -1; 
    int index_q1 = -1; 
    int index_q2 = -1; 
    int index_q3 = -1; 
    int index_q4 = -1; 
 
 
    index_Whad1 = TruthIdx_Wplus; 
    index_tophad1 = TruthIdx_t; 
    index_bhad1 = TruthIdx_b; 
    index_q1 = TruthIdx_QfromWplus; 
    index_q2 = TruthIdx_QbarfromWplus; 


    index_Whad2 = TruthIdx_Wminus; 
    index_tophad2 = TruthIdx_tbar; 
    index_bhad2 = TruthIdx_bbar; 
    index_q3 = TruthIdx_QfromWminus; 
    index_q4 = TruthIdx_QbarfromWminus;

    //Create new temp TLorentzVector
    TLorentzVector * tlv_tmp = new TLorentzVector(0,0,0,0);
    tlv_tmp->SetPtEtaPhiM(mc_pt->at(index_bhad1) / 1000.,
                      mc_eta->at(index_bhad1), 
                      mc_phi->at(index_bhad1), 
                      mc_m->at(index_bhad1) / 1000.);
    fParticlesTruth->AddParticle(tlv_tmp, KLFitter::Particles::kParton, "hadronic b quark 1");   
    tlv_tmp->SetPtEtaPhiM(mc_pt->at(index_bhad2) / 1000.,
                      mc_eta->at(index_bhad2), 
                      mc_phi->at(index_bhad2), 
                      mc_m->at(index_bhad2) / 1000.);
    fParticlesTruth->AddParticle(tlv_tmp, KLFitter::Particles::kParton, "hadronic b quark 2");    
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
    tlv_tmp->SetPtEtaPhiM(mc_pt->at(index_q3) / 1000., 
                      mc_eta->at(index_q3), 
                      mc_phi->at(index_q3), 
                      mc_m->at(index_q3) / 1000.);
    fParticlesTruth->AddParticle(tlv_tmp, KLFitter::Particles::kParton, "light quark 3"); 
    tlv_tmp->SetPtEtaPhiM(mc_pt->at(index_q4) / 1000., 
                      mc_eta->at(index_q4), 
                      mc_phi->at(index_q4), 
                      mc_m->at(index_q4) / 1000.);
    fParticlesTruth->AddParticle(tlv_tmp, KLFitter::Particles::kParton, "light quark 4");        
 

    if (index_tophad1 >= 0) {
      tlv_tmp->SetPtEtaPhiM(mc_pt->at(index_tophad1) / 1000., 
                        mc_eta->at(index_tophad1), 
                        mc_phi->at(index_tophad1), 
                        mc_m->at(index_tophad1) / 1000.);
      fParticlesTruth->AddParticle(tlv_tmp, KLFitter::Particles::kParton, "hadronic top 1"); 
    }

    if (index_tophad2 >= 0) {
      tlv_tmp->SetPtEtaPhiM(mc_pt->at(index_tophad2) / 1000., 
                        mc_eta->at(index_tophad2), 
                        mc_phi->at(index_tophad2), 
                        mc_m->at(index_tophad2) / 1000.);
      fParticlesTruth->AddParticle(tlv_tmp, KLFitter::Particles::kParton, "hadronic top 2"); 
    }
    //free memory
    delete tlv_tmp;

  // no error 
  return 1;
}
// --------------------------------------------------------- 

bool KLFitter::InterfaceD3PD_Allhadronic::OriginatesFromPDG(int truthIdx,long pdg)
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

int KLFitter::InterfaceD3PD_Allhadronic::TruthMapper(){


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
      else if (TruthHERWIGFlag && (pdg==-11 || pdg==-13 || pdg==-15)) { //lplus (HERWIG)
        Truth_WplusHad=false;
      }
      else if (TruthHERWIGFlag && (pdg==12 || pdg==14 || pdg==16) &&  OriginatesFromPDG(i,24)) { //neutrino (HERWIG)
        Truth_WplusHad=false;
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
      else if (TruthHERWIGFlag && (pdg==11 || pdg==13 || pdg==15)) { //lminus (HERWIG)
        Truth_WminusHad=false;
      }
      else if (TruthHERWIGFlag && (pdg==-12 || pdg==-14 || pdg==-16)) { //anti neutrino (HERWIG)
        Truth_WminusHad=false;
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
 
bool KLFitter::InterfaceD3PD_Allhadronic::IsProperMCEvent(){

   bool sane = (TruthIdx_t!=-1 && TruthIdx_tbar!=-1 && TruthIdx_b!=-1 && TruthIdx_bbar!=-1 && TruthIdx_Wplus!=-1 && TruthIdx_Wminus!=-1 &&  //ttbar->W+W-bbbar
   ((Truth_WplusHad && Truth_WminusHad && TruthIdx_QfromWplus!=-1 && TruthIdx_QbarfromWplus!=-1 && TruthIdx_QfromWminus!=-1 && TruthIdx_QbarfromWminus!=-1)));

    return sane;
}


// --------------------------------------------------------- 

