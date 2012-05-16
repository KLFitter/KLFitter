#include "LikelihoodTopLeptonJetsUDSep.h" 
#include "ResolutionBase.h"
#include "Particles.h"
#include "Permutations.h"
#include "PhysicsConstants.h"
#include "DetectorBase.h"

#include <iostream> 
#include <algorithm> 

#include <BAT/BCMath.h> 

// --------------------------------------------------------- 
KLFitter::LikelihoodTopLeptonJetsUDSep::LikelihoodTopLeptonJetsUDSep() : KLFitter::LikelihoodTopLeptonJets::LikelihoodTopLeptonJets(),
	fLJetSeparationMethod(KLFitter::LikelihoodTopLeptonJetsUDSep::kNone)

{

  // define model particles 
  this->DefineModelParticles(); 

  // define parameters 
  this->DefineParameters(); 
}

// --------------------------------------------------------- 
KLFitter::LikelihoodTopLeptonJetsUDSep::~LikelihoodTopLeptonJetsUDSep()
{
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTopLeptonJetsUDSep::DefineModelParticles()
{
  // check if model particles and lorentz vector container exist and delete
  if (fParticlesModel) {
    delete fParticlesModel; 
    fParticlesModel = 0;
  }

  // create the particles of the model 
  fParticlesModel = new KLFitter::Particles(); 

  // add model particles
  //create dummy TLorentzVector
  TLorentzVector * dummy = new TLorentzVector(0,0,0,0); // 4-vector
  fParticlesModel->AddParticle(dummy,
                               KLFitter::Particles::kParton, // type 
                               "hadronic b quark",           // name 
                               0,                            // index of corresponding particle 
                               KLFitter::Particles::kB);     // b jet (truth)

  fParticlesModel->AddParticle(dummy,
                               KLFitter::Particles::kParton, 
                               "leptonic b quark",
                               1,                            // index of corresponding particle 
                               KLFitter::Particles::kB);     // b jet (truth)

  fParticlesModel->AddParticle(dummy,
                               KLFitter::Particles::kParton,
                               "light up type quark",
                               2,                            // index of corresponding particle 
                               KLFitter::Particles::kLightUp); // light up type jet (truth)

  fParticlesModel->AddParticle(dummy,
                               KLFitter::Particles::kParton,
                               "light down type quark",
                               3,                            // index of corresponding particle
                               KLFitter::Particles::kLightDown); // light down type jet (truth)
        
  if (fTypeLepton == kElectron) {
    fParticlesModel->AddParticle(dummy,
                                 KLFitter::Particles::kElectron,
                                 "electron"); 
  }
  else if (fTypeLepton == kMuon) {
    fParticlesModel->AddParticle(dummy,
                                 KLFitter::Particles::kMuon,
                                 "muon"); 
  }

  fParticlesModel->AddParticle(dummy,
                               KLFitter::Particles::kNeutrino, 
                               "neutrino"); 
  
  fParticlesModel->AddParticle(dummy,
                               KLFitter::Particles::kBoson, 
                               "hadronic W"); 
  
  fParticlesModel->AddParticle(dummy,
                               KLFitter::Particles::kBoson,
                               "leptonic W"); 

  fParticlesModel->AddParticle(dummy,
                               KLFitter::Particles::kParton,
                               "hadronic top");

  fParticlesModel->AddParticle(dummy,
                               KLFitter::Particles::kParton,
                               "leptonic top");

  //free memory
  delete dummy; 

  // no error 
  return 1;
}

// --------------------------------------------------------- 
void KLFitter::LikelihoodTopLeptonJetsUDSep::DefineParameters()
{
  // rename light quark parameters
  this->GetParameter("energy light quark 1")->SetName("energy light up type quark");
  this->GetParameter("energy light quark 2")->SetName("energy light down type quark");
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTopLeptonJetsUDSep::RemoveInvariantParticlePermutations()
{
  // error code 
  int err = 1; 

  KLFitter::Particles::ParticleType ptype = KLFitter::Particles::kParton;
  std::vector<int> indexVector_Jets;
  //remove invariant jet permutations of notevent jets
  KLFitter::Particles* particles = (*fPermutations)->Particles();
  indexVector_Jets.clear();
  for (int iPartons = 4; iPartons < particles->NPartons(); iPartons++)
    indexVector_Jets.push_back(iPartons);
  err *= (*fPermutations)->InvariantParticlePermutations(ptype, indexVector_Jets);

  // remove the permutation from the other lepton
  if (fTypeLepton == kElectron)
    {
      ptype = KLFitter::Particles::kMuon;
      std::vector<int> indexVector_Muons;
      for (int iMuon = 0; iMuon < particles->NMuons(); iMuon++)
        indexVector_Muons.push_back(iMuon);
      err *= (*fPermutations)->InvariantParticlePermutations(ptype, indexVector_Muons); 
    }
  else if (fTypeLepton == kMuon)
    {
      ptype = KLFitter::Particles::kElectron;
      std::vector<int> indexVector_Electrons;
      for (int iElectron = 0; iElectron < particles->NElectrons(); iElectron++)
        indexVector_Electrons.push_back(iElectron);
      err *= (*fPermutations)->InvariantParticlePermutations(ptype, indexVector_Electrons); 
    }

  // return error code 
  return err; 
}

// --------------------------------------------------------- 
double KLFitter::LikelihoodTopLeptonJetsUDSep::LogEventProbability()
{
  double logprob = 0; 
  if (fBTagMethod != kNotag) {
    double logprobbtag = LogEventProbabilityBTag();
    if (logprobbtag <= -1e99) return -1e99;
    logprob += logprobbtag;
  }
  if (fLJetSeparationMethod != kNone){
    double logprobljetweight = LogEventProbabilityLJetReweight();
    if (logprobljetweight <= -1e99) return -1e99;
    logprob += logprobljetweight;
  }

  // use integrated value of LogLikelihood (default)
  if (fFlagIntegrate)
    logprob += log(GetNormalization()); 
  else
    logprob += LogLikelihood( GetBestFitParameters() ); 
  return logprob; 
}

// --------------------------------------------------------- 
double KLFitter::LikelihoodTopLeptonJetsUDSep::LogEventProbabilityLJetReweight()
{
//	  std::cout <<  " KDEBUG! Extraweight " << std::endl;
  double logprob = 0; 
switch (fLJetSeparationMethod){

case kPermReweight:

  if (!(fUpJetPtHisto && fDownJetPtHisto&& fBJetPtHisto && fUpJetTagWeightHisto && fDownJetTagWeightHisto && fBJetTagWeightHisto)) {
    std::cout <<  " KLFitter::LikelihoodTopLeptonJetsUDSep::LogEventProbabilityLJetReweight() : Histograms were not set properly! " << std::endl;
    return -1e99; 
  }


    for (int i = 0; i < fParticlesModel->NPartons(); ++i){
      // get index of corresponding measured particle. 

      int index = fParticlesModel->JetIndex(i); 

	if (index<0) { 
	  continue; 
	}
        if (!((*fParticlesPermuted)->BTagWeightSet(index))){
	  std::cout <<  " KLFitter::LikelihoodTopLeptonJetsUDSep::LogEventProbabilityLJetReweight() : bTag weight for particle was not set ! " << std::endl;
	  return -1e99;
        }
      KLFitter::Particles::TrueFlavorType trueFlavor = fParticlesModel->TrueFlavor(i);
      if(trueFlavor==KLFitter::Particles::kB) {
  	logprob += log(BJetPt((*fParticlesPermuted)->Parton(index)->Pt())); 	  
  	logprob += log(BJetTagWeight((*fParticlesPermuted)->BTagWeight(index))); 
//std::cout<<"DEBUG! adding pT weight for b: "<<BJetPt((*fParticlesPermuted)->Parton(index)->Pt())<<std::endl;
//std::cout<<"DEBUG! adding tag weight for b: "<<BJetTagWeight((*fParticlesPermuted)->BTagWeight(index))<<std::endl;
      }
      if(trueFlavor==KLFitter::Particles::kLightUp) {
  	logprob += log(UpJetPt((*fParticlesPermuted)->Parton(index)->Pt())); 	  
  	logprob += log(UpJetTagWeight((*fParticlesPermuted)->BTagWeight(index))); 
//std::cout<<"DEBUG! adding pT weight for b: "<<UpJetPt((*fParticlesPermuted)->Parton(index)->Pt())<<std::endl;
//std::cout<<"DEBUG! adding tag weight for b: "<<UpJetTagWeight((*fParticlesPermuted)->BTagWeight(index))<<std::endl;
      }
      if(trueFlavor==KLFitter::Particles::kLightDown) {
  	logprob += log(DownJetPt((*fParticlesPermuted)->Parton(index)->Pt())); 
 	logprob += log(DownJetTagWeight((*fParticlesPermuted)->BTagWeight(index))); 
//std::cout<<"DEBUG! adding pT weight for b: "<<DownJetPt((*fParticlesPermuted)->Parton(index)->Pt())<<std::endl;
//std::cout<<"DEBUG! adding tag weight for b: "<<DownJetTagWeight((*fParticlesPermuted)->BTagWeight(index))<<std::endl;
      }
    }
  return logprob; 
break;
//////////////////////////////////////////
case kPermReweight2D:
  if (!(fUpJet2DWeightHisto && fDownJet2DWeightHisto && fBJet2DWeightHisto)) {
    std::cout <<  " KLFitter::LikelihoodTopLeptonJetsUDSep::LogEventProbabilityLJetReweight() : 2D Histograms were not set properly! " << std::endl;
    return -1e99; 
  }

    for (int i = 0; i < fParticlesModel->NPartons(); ++i){
      // get index of corresponding measured particle. 

      int index = fParticlesModel->JetIndex(i); 

	if (index<0) { 
	  continue; 
	}
        if (!((*fParticlesPermuted)->BTagWeightSet(index))){
	  std::cout <<  " KLFitter::LikelihoodTopLeptonJetsUDSep::LogEventProbabilityLJetReweight() : bTag weight for particle was not set ! " << std::endl;
	  return -1e99;
        }
      KLFitter::Particles::TrueFlavorType trueFlavor = fParticlesModel->TrueFlavor(i);
      if(trueFlavor==KLFitter::Particles::kB) {
	logprob += log(BJetProb((*fParticlesPermuted)->BTagWeight(index), (*fParticlesPermuted)->Parton(index)->Pt()));
//std::cout<<"DEBUG! adding prob weight for b: "<<BJetProb((*fParticlesPermuted)->BTagWeight(index), (*fParticlesPermuted)->Parton(index)->Pt())<<std::endl;

      }
      if(trueFlavor==KLFitter::Particles::kLightUp) {
	logprob += log(UpJetProb((*fParticlesPermuted)->BTagWeight(index), (*fParticlesPermuted)->Parton(index)->Pt()));
//std::cout<<"DEBUG! adding prob weight for up: "<<UpJetProb((*fParticlesPermuted)->BTagWeight(index), (*fParticlesPermuted)->Parton(index)->Pt())<<std::endl;
      }
      if(trueFlavor==KLFitter::Particles::kLightDown) {
	logprob += log(DownJetProb((*fParticlesPermuted)->BTagWeight(index), (*fParticlesPermuted)->Parton(index)->Pt()));
//std::cout<<"DEBUG! adding prob weight for down: "<<DownJetProb((*fParticlesPermuted)->BTagWeight(index), (*fParticlesPermuted)->Parton(index)->Pt())<<std::endl;
      }
    }
  return logprob; 
break;

default:
  return logprob; 
break;

}
}

// --------------------------------------------------------- 
double KLFitter::LikelihoodTopLeptonJetsUDSep::LogEventProbabilityBTag()
{
  double logprob = 0; 

    double probbtag = 1; 
    
    if(fBTagMethod == kVeto){
      // loop over all model particles.  calculate the overall b-tagging
      // probability which is the product of all probabilities. 
      for (int i = 0; i < fParticlesModel->NPartons(); ++i){
        // get index of corresponding measured particle.                                                                   
        int index = fParticlesModel->JetIndex(i);
        if (index < 0)
          continue;

        KLFitter::Particles::TrueFlavorType trueFlavor = fParticlesModel->TrueFlavor(i);
        bool isBTagged = fParticlesModel->IsBTagged(i);
	if (((trueFlavor == KLFitter::Particles::kLightUp)||(trueFlavor == KLFitter::Particles::kLightDown)) && isBTagged == true)
          probbtag = 0.;
      }
      
      if (probbtag > 0)
	logprob += log(probbtag); 
      else
	return -1e99; 
    }
    else if (fBTagMethod == kWorkingPoint){
      for (int i = 0; i < fParticlesModel->NPartons(); ++i){
        // get index of corresponding measured particle.                                                                   
        int index = fParticlesModel->JetIndex(i);
        if (index < 0)
          continue;

        KLFitter::Particles::TrueFlavorType trueFlavor = fParticlesModel->TrueFlavor(i);
        bool isBTagged = fParticlesModel->IsBTagged(i);
        double efficiency = fParticlesModel->BTaggingEfficiency(i);
        double rejection = fParticlesModel->BTaggingRejection(i);
	if(rejection < 0 || efficiency < 0){
	  std::cout <<  " KLFitter::LikelihoodBase::LogEventProbability() : Your working points are not set properly! Returning 0 probability " << std::endl;
	  return -1e99;
	}

	if(((trueFlavor == KLFitter::Particles::kLightUp)||(trueFlavor == KLFitter::Particles::kLightDown)) && isBTagged)
          logprob += log(1./rejection);
	else if(((trueFlavor == KLFitter::Particles::kLightUp)||(trueFlavor == KLFitter::Particles::kLightDown)) && !isBTagged)
          logprob += log(1 - 1./rejection);
	else if(trueFlavor == KLFitter::Particles::kB && isBTagged)
          logprob += log(efficiency);
	else if(trueFlavor == KLFitter::Particles::kB && !isBTagged)
          logprob += log(1 - efficiency);
	else
          std::cout << " KLFitter::LikelihoodBase::LogEventProbability() : b-tagging association failed! " << std::endl;
      }            
    }

  return logprob; 
}



// --------------------------------------------------------- 
double KLFitter::LikelihoodTopLeptonJetsUDSep::UpJetPt(double pt)
{
return fUpJetPtHisto->GetBinContent(fUpJetPtHisto->GetXaxis()->FindBin(pt));
}

// --------------------------------------------------------- 
double KLFitter::LikelihoodTopLeptonJetsUDSep::DownJetPt(double pt)
{
return fDownJetPtHisto->GetBinContent(fDownJetPtHisto->GetXaxis()->FindBin(pt));
}

// --------------------------------------------------------- 
double KLFitter::LikelihoodTopLeptonJetsUDSep::BJetPt(double pt)
{
return fBJetPtHisto->GetBinContent(fBJetPtHisto->GetXaxis()->FindBin(pt));
}

// --------------------------------------------------------- 
double KLFitter::LikelihoodTopLeptonJetsUDSep::UpJetTagWeight(double tagweight)
{
return fUpJetTagWeightHisto->GetBinContent(fUpJetTagWeightHisto->GetXaxis()->FindBin(tagweight));
}

// --------------------------------------------------------- 
double KLFitter::LikelihoodTopLeptonJetsUDSep::DownJetTagWeight(double tagweight)
{
return fDownJetTagWeightHisto->GetBinContent(fDownJetTagWeightHisto->GetXaxis()->FindBin(tagweight));
}

// --------------------------------------------------------- 
double KLFitter::LikelihoodTopLeptonJetsUDSep::BJetTagWeight(double tagweight)
{
return fBJetTagWeightHisto->GetBinContent(fBJetTagWeightHisto->GetXaxis()->FindBin(tagweight));
}

// --------------------------------------------------------- 
double KLFitter::LikelihoodTopLeptonJetsUDSep::UpJetProb(double tagweight, double pt)
{
return fUpJet2DWeightHisto->GetBinContent(fUpJet2DWeightHisto->GetXaxis()->FindBin(tagweight), fUpJet2DWeightHisto->GetYaxis()->FindBin(pt));
}

// --------------------------------------------------------- 
double KLFitter::LikelihoodTopLeptonJetsUDSep::DownJetProb(double tagweight, double pt)
{
return fDownJet2DWeightHisto->GetBinContent(fDownJet2DWeightHisto->GetXaxis()->FindBin(tagweight), fDownJet2DWeightHisto->GetYaxis()->FindBin(pt));
}

// --------------------------------------------------------- 
double KLFitter::LikelihoodTopLeptonJetsUDSep::BJetProb(double tagweight, double pt)
{
return fBJet2DWeightHisto->GetBinContent(fBJet2DWeightHisto->GetXaxis()->FindBin(tagweight), fBJet2DWeightHisto->GetYaxis()->FindBin(pt));
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTopLeptonJetsUDSep::LHInvariantPermutationPartner(int iperm, int nperms, int &switchpar1, int &switchpar2)
{
  int partnerid = -1;
  int cache = iperm%6; 	
  switch (nperms)
  {
  case 24: 	
	if ((iperm)%2) {
		partnerid = (iperm -1);
	} else {
		partnerid = (iperm+1);
	}
  break;

  case 120:
	if (cache > 2) {
		partnerid = (iperm -3);
	} else {
		partnerid = (iperm+3);
	}
  break;
	
  default: partnerid = -1;
  }
switchpar1 = 2;
switchpar2 = 3;
return partnerid;
}
