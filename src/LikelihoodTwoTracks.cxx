/*
 * Copyright (c) 2009--2018, the KLFitter developer team
 *
 * This file is part of KLFitter.
 *
 * KLFitter is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at
 * your option) any later version.
 *
 * KLFitter is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with KLFitter. If not, see <http://www.gnu.org/licenses/>.
 */


#include "KLFitter/LikelihoodTwoTracks.h" 
#include "KLFitter/ResolutionBase.h"
#include "KLFitter/Particles.h"
#include "KLFitter/Permutations.h"
#include "KLFitter/PhysicsConstants.h"
#include "KLFitter/DetectorBase.h"

#include <iostream> 
#include <algorithm> 

#include <BAT/BCMath.h> 
#include "BAT/BCParameter.h"

// --------------------------------------------------------- 
KLFitter::LikelihoodTwoTracks::LikelihoodTwoTracks() : KLFitter::LikelihoodBase::LikelihoodBase()
                                                            , pion_mass(139.57)
                                                            , kshort_mass(497.6)
                                                            , kshort_width(7.351e-12)
{

  // define model particles 
  this->DefineModelParticles(); 

  // define parameters 
  this->DefineParameters(); 

  fFlagIntegrate = 0;
}

// --------------------------------------------------------- 
KLFitter::LikelihoodTwoTracks::~LikelihoodTwoTracks()  = default;


double KLFitter::LikelihoodTwoTracks::Log3DGaus(double x0, double x1, double x2, double mean0, double mean1, double mean2, double sigma00, double sigma10, double sigma11, double sigma20, double sigma21, double sigma22)
{
    
    double det = sigma00*sigma11*sigma22+2*(sigma10*sigma21*sigma20)-sigma11*sigma20*sigma20-sigma00*sigma21*sigma21-sigma22*sigma10*sigma10;
    
    double inverse = (1./det)*(
    (x0-mean0)*(x0-mean0)*(sigma11*sigma22-sigma21*sigma21)+
    2*(x0-mean0)*(x1-mean1)*(sigma20*sigma21-sigma10*sigma22)+
    2*(x0-mean0)*(x2-mean2)*(sigma10*sigma21-sigma11*sigma20)+
    (x1-mean1)*(x1-mean1)*(sigma00*sigma22-sigma20*sigma20)+
    2*(x1-mean1)*(x2-mean2)*(sigma10*sigma20-sigma00*sigma21)+
    (x2-mean2)*(x2-mean2)*(sigma00*sigma11-sigma10*sigma10)
    );

    return -0.5*(log(det)+inverse +3+log(2*M_PI));

}


// --------------------------------------------------------- 
int KLFitter::LikelihoodTwoTracks::DefineModelParticles()
{

  // create the particles of the model
  fParticlesModel.reset(new KLFitter::Particles{});

  // add model particles
  //create dummy TLorentzVector
  TLorentzVector * dummy = new TLorentzVector(0,0,0,0); // 4-vector
  std::vector<double> * moredummy = new std::vector<double>(0);

  fParticlesModel->AddParticle(dummy,
                               KLFitter::Particles::kBoson, 
                               "Kshort"); 
  
  fParticlesModel->AddParticle(dummy,
                               KLFitter::Particles::kTrack,
                               "pi plus",0,moredummy); 

  fParticlesModel->AddParticle(dummy,
                               KLFitter::Particles::kTrack,
                               "pi minus",1,moredummy); 

  //free memory
  delete dummy; 
  delete moredummy;

  // no error 
  return 1;
}

// --------------------------------------------------------- 
void KLFitter::LikelihoodTwoTracks::DefineParameters()
{
  // add parameters of model
  AddParameter("phi pi plus",                       -M_PI, M_PI);                             // parPiPlusPhi
  AddParameter("theta pi plus",                       0.0, M_PI);                             // parPiPlusTheta
  AddParameter("p pi plus",                       0.0, 100000.0);                             // parPiPlusP
  AddParameter("phi pi minus",                      -M_PI, M_PI);                             // parPiMinusPhi
  AddParameter("theta pi minus",                      0.0, M_PI);                             // parPiMinusTheta
  AddParameter("p pi minus",                      0.0, 100000.0);                             // parPiMinusP
  AddParameter("m Kshort",   kshort_mass-250., kshort_mass+250.);                             // parKShortM
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTwoTracks::CalculateLorentzVectors(std::vector <double> const& parameters)
{

/*
  static double t1_fit_phi;
  static double t1_fit_theta;
  static double t1_fit_p;
  static double t1_fit_m;

  static double t2_fit_phi;
  static double t2_fit_theta;
  static double t2_fit_p;
  static double t2_fit_m;

  static double ks_fit_m;*/

  TLorentzVector * t1 = new TLorentzVector();
  TLorentzVector * t2 = new TLorentzVector();
  TLorentzVector * Ks = new TLorentzVector();

  t1->SetPtEtaPhiM(sin(parameters[parPiPlusTheta])*parameters[parPiPlusP],-log(tan(parameters[parPiPlusTheta]/2.)),parameters[parPiPlusPhi],pion_mass);

  t2->SetPtEtaPhiM(sin(parameters[parPiMinusTheta])*parameters[parPiMinusP],-log(tan(parameters[parPiMinusTheta]/2.)),parameters[parPiMinusPhi],pion_mass);

  (*Ks) = (*t1) + (*t2); 

  t1_fit_phi = parameters[parPiPlusPhi];
  t1_fit_theta = parameters[parPiPlusTheta];
  t1_fit_p = parameters[parPiPlusP];
  t1_fit_m = pion_mass;
  t2_fit_phi = parameters[parPiMinusPhi];
  t2_fit_theta = parameters[parPiMinusTheta];
  t2_fit_p = parameters[parPiMinusP];
  t2_fit_m = pion_mass;
  ks_fit_m = Ks->M();

  delete t1; delete t2; delete Ks;  

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTwoTracks::Initialize()
{
  // error code 
  int err = 1; 

  // save the current permuted particles
  err *= SavePermutedParticles();

  // set initial values
  // (only for Markov chains - initial parameters for other minimisation methods are set in Fitter.cxx)
  SetInitialParameters( GetInitialParameters() );       

  // return error code 
  return err; 
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTwoTracks::RemoveInvariantParticlePermutations()
{
  // error code 
  int err = 1; 

  // remove the permutation from the second and the third jet
  KLFitter::Particles::ParticleType ptype = KLFitter::Particles::kTrack;
  std::vector<int> indexVector_Tracks;
  indexVector_Tracks.push_back(0);
  indexVector_Tracks.push_back(1);
  err *= (*fPermutations)->InvariantParticlePermutations(ptype, indexVector_Tracks); 
        
  // return error code 
  return err; 
}


// --------------------------------------------------------- 
double KLFitter::LikelihoodTwoTracks::LogLikelihood(const std::vector<double> & parameters)
{
  // calculate 4-vectors 
  CalculateLorentzVectors(parameters); 

  // define log of likelihood 
  double logprob(0.); 

  logprob += Log3DGaus(t1_fit_phi,t1_fit_theta,1/t1_fit_p,t1_meas_phi,t1_meas_theta,1/t1_meas_p,t1_meas_sigma00,t1_meas_sigma10,t1_meas_sigma11,t1_meas_sigma20,t1_meas_sigma21,t1_meas_sigma22);

  logprob += Log3DGaus(t2_fit_phi,t2_fit_theta,1/t2_fit_p,t2_meas_phi,t2_meas_theta,1/t2_meas_p,t2_meas_sigma00,t2_meas_sigma10,t2_meas_sigma11,t2_meas_sigma20,t2_meas_sigma21,t2_meas_sigma22);

  logprob += BCMath::LogBreitWignerRel(ks_fit_m,kshort_mass,kshort_width);

  // return log of likelihood 
  return logprob; 
}

// --------------------------------------------------------- 
std::vector<double> KLFitter::LikelihoodTwoTracks::GetInitialParameters()
{
  std::vector<double> values(GetNParameters());


  TLorentzVector * t1 = new TLorentzVector();
  TLorentzVector * t2 = new TLorentzVector();
  TLorentzVector * Ks = new TLorentzVector();

  t1->SetPtEtaPhiM(sin(t1_meas_theta)*t1_fit_p,-log(tan(t1_meas_theta/2.)),t1_meas_phi,pion_mass);

  t2->SetPtEtaPhiM(sin(t2_meas_theta)*t1_fit_p,-log(tan(t2_meas_theta/2.)),t2_meas_phi,pion_mass);

  (*Ks) = (*t1) + (*t2); 

  values[parPiPlusPhi] = t1_meas_phi;
  values[parPiPlusTheta] = t1_meas_theta;
  values[parPiPlusP]  = t1_meas_p;
  values[parPiMinusPhi]  = t2_meas_phi;
  values[parPiMinusTheta] = t2_meas_theta;
  values[parPiMinusP]  = t2_meas_p;
  values[parKShortM] = Ks->M();


  delete t1; delete t2; delete Ks;

  // return the vector
  return values;
}


// --------------------------------------------------------- 
int KLFitter::LikelihoodTwoTracks::SavePermutedParticles() {

  t1_meas_phi = (*fParticlesPermuted)->Track(0)->Phi();
  t1_meas_theta = (*fParticlesPermuted)->Track(0)->Theta();
  t1_meas_p = (*fParticlesPermuted)->Track(0)->P();

  t1_meas_sigma00 =(((*fParticlesPermuted)->Uncertainties(0,KLFitter::Particles::kTrack)))[0]; 
  t1_meas_sigma10 =(((*fParticlesPermuted)->Uncertainties(0,KLFitter::Particles::kTrack)))[1];
  t1_meas_sigma11 =(((*fParticlesPermuted)->Uncertainties(0,KLFitter::Particles::kTrack)))[2];
  t1_meas_sigma20 =(((*fParticlesPermuted)->Uncertainties(0,KLFitter::Particles::kTrack)))[3];
  t1_meas_sigma21 =(((*fParticlesPermuted)->Uncertainties(0,KLFitter::Particles::kTrack)))[4];
  t1_meas_sigma22 =(((*fParticlesPermuted)->Uncertainties(0,KLFitter::Particles::kTrack)))[5];

  t2_meas_phi = (*fParticlesPermuted)->Track(1)->Phi();
  t2_meas_theta = (*fParticlesPermuted)->Track(1)->Theta();
  t2_meas_p = (*fParticlesPermuted)->Track(1)->P();

  t2_meas_sigma00 =(((*fParticlesPermuted)->Uncertainties(1,KLFitter::Particles::kTrack)))[0]; 
  t2_meas_sigma10 =(((*fParticlesPermuted)->Uncertainties(1,KLFitter::Particles::kTrack)))[1];
  t2_meas_sigma11 =(((*fParticlesPermuted)->Uncertainties(1,KLFitter::Particles::kTrack)))[2];
  t2_meas_sigma20 =(((*fParticlesPermuted)->Uncertainties(1,KLFitter::Particles::kTrack)))[3];
  t2_meas_sigma21 =(((*fParticlesPermuted)->Uncertainties(1,KLFitter::Particles::kTrack)))[4];
  t2_meas_sigma22 =(((*fParticlesPermuted)->Uncertainties(1,KLFitter::Particles::kTrack)))[5];


  // no error
  return 1;
}


// --------------------------------------------------------- 

int KLFitter::LikelihoodTwoTracks::BuildModelParticles() {	
if (GetBestFitParameters().size() > 0) CalculateLorentzVectors(GetBestFitParameters());

  TLorentzVector * t1 = fParticlesModel->Track(0);
  TLorentzVector * t2 = fParticlesModel->Track(1);
  TLorentzVector * Ks = fParticlesModel->Boson(0);


  t1->SetPtEtaPhiM(sin(t1_fit_theta)*t1_fit_p,-log(tan(t1_fit_theta/2.)),t1_fit_phi,t1_fit_m);

  t2->SetPtEtaPhiM(sin(t2_fit_theta)*t1_fit_p,-log(tan(t2_fit_theta/2.)),t2_fit_phi,t2_fit_m);

  (*Ks) = (*t1) + (*t2); 

  // no error
  return 1;
}
// --------------------------------------------------------- 

// --------------------------------------------------------- 
std::vector<double> KLFitter::LikelihoodTwoTracks::LogLikelihoodComponents(std::vector<double> parameters)
{
std::vector<double> vecci;

  // calculate 4-vectors 
  CalculateLorentzVectors(parameters); 

  vecci.push_back(Log3DGaus(t1_fit_phi,t1_fit_theta,t1_fit_p,t1_meas_phi,t1_meas_theta,t1_meas_p,t1_meas_sigma00,t1_meas_sigma10,t1_meas_sigma11,t1_meas_sigma20,t1_meas_sigma21,t1_meas_sigma22)); //comp0
  
  vecci.push_back(Log3DGaus(t2_fit_phi,t2_fit_theta,t2_fit_p,t2_meas_phi,t2_meas_theta,t2_meas_p,t2_meas_sigma00,t2_meas_sigma10,t2_meas_sigma11,t2_meas_sigma20,t2_meas_sigma21,t2_meas_sigma22)); //comp1

  vecci.push_back(BCMath::LogBreitWignerRel(ks_fit_m,kshort_mass,kshort_width));  //comp2

  // return log of likelihood 
  return vecci; 
}

