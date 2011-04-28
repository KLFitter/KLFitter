#include "LikelihoodBase.h" 
#include "PhysicsConstants.h" 
#include "Permutations.h" 
#include "DetectorBase.h" 

#include "TRandom3.h"

#include "BAT/BCLog.h" 

#include <iostream> 
#include <string> 

// --------------------------------------------------------- 
KLFitter::LikelihoodBase::LikelihoodBase(Particles** particles) : BCModel(), 
                                                                  fParticlesPermuted(particles), 
                                                                  fPermutations(0),
                                                                  fParticlesModel(0),
                                                                  fPhysicsConstants(new KLFitter::PhysicsConstants()),  
                                                                  fDetector(0),
                                                                  fEventProbability(std::vector<double>(0)),
                                                                  fCutBTag(1e4),
                                                                  fFlagIntegrate(0),
                                                                  fFlagIsNan(false),
                                                                  Setbtagging(kNotag),
								  fbtagRej(-1),
								  fbtagEff(-1)
  
{
  BCLog::SetLogLevel(BCLog::nothing);
  MCMCGetTRandom3()->SetSeed(123456789);
}

// --------------------------------------------------------- 
KLFitter::LikelihoodBase::~LikelihoodBase()
{
  if (fParticlesModel)
    delete fParticlesModel;

	if (fPhysicsConstants)
		delete fPhysicsConstants;
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodBase::SetPhysicsConstants(KLFitter::PhysicsConstants* physicsconstants)
{
  fPhysicsConstants = physicsconstants; 

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodBase::SetInitialParameters(std::vector<double> const& parameters)
{
  // check number of parameters 
  if (int(parameters.size()) != NParameters())
    {
      std::cout << "KLFitter::SetInitialPosition(). Length of vector does not equal the number of parameters." << std::endl; 
      return 0; 
    }

  // set starting point for MCMC 
  MCMCSetInitialPositions(parameters);          

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodBase::SetParameterRange(int index, double parmin, double parmax)
{
  // check index
  if (index < 0 || index >= NParameters())
    {
      std::cout << " KLFitter::Combinatorics::SetParameterRange(). Index out of range." << std::endl; 
      return 0; 
    }

  // set parameter ranges in BAT 
  GetParameter(index)->SetLowerLimit(parmin); 
  GetParameter(index)->SetUpperLimit(parmax); 
  fMCMCBoundaryMin[index] = parmin; 
  fMCMCBoundaryMax[index] = parmax; 

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodBase::SetDetector(KLFitter::DetectorBase** detector)
{
  // set pointer to pointer of detector 
  fDetector = detector; 

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodBase::SetParticlesPermuted(KLFitter::Particles** particles)
{
  // set pointer to pointer of permuted particles
  fParticlesPermuted  = particles; 

  // no error
  return 1;
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodBase::SetPermutations(KLFitter::Permutations** permutations)
{
  // error code 
  int err = 1; 

  // set pointer to pointer of permutation object
  fPermutations = permutations; 

  // return error code
  return err; 
}

// --------------------------------------------------------- 
double KLFitter::LikelihoodBase::ParMin(int index)
{
  // check index
  if (index < 0 || index >= NParameters()) {
    std::cout << " KLFitter::Combinatorics::ParMin(). Index out of range." << std::endl; 
    return 0; 
  }

  // return parameter range from BAT
  return GetParameter(index)->GetLowerLimit(); 
}

// --------------------------------------------------------- 
double KLFitter::LikelihoodBase::ParMax(int index)
{
  // check index
  if (index < 0 || index >= NParameters())
    {
      std::cout << " KLFitter::Combinatorics::ParMax(). Index out of range." << std::endl; 
      return 0; 
    }

  // return parameter range from BAT
  return GetParameter(index)->GetUpperLimit(); 
}

// --------------------------------------------------------- 
double KLFitter::LikelihoodBase::LogEventProbability()
{
  double logprob = 0; 

  if (Setbtagging != kNotag) {
    //                  double probbtag = BTaggingProbability();
    double probbtag = 1; 
    
    if(Setbtagging == kVeto){

      // loop over all model particles.  calculate the overall b-tagging
      // probability which is the product of all probabilities. 
      for (int i = 0; i < fParticlesModel->NPartons(); ++i){

	// get index of corresponding measured particle. 
	int index = fParticlesModel->JetIndex(i); 
	
	if (index<0) { 
	  continue; 
	}
	if (fParticlesModel->BTaggingProbability(i)==0 &&
	    fParticlesModel->FlavorTag(i)==1) probbtag=0;
      }
      
      if (probbtag > 0)
	logprob += log(probbtag); 
      else
	return -1e99; 
    }
    else if (Setbtagging == kWorkingPoint){
      for (int i = 0; i < fParticlesModel->NPartons(); ++i){

	// get index of corresponding measured particle. 
	int index = fParticlesModel->JetIndex(i); 
	
	if (index<0) { 
	  continue; 
	}
	
	probbtag = -1;

	if(fParticlesModel->BTaggingProbability(i)==0 &&
	   fParticlesModel->FlavorTag(i)==1) probbtag = 0; // Light Tagged
	else if(fParticlesModel->BTaggingProbability(i)==0 &&
	   fParticlesModel->FlavorTag(i)==0) probbtag = 1; // Light Not tagged	
	else if(fParticlesModel->BTaggingProbability(i)==1 &&
	   fParticlesModel->FlavorTag(i)==1) probbtag = 2; // b Tagged
	else if(fParticlesModel->BTaggingProbability(i)==1 &&
	   fParticlesModel->FlavorTag(i)==0) probbtag = 3; // b Not tagged
	else std::cout << " KLFitter::LikelihoodBase::LogEventProbability() : b-tagging association failed! " << std::endl;
	  
	if(fbtagRej == -1 || fbtagEff == -1){
	  std::cout <<  " KLFitter::LikelihoodBase::LogEventProbability() : Your working points are not set properly! Returning 0 probability " << std::endl;
	  return -1e99;
	}

	if(probbtag == 0) logprob += log(1./fbtagRej);
	if(probbtag == 1) logprob += log(1 - 1./fbtagRej);
	if(probbtag == 2) logprob += log(fbtagEff);
	if(probbtag == 3) logprob += log(1 - fbtagEff);	
      }            
    }
  }

  // use integrated value of LogLikelihood (default)
  if (fFlagIntegrate)
    logprob += log(GetNormalization()); 
  else
    logprob += LogLikelihood( GetBestFitParameters() ); 
  
  return logprob; 
}

// --------------------------------------------------------- 
double KLFitter::LikelihoodBase::BTaggingProbability()
{
  double prob = 1.0; 

  // get number of partons
  unsigned int npartons = (*fParticlesPermuted)->NPartons();

  // loop over all measured particles. important: assumes that index
  // of measured particle is the same as the index of the model
  // particle
  for (unsigned int i = 0; i < npartons; ++i) {

    // get index of corresponding measured particle. 
    int index = fParticlesModel->JetIndex(i); 
    
    if (index<0) { 
      continue; 
    }

    if (fParticlesModel->BTaggingProbability(i) > 0.5)
      prob *= (*fParticlesPermuted)->BTaggingProbability(i); 
    else
      prob *= (1.0 - (*fParticlesPermuted)->BTaggingProbability(i)); 
  }

  // return log of probability 
  return prob; 
}

// --------------------------------------------------------- 
void KLFitter::LikelihoodBase::CalculateFlavorTags()
{
  // get number of partons 
  unsigned int npartons = fParticlesModel->NPartons();

  // loop over all model particles. 
  for (unsigned int i = 0; i < npartons; ++i)
    {
      // get index of corresponding measured particle. 
      int index = fParticlesModel->JetIndex(i); 

      if (index<0) { 
        continue; 
      }

      if ((*fParticlesPermuted)->BTaggingProbability(index) > fCutBTag) {
        fParticlesModel->SetFlavorTag(i, 1); 
        (*fParticlesPermuted)->SetFlavorTag(index, 1); 
      }
      else {
        fParticlesModel->SetFlavorTag(i, 0); 
        (*fParticlesPermuted)->SetFlavorTag(index, 0); 
      }
    }
}
