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
								  fMyParticlesTruth(0),
                                                                  fPhysicsConstants(new KLFitter::PhysicsConstants()),  
                                                                  fDetector(0),
                                                                  fEventProbability(std::vector<double>(0)),
                                                                  fFlagIntegrate(0),
                                                                  fFlagIsNan(false),
                                                                  fBTagMethod(kNotag)

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
int KLFitter::LikelihoodBase::SetInitialParametersNChains(std::vector<double> const& parameters, uint nchains)
{
  // check number of parameters 
  if (int(parameters.size()) != NParameters())
    {
      std::cout << "KLFitter::SetInitialPosition(). Length of vector does not equal the number of parameters." << std::endl; 
      return 0; 
    }
  
  // set starting point for MCMC 
  std::vector< std::vector<double> > par(0.);

  for(uint i=0; i< nchains; ++i) {
    par.push_back(parameters);
  }

  MCMCSetInitialPositions(par);

  MCMCSetFlagInitialPosition(2);

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
int KLFitter::LikelihoodBase::SetMyParticlesTruth(KLFitter::Particles** particles)
{
  // set pointer to pointer of truth particles
  fMyParticlesTruth  = particles; 

  //std::cout << "set particlestruth inside likelihoodbase!" << std::endl;

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

  if (fBTagMethod != kNotag) {
    double logprobbtag = LogEventProbabilityBTag();
    if (logprobbtag <= -1e99) return -1e99;
    logprob += logprobbtag;
  }

  // use integrated value of LogLikelihood (default)
  if (fFlagIntegrate)
    logprob += log(GetNormalization()); 
  else
    logprob += LogLikelihood( GetBestFitParameters() ); 
  
  return logprob; 
}



// --------------------------------------------------------- 
double KLFitter::LikelihoodBase::LogEventProbabilityBTag()
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
	if (trueFlavor == KLFitter::Particles::kLight && isBTagged == true)
          probbtag = 0.;
      }
      
      if (probbtag > 0)
	logprob += log(probbtag); 
      else
	return -1e99; 
    }
    else if(fBTagMethod == kVetoLight){
      // loop over all model particles.  calculate the overall b-tagging
      // probability which is the product of all probabilities. 
      for (int i = 0; i < fParticlesModel->NPartons(); ++i){
        // get index of corresponding measured particle.                                                                   
        int index = fParticlesModel->JetIndex(i);
        if (index < 0)
          continue;

        KLFitter::Particles::TrueFlavorType trueFlavor = fParticlesModel->TrueFlavor(i);
        bool isBTagged = fParticlesModel->IsBTagged(i);
	if (trueFlavor == KLFitter::Particles::kB && isBTagged == false)
          probbtag = 0.;
      }
      
      if (probbtag > 0)
	logprob += log(probbtag); 
      else
	return -1e99; 
    }
    else if(fBTagMethod == kVetoBoth){
      // loop over all model particles.  calculate the overall b-tagging
      // probability which is the product of all probabilities. 
      for (int i = 0; i < fParticlesModel->NPartons(); ++i){
        // get index of corresponding measured particle.                                                                   
        int index = fParticlesModel->JetIndex(i);
        if (index < 0)
          continue;

        KLFitter::Particles::TrueFlavorType trueFlavor = fParticlesModel->TrueFlavor(i);
        bool isBTagged = fParticlesModel->IsBTagged(i);
	if (trueFlavor == KLFitter::Particles::kLight && isBTagged == true)
          probbtag = 0.;
	if (trueFlavor == KLFitter::Particles::kB && isBTagged == false)
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

	if(trueFlavor == KLFitter::Particles::kLight && isBTagged)
          logprob += log(1./rejection);
	else if(trueFlavor == KLFitter::Particles::kLight && !isBTagged)
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
void KLFitter::LikelihoodBase::PropagateBTaggingInformation()
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

      fParticlesModel->SetIsBTagged(         index, (*fParticlesPermuted)->IsBTagged(index));
      fParticlesModel->SetBTaggingEfficiency(index, (*fParticlesPermuted)->BTaggingEfficiency(index));
      fParticlesModel->SetBTaggingRejection( index, (*fParticlesPermuted)->BTaggingRejection(index));
    }
}
// ---------------------------------------------------------.

 std::vector <double> KLFitter::LikelihoodBase::GetBestFitParameters()
{
if (fCachedParameters.size() > 0) {
	return fCachedParameters;
} 
else return BCModel::GetBestFitParameters();

}



// ---------------------------------------------------------.

 std::vector <double> KLFitter::LikelihoodBase::GetBestFitParameterErrors()
{

if (fCachedParameterErrors.size() > 0) {
	return fCachedParameterErrors;
}
else return BCModel::GetBestFitParameterErrors();

}

// ---------------------------------------------------------

int KLFitter::LikelihoodBase::SetParametersToCache(int iperm, int nperms)
{
//set correct size of cachevector
if (iperm==0) {
	fCachedParametersVector.clear();
	fCachedParametersVector.assign(nperms, std::vector<double>(NParameters(), 0));

	fCachedParameterErrorsVector.clear();
	fCachedParameterErrorsVector.assign(nperms, std::vector<double>(NParameters(), 0));

	fCachedNormalizationVector.clear();
	fCachedNormalizationVector.assign(nperms, 0.);
}	

if ((iperm>(int)fCachedParametersVector.size())||(iperm>(int)fCachedParameterErrorsVector.size())) {
	std::cout<<"KLFitter::LikelihoodBase::SetParametersToCache: iperm > size of fCachedParametersVector or fCachedParameterErrorsVector!"<<std::endl;
	return 0;
}
fCachedParametersVector.at(iperm) = BCModel::GetBestFitParameters();
fCachedParameterErrorsVector.at(iperm) = BCModel::GetBestFitParameterErrors();
fCachedNormalizationVector.at(iperm) = BCModel::GetNormalization();

int switchpar1 = -1;
int switchpar2 = -1;
double switchcache = 0;
int partner = LHInvariantPermutationPartner(iperm, nperms, switchpar1, switchpar2);

if (partner > iperm) {

	if (((int)fCachedParametersVector.size() > partner)&&((int)fCachedParameterErrorsVector.size() > partner)){
		fCachedParametersVector.at(partner) = BCModel::GetBestFitParameters();
		switchcache = fCachedParametersVector.at(partner).at(switchpar1);
		fCachedParametersVector.at(partner).at(switchpar1) = fCachedParametersVector.at(partner).at(switchpar2);
		fCachedParametersVector.at(partner).at(switchpar2) = switchcache;

		fCachedParameterErrorsVector.at(partner) = BCModel::GetBestFitParameterErrors();
		switchcache = fCachedParameterErrorsVector.at(partner).at(switchpar1);
		fCachedParameterErrorsVector.at(partner).at(switchpar1) = fCachedParameterErrorsVector.at(partner).at(switchpar2);
		fCachedParameterErrorsVector.at(partner).at(switchpar2) = switchcache;

		fCachedNormalizationVector.at(partner) = BCModel::GetNormalization();


	} else {
		std::cout<<"KLFitter::LikelihoodBase::SetParametersToCache: size of fCachedParametersVector too small!"<<std::endl;
	}	
}
GetParametersFromCache(iperm);

return 1;
}
// ---------------------------------------------------------

int KLFitter::LikelihoodBase::GetParametersFromCache(int iperm)
{

	if (((int)fCachedParametersVector.size() > iperm)&&((int)fCachedParameterErrorsVector.size() > iperm)) {
		fCachedParameters = fCachedParametersVector.at(iperm);
		fCachedParameterErrors = fCachedParameterErrorsVector.at(iperm);
		fCachedNormalization = fCachedNormalizationVector.at(iperm);
	} else {
		std::cout<<"KLFitter::LikelihoodBase::GetParametersFromCache: size of fCachedParametersVector,  fCachedParameterErrorsVector or fCachedNormalizationVector too small!"<<std::endl;
	}
return 1;
}

// ---------------------------------------------------------.

 double KLFitter::LikelihoodBase::GetNormalization()
{

if (fCachedNormalizationVector.size() > 0) {
	return fCachedNormalization;
}
else {
return BCModel::GetNormalization();
}
}
// ---------------------------------------------------------.

 int KLFitter::LikelihoodBase::ResetCache()
{
fCachedParameters.clear();
fCachedParameterErrors.clear();

fCachedNormalization = 0.;

return 1;
}
// ---------------------------------------------------------.

 double KLFitter::LikelihoodBase::GetBestFitParameter(unsigned int index)
{
if (fCachedParameters.size() > 0) {
	return fCachedParameters.at(index);
}
else return BCModel::GetBestFitParameter(index);

}



// ---------------------------------------------------------.

 double KLFitter::LikelihoodBase::GetBestFitParameterError(unsigned int index)
{

if (fCachedParameterErrors.size() > 0) {
	return fCachedParameterErrors.at(index);
}
else return BCModel::GetBestFitParameterError(index);

}
