#include "LikelihoodBase.h" 
#include "BAT/BCLog.h" 
#include <iostream> 
#include <string> 

// --------------------------------------------------------- 
KLFitter::LikelihoodBase::LikelihoodBase() : BCModel()
{
	fParticlesPermuted = 0; 
	fParticlesModel = 0; 
	fDetector = 0; 
	fLorentzVectorContainer = 0; 
	fFlagBTagging = 0; 
	fFlagIntegrate = 0; 
	fPhysicsConstants = new KLFitter::PhysicsConstants(); 
	fCutBTag = 1e4;
	BCLog::SetLogLevel(BCLog::nothing);
}

// --------------------------------------------------------- 
KLFitter::LikelihoodBase::LikelihoodBase(Particles ** particles) : BCModel() 
{
	fParticlesPermuted = particles; 
	fParticlesModel = 0; 
	fDetector = 0; 
	fFlagBTagging = 0; 
	fFlagIntegrate = 0; 
	fPhysicsConstants = new KLFitter::PhysicsConstants(); 
	fCutBTag = 1e4;
	BCLog::SetLogLevel(BCLog::nothing);}

// --------------------------------------------------------- 
KLFitter::LikelihoodBase::~LikelihoodBase()
{
 	if (fPhysicsConstants)
		delete fPhysicsConstants; 
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodBase::SetPhysicsConstants(KLFitter::PhysicsConstants * physicsconstants)
{
	fPhysicsConstants = physicsconstants; 

	// no error 
	return 1; 
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodBase::SetInitialParameters(std::vector <double> parameters)
{
	// check number of parameters 
	if (int(parameters.size()) != this -> NParameters())
		{
			std::cout << "KLFitter::SetInitialPosition(). Length of vector does not equal the number of parameters." << std::endl; 
			return 0; 
		}

	// set starting point for MCMC 
	this -> MCMCSetInitialPositions(parameters); 		

	// no error 
	return 1; 
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodBase::SetParameterRange(int index, double parmin, double parmax)
{
	// check index
	if (index < 0 || index >= this -> NParameters())
		{
			std::cout << " KLFitter::Combinatorics::SetParameterRange(). Index out of range." << std::endl; 
			return 0; 
		}

	// set parameter ranges in BAT 
	this -> GetParameter(index) -> SetLowerLimit(parmin); 
	this -> GetParameter(index) -> SetUpperLimit(parmax); 
	fMCMCBoundaryMin[index] = parmin; 
	fMCMCBoundaryMax[index] = parmax; 

	// no error 
	return 1; 
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodBase::SetDetector(KLFitter::DetectorBase ** detector)
{
	// set pointer to pointer of detector 
	fDetector = detector; 

	// no error 
	return 1; 
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodBase::SetParticlesPermuted(KLFitter::Particles ** particles)
{
	// set pointer to pointer of permuted particles
	fParticlesPermuted  = particles; 

	// no error
	return 1;
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodBase::SetPermutations(KLFitter::Permutations ** permutations)
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
	if (index < 0 || index >= this -> NParameters())
		{
			std::cout << " KLFitter::Combinatorics::ParMin(). Index out of range." << std::endl; 
			return 0; 
		}

	// return parameter range from BAT
	return this -> GetParameter(index) -> GetLowerLimit(); 
}

// --------------------------------------------------------- 
double KLFitter::LikelihoodBase::ParMax(int index)
{
	// check index
	if (index < 0 || index >= this -> NParameters())
		{
			std::cout << " KLFitter::Combinatorics::ParMax(). Index out of range." << std::endl; 
			return 0; 
		}

	// return parameter range from BAT
	return this -> GetParameter(index) -> GetUpperLimit(); 
}

// --------------------------------------------------------- 
double KLFitter::LikelihoodBase::LogEventProbability()
{
	double logprob = 0; 

	if (fFlagBTagging)
		{
			//			double probbtag = this -> BTaggingProbability();
			double probbtag = 1; 
			
			// loop over all model particles. 
			for (int i = 0; i < fParticlesModel -> NPartons(); ++i)
				if (fParticlesModel->BTaggingProbability(i)==0 &&
						fParticlesModel->FlavorTag(i)==1)
					probbtag=0;
			
			if (probbtag > 0)
				logprob += log(probbtag); 
			else
				return -1e99; 
		}

	if (fFlagIntegrate)
		{
			logprob += log(this -> GetNormalization()); 
		}
	else
		{
			logprob += this -> LogLikelihood( this -> GetBestFitParameters() ); 
		}

	return logprob; 
}

// --------------------------------------------------------- 
double KLFitter::LikelihoodBase::BTaggingProbability()
{
	double prob = 1.0; 

	// loop over all measured particles. important: assumes that index
	// of measured particle is the same as the index of the model
	// particle
	for (int i = 0; i < (*fParticlesPermuted) -> NPartons(); ++i)
		{
			if (fParticlesModel -> BTaggingProbability(i) > 0.5)
				prob *= (*fParticlesPermuted) -> BTaggingProbability(i); 
			else
				prob *= (1.0 - (*fParticlesPermuted) -> BTaggingProbability(i)); 
		}

	// return log of probability 
	return prob; 
}

// --------------------------------------------------------- 
void KLFitter::LikelihoodBase::CalculateFlavorTags()
{
	// loop over all model particles. 
	for (int i = 0; i < fParticlesModel -> NPartons(); ++i)
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
