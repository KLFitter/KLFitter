#include "Fitter.h"
#include "DetectorBase.h" 
#include "LikelihoodBase.h" 
#include "Permutations.h" 
#include "Particles.h" 

#include <iostream> 

// --------------------------------------------------------- 
KLFitter::Fitter::Fitter()
{
  fDetector = 0; 
  fLikelihood = 0; 
  fParticles = 0; 
  ETmiss_x = 0.;
  ETmiss_y = 0.;
  SumET = 0.;
  fParticlesPermuted = 0; 
  fMyParticlesTruth = 0;
  fPermutations = new KLFitter::Permutations(&fParticles, &fParticlesPermuted);         
  fMinuitStatus = 0; 
  fConvergenceStatus = 0;
  fTurnOffSA = false;
  fMinimizationMethod = kMinuit;
}

// --------------------------------------------------------- 
KLFitter::Fitter::~Fitter()
{
  if (fPermutations)
    delete fPermutations; 
}

// --------------------------------------------------------- 
int KLFitter::Fitter::SetParticles(KLFitter::Particles * particles)
{
  fParticles = particles; 

  // reset old table of permutations
  if (fPermutations)
    fPermutations->Reset(); 

  // create table of permutations 
  fPermutations->CreatePermutations();  

  // remove invariant permutations if likelihood exists 
  if (fLikelihood)
    fLikelihood->RemoveInvariantParticlePermutations();
  
  // remove forbidden permutations
  if (fLikelihood)
    fLikelihood->RemoveForbiddenParticlePermutations(); 
  
  // set first permutation 
  if (!fPermutations->SetPermutation(0))
    return 0; 

  // get new permutation 
  fParticlesPermuted = fPermutations->ParticlesPermuted(); 

  // no error 
  return 1; 
}
// --------------------------------------------------------- 
int KLFitter::Fitter::SetMyParticlesTruth(KLFitter::Particles * particles)
{
  fMyParticlesTruth = particles;

  // set pointer to truth particles 
  
  fLikelihood->SetMyParticlesTruth(&fMyParticlesTruth); 
  
  //std:: cout << "fMyParticlesTruth set in Fitter!!!!!" << std::endl;

  // no error
  return 1;
}
// --------------------------------------------------------- 
int KLFitter::Fitter::SetET_miss_XY_SumET(double etx, double ety, double sumet)
{
  // set missing ET x and y component and sumET
  ETmiss_x = etx;
  ETmiss_y = ety;
  SumET = sumet;
  // no error
  return 1;
}

// --------------------------------------------------------- 
int KLFitter::Fitter::SetDetector(KLFitter::DetectorBase * detector)
{
  // set detector 
  fDetector = detector; 

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::Fitter::SetLikelihood(KLFitter::LikelihoodBase * likelihood)
{
  // set likelihood 
  fLikelihood = likelihood; 
        
  // set pointer to pointer of detector 
  fLikelihood->SetDetector(&fDetector); 

  // set pointer to pointer of permutation object 
  fLikelihood->SetPermutations(&fPermutations); 

  // set pointer to permuted particles 
  fLikelihood->SetParticlesPermuted(&fParticlesPermuted); 

  
  
  // remove invariant permutations if particles are defined alreday
  if (fParticles)
    fLikelihood->RemoveInvariantParticlePermutations(); 

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::Fitter::Fit(int index)
{
  fLikelihood->ResetCache();
  ResetCache();

  // check status
  if (!Status())
    return 0; 

  // set permutation 
  if (!fPermutations->SetPermutation(index))
    return 0; 

  // get new permutation 
  fParticlesPermuted = fPermutations->ParticlesPermuted(); 

  // set missing ET x and y components and the sumET
  fLikelihood->SetET_miss_XY_SumET(ETmiss_x, ETmiss_y, SumET);

  // initialize likelihood (likelihood MUST be initialized after
  // setting the missing ET, because AdjustParameterRanges() might
  // make use of the missing ET information !!!)
  fLikelihood->Initialize();
  fLikelihood->SetFlagIsNan(false);

  // set flavor tags if b-tagging is on
  if (fLikelihood->GetBTagging() != LikelihoodBase::kNotag) {
    fLikelihood->PropagateBTaggingInformation();
  }

  // perform fitting 
  // Check if LH is invariant
  int dummy;
  int nperms = fPermutations->NPermutations();
  int partnerindex = fLikelihood->LHInvariantPermutationPartner(index, nperms, dummy, dummy);


  // check if permutation is LH invariant and has already been calculated
  if ((partnerindex > -1)&&(partnerindex < index)){
	fLikelihood->GetParametersFromCache(index);
	GetFitStatusFromCache(index);

  } else {

    // Markov Chain MC
    if (fMinimizationMethod == kMarkovChainMC) {
      // comment out BCLog lines if interested in more BAT details
      //BCLog * log = new BCLog();
      //log->OpenLog("help.txt",BCLog::detail,BCLog::detail);
      fLikelihood->MCMCSetFlagFillHistograms(true); 
      fLikelihood->MCMCSetNChains(5); 
      fLikelihood->MCMCSetNIterationsRun(20000); 
      fLikelihood->MCMCSetNIterationsMax(1000000); 
      fLikelihood->MCMCSetNIterationsUpdate(100); 
      fLikelihood->MarginalizeAll();
      //log->CloseLog();
      //delete log;
    }
    // simulated annealing
    else if (fMinimizationMethod == kSimulatedAnnealing) {
      fLikelihood->SetOptimizationMethod( BCIntegrate::kOptSA );
      fLikelihood->SetSAT0(10);
      fLikelihood->SetSATmin(0.001);
      fLikelihood->FindMode( fLikelihood->GetInitialParameters() );
    }
    // MINUIT
    else if (fMinimizationMethod == kMinuit) {
      fLikelihood->SetOptimizationMethod( BCIntegrate::kOptMinuit); 
      fLikelihood->FindMode( fLikelihood->GetInitialParameters() );
      //    fLikelihood->FindMode( fLikelihood->GetBestFitParameters() ); 

      fMinuitStatus = fLikelihood->GetMinuitErrorFlag(); 

      // check if any parameter is at its borders->set MINUIT flag to 500
      if ( fMinuitStatus == 0)
        {
          std::vector<double> BestParameters = fLikelihood->GetBestFitParameters();
          BCParameterSet * ParameterSet = fLikelihood->GetParameterSet();
          for (unsigned int iPar = 0; iPar < BestParameters.size(); iPar++)
            {
              if ( (*ParameterSet)[iPar]->IsAtLimit(BestParameters[iPar]) )
                {
                  fMinuitStatus = 500;
                }
            }
        }
      if(fLikelihood->GetFlagIsNan()==true)
        {
          fMinuitStatus=508;
        }

      // re-run if Minuit status bad 
      if (fMinuitStatus != 0)
        {
	  fLikelihood->ResetCache();
          // print to screen
          //        std::cout << "KLFitter::Fit(). Minuit did not find proper minimum. Rerun with Simulated Annealing."<<std::endl; 
          if (!fTurnOffSA) {
            fLikelihood->SetFlagIsNan(false);
            fLikelihood->SetOptimizationMethod( BCIntegrate::kOptSA );
            fLikelihood->FindMode( fLikelihood->GetInitialParameters() );
          }
                          
          fLikelihood->SetOptimizationMethod( BCIntegrate::kOptMinuit); 
          fLikelihood->FindMode( fLikelihood->GetBestFitParameters() ); 
          fMinuitStatus = fLikelihood->GetMinuitErrorFlag(); 
        }   

      fConvergenceStatus = 0;
      if (fMinuitStatus == 4)
        fConvergenceStatus |= MinuitDidNotConvergeMask;
    }

          
    // check if any parameter is at its borders->set MINUIT flag to 501
    if ( fMinuitStatus == 0)
      {
        std::vector<double> BestParameters = fLikelihood->GetBestFitParameters();
        BCParameterSet * ParameterSet = fLikelihood->GetParameterSet();
        for (unsigned int iPar = 0; iPar < BestParameters.size(); iPar++)
          {
            if ( (*ParameterSet)[iPar]->IsAtLimit(BestParameters[iPar]) )
              {
                fMinuitStatus = 501;
                fConvergenceStatus |= AtLeastOneFitParameterAtItsLimitMask;
              }
          }
      }
    if(fLikelihood->GetFlagIsNan()==true)
      {
        fMinuitStatus=509;
        fConvergenceStatus |= FitAbortedDueToNaNMask;
      }
    else {
      // check if TF problem
      if (! fLikelihood->NoTFProblem(fLikelihood->GetBestFitParameters())) {
        fMinuitStatus = 510;
        fConvergenceStatus |= InvalidTransferFunctionAtConvergenceMask;
    }
    }
  
    // calculate integral 
    if (fLikelihood->FlagIntegrate())
      {
        fLikelihood->SetIntegrationMethod(BCIntegrate::kIntCuba); 
      fLikelihood->Normalize();
      }

    // caching parameters
    fLikelihood->SetParametersToCache(index, nperms);
    SetFitStatusToCache(index, nperms);

  } //end of fitting "else"

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::Fitter::Fit()
{
  // check status
  if (!Status())
    return 0; 

  // get number of permutations 
  int npermutations = fPermutations->NPermutations(); 

  // loop over all permutations
  for (int ipermutation = 0; ipermutation < npermutations; ++ipermutation)
    {
      // set permutation 
      if (!fPermutations->SetPermutation(ipermutation))
        return 0; 
                        
      // get new permutation 
      fParticlesPermuted = fPermutations->ParticlesPermuted(); 
                        
      // initialize likelihood 
      fLikelihood->Initialize(); 

      // perform fitting 
      fLikelihood->MCMCSetNChains(5); 
      fLikelihood->MCMCSetNIterationsRun(2000);
      fLikelihood->MCMCSetNIterationsMax(1000);
      fLikelihood->MCMCSetNIterationsUpdate(100);                       
      fLikelihood->MarginalizeAll();
      fLikelihood->FindModeMinuit(fLikelihood->GetBestFitParameters(), -1); 
      fMinuitStatus = fLikelihood->GetMinuitErrorFlag(); 
    }

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::Fitter::Status()
{
  // check if measured particles exist 
  if (!fParticles)
    {
      std::cout << "KLFitter::Fitter::Status(). Set of measured particles not defined." << std::endl;
      return 0; 
    } 

  // check if detector exists 
  if (!fDetector)
    {
      std::cout << "KLFitter::Fitter::Status(). No detector defined." << std::endl; 
      return 0; 
    }

  // check detector 
  if (!fDetector->Status())
    { 
      return 0; 
    }

  // no error
  return 1; 
}

// --------------------------------------------------------- 

int KLFitter::Fitter::GetFitStatusFromCache(int iperm)
{
	if (((int)fCachedConvergenceStatusVector.size() > iperm)&&((int)fCachedMinuitStatusVector.size() > iperm)) {
		fConvergenceStatus = fCachedConvergenceStatusVector.at(iperm);
		fMinuitStatus = fCachedMinuitStatusVector.at(iperm);
	} else {
		std::cout<<"KLFitter::Fitter::GetFitStatusFromCache: size of fCachedConvergenceStatusVector or fCachedMinuitStatusVector too small!"<<std::endl;
	}

return 1;
}

// ---------------------------------------------------------

int KLFitter::Fitter::SetFitStatusToCache(int iperm, int nperms)
{

if (iperm==0) {
	fCachedMinuitStatusVector.clear();
	fCachedMinuitStatusVector.assign(nperms, -1);

	fCachedConvergenceStatusVector.clear();
	fCachedConvergenceStatusVector.assign(nperms, -1);
}	

if ((iperm>(int)fCachedMinuitStatusVector.size())||(iperm>(int)fCachedConvergenceStatusVector.size())) {
	std::cout<<"KLFitter::Fitter::SetFitStatusToCache: iperm > size of fCachedMinuitStatusVector or fCachedConvergenceStatusVector!"<<std::endl;
	return 0;
}
fCachedMinuitStatusVector.at(iperm) = fMinuitStatus;
fCachedConvergenceStatusVector.at(iperm) = fConvergenceStatus;

int dummy;
int partner = fLikelihood->LHInvariantPermutationPartner(iperm, nperms, dummy, dummy);

if (partner > iperm) {

	if (((int)fCachedMinuitStatusVector.size() > partner)&&((int)fCachedConvergenceStatusVector.size() > partner)){

		fCachedMinuitStatusVector.at(partner) = fMinuitStatus;
		fCachedConvergenceStatusVector.at(partner) = fConvergenceStatus;

	} else {
		std::cout<<"KLFitter::Fitter::SetFitStatusToCache: size of fCachedMinuitStatusVector or fCachedConvergenceStatusVector too small!"<<std::endl;
	}	
}	

GetFitStatusFromCache(iperm);
return 1;
}

// ---------------------------------------------------------.

 int KLFitter::Fitter::ResetCache()
{
fMinuitStatus = -1;
fConvergenceStatus = -1;

return 1;
}
