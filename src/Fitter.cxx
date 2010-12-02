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
  fParticlesPermuted = 0; 
  fPermutations = new KLFitter::Permutations(&fParticles, &fParticlesPermuted);         
  fMinuitStatus = 0; 
  fTurnOffSA = false;
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

  // set first permutation 
  if (!fPermutations->SetPermutation(0))
    return 0; 

  // get new permutation 
  fParticlesPermuted = fPermutations->ParticlesPermuted(); 

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::Fitter::SetET_miss_XY(double etx, double ety)
{
  // set missing ET x and y component
  ETmiss_x = etx;
  ETmiss_y = ety;

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
  // check status
  if (!Status())
    return 0; 

  // set permutation 
  if (!fPermutations->SetPermutation(index))
    return 0; 

  // get new permutation 
  fParticlesPermuted = fPermutations->ParticlesPermuted(); 

  // set missing ET x and y components
  fLikelihood->SetET_miss_XY(ETmiss_x, ETmiss_y);

  // initialize likelihood (likelihood MUST be initialized after
  // setting the missing ET, because AdjustParameterRanges() might
  // make use of the missing ET information !!!)
  fLikelihood->Initialize();
  fLikelihood->SetFlagIsNan(false);

  // set flavor tags if b-tagging is on
  if (fLikelihood->FlagBTagging()) {
    fLikelihood->CalculateFlavorTags(); // not to be written with "ou" 
  }

  // perform fitting 
  //    fLikelihood->MCMCSetFlagFillHistograms(false); 
  //    fLikelihood->MCMCSetNChains(5); 
  //    fLikelihood->MCMCSetNIterationsRun(2000); 
  //    fLikelihood->MCMCSetNIterationsMax(1000); 
  //    fLikelihood->MCMCSetNIterationsUpdate(100); 
  //    fLikelihood->MarginalizeAll();

  // simulated annealing
  //    fLikelihood->SetOptimizationMethod( BCIntegrate::kOptSA );
  //    fLikelihood->SetSAT0(10);
  //    fLikelihood->SetSATmin(0.001);
        
  //    fLikelihood->FindMode( fLikelihood->GetInitialParameters() );

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
          if ( ParameterSet->at(iPar)->IsAtLimit(BestParameters.at(iPar)) )
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
        
  // check if any parameter is at its borders->set MINUIT flag to 501
  if ( fMinuitStatus == 0)
    {
      std::vector<double> BestParameters = fLikelihood->GetBestFitParameters();
      BCParameterSet * ParameterSet = fLikelihood->GetParameterSet();
      for (unsigned int iPar = 0; iPar < BestParameters.size(); iPar++)
        {
          if ( ParameterSet->at(iPar)->IsAtLimit(BestParameters.at(iPar)) )
            {
              fMinuitStatus = 501;
            }
        }
    }
  if(fLikelihood->GetFlagIsNan()==true)
    {
      fMinuitStatus=509;
    }
  else {
    // check if TF problem
    if (! fLikelihood->NoTFProblem(fLikelihood->GetBestFitParameters()))
      fMinuitStatus = 510;
  }

  // check b-tagging and calculate probability
  if (fLikelihood->FlagBTagging()) // && fLikelihood->FlagIntegrate())
    if (fLikelihood->BTaggingProbability() == 0)
      return 1; 

  // calculate integral 
  if (fLikelihood->FlagIntegrate())
    {
      fLikelihood->SetIntegrationMethod(BCIntegrate::kIntCuba); 
      fLikelihood->Normalize();
    }

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

