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

#include "KLFitter/Fitter.h"

#include <iostream>

#include "BAT/BCParameter.h"
#include "KLFitter/DetectorBase.h"
#include "KLFitter/LikelihoodBase.h"
#include "KLFitter/ParticleCollection.h"
#include "KLFitter/Permutations.h"

// ---------------------------------------------------------
KLFitter::Fitter::Fitter()
  : fDetector(nullptr)
  , fParticles(nullptr)
  , ETmiss_x(0.)
  , ETmiss_y(0.)
  , SumET(0.)
  , fParticlesPermuted(nullptr)
  , fMyParticlesTruth(nullptr)
  , fLikelihood(nullptr)
  , fPermutations(std::unique_ptr<KLFitter::Permutations>(new KLFitter::Permutations{&fParticles, &fParticlesPermuted}))
  , fMinuitStatus(0)
  , fConvergenceStatus(0)
  , fTurnOffSA(false)
  , fMinimizationMethod(kMinuit) {
  // empty
}

// ---------------------------------------------------------
KLFitter::Fitter::~Fitter() = default;

// ---------------------------------------------------------
int KLFitter::Fitter::SetParticles(KLFitter::ParticleCollection * particles, int nPartonsInPermutations) {
  fParticles = particles;

  // reset old table of permutations
  if (fPermutations)
    fPermutations->Reset();

  // create table of permutations
  fPermutations->CreatePermutations(nPartonsInPermutations);

  // remove invariant permutations if likelihood exists
  if (fLikelihood)
    fLikelihood->RemoveInvariantParticlePermutations();

  // remove forbidden permutations
  if (fLikelihood)
    fLikelihood->RemoveForbiddenParticlePermutations();

  // check if any permutations are left
  if (fPermutations->NPermutations() == 0) {
    std::cout << "KLFitter::Fitter::SetParticles(). No permutations left to fit. Are you vetoing?" << std::endl;
    return 0;
  }

  // set first permutation
  if (!fPermutations->SetPermutation(0))
    return 0;

  // get new permutation
  fParticlesPermuted = fPermutations->ParticlesPermuted();

  // no error
  return 1;
}

// ---------------------------------------------------------
int KLFitter::Fitter::SetMyParticlesTruth(KLFitter::ParticleCollection * particles) {
  fMyParticlesTruth = particles;

  // no error
  return 1;
}

// ---------------------------------------------------------
int KLFitter::Fitter::SetET_miss_XY_SumET(double etx, double ety, double sumet) {
  // set missing ET x and y component and sumET
  ETmiss_x = etx;
  ETmiss_y = ety;
  SumET = sumet;
  // no error
  return 1;
}

// ---------------------------------------------------------
int KLFitter::Fitter::SetDetector(KLFitter::DetectorBase * detector) {
  // set detector
  fDetector = detector;

  // If likelihood != nullptr, request the resolution functions now.
  if (fLikelihood) fLikelihood->RequestResolutionFunctions();

  // Return the status of the detector.
  return fDetector->Status();
}

// ---------------------------------------------------------
int KLFitter::Fitter::SetLikelihood(KLFitter::LikelihoodBase * likelihood) {
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

  // If detector != nullptr, return its status right now.
  if (fDetector) return fDetector->Status();

  // no error
  return 1;
}

// ---------------------------------------------------------
int KLFitter::Fitter::Fit(int index) {
  fLikelihood->ResetCache();
  fLikelihood->ResetResults();
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
  int partnerindex = fLikelihood->LHInvariantPermutationPartner(index, nperms, &dummy, &dummy);

  // check if permutation is LH invariant and has already been calculated
  if ((partnerindex > -1)&&(partnerindex < index)) {
    fLikelihood->GetParametersFromCache(index);
    GetFitStatusFromCache(index);
  } else {
    // Markov Chain MC
    if (fMinimizationMethod == kMarkovChainMC) {
      fLikelihood->MCMCSetFlagFillHistograms(true);
      fLikelihood->MCMCSetNChains(5);
      fLikelihood->MCMCSetNIterationsRun(20000);
      fLikelihood->MCMCSetNIterationsMax(1000000);
      fLikelihood->MCMCSetNIterationsUpdate(100);
      fLikelihood->MarginalizeAll();
    } else if (fMinimizationMethod == kSimulatedAnnealing) {
      // simulated annealing
      fLikelihood->SetOptimizationMethod(BCIntegrate::kOptSimAnn);
      fLikelihood->SetSAT0(10);
      fLikelihood->SetSATmin(0.001);
      fLikelihood->FindMode(fLikelihood->GetInitialParameters());
    } else if (fMinimizationMethod == kMinuit) {
      // MINUIT
      fLikelihood->SetOptimizationMethod(BCIntegrate::kOptMinuit);
      fLikelihood->FindMode(fLikelihood->GetInitialParameters());

      fMinuitStatus = fLikelihood->GetMinuitErrorFlag();

      // check if any parameter is at its borders->set MINUIT flag to 500
      if (fMinuitStatus == 0) {
        std::vector<double> BestParameters = fLikelihood->GetBestFitParameters();
        for (unsigned int iPar = 0; iPar < fLikelihood->GetNParameters(); iPar++) {
          if (fLikelihood->GetParameter(0)->IsAtLimit(BestParameters[iPar])) {
            fMinuitStatus = 500;
          }
        }
      }
      if (fLikelihood->GetFlagIsNan()== true) {
        fMinuitStatus = 508;
      }

      // re-run if Minuit status bad
      if (fMinuitStatus != 0) {
        fLikelihood->ResetCache();
        fLikelihood->ResetResults();
        if (!fTurnOffSA) {
          fLikelihood->SetFlagIsNan(false);
          fLikelihood->SetOptimizationMethod(BCIntegrate::kOptSimAnn);
          fLikelihood->FindMode(fLikelihood->GetInitialParameters());
        }

        fLikelihood->SetOptimizationMethod(BCIntegrate::kOptMinuit);
        fLikelihood->FindMode(fLikelihood->GetBestFitParameters());
        fMinuitStatus = fLikelihood->GetMinuitErrorFlag();
      }

      fConvergenceStatus = 0;
      if (fMinuitStatus == 4)
        fConvergenceStatus |= MinuitDidNotConvergeMask;
    }

    // check if any parameter is at its borders->set MINUIT flag to 501
    if (fMinuitStatus == 0) {
      std::vector<double> BestParameters = fLikelihood->GetBestFitParameters();
      for (unsigned int iPar = 0; iPar < fLikelihood->GetNParameters(); iPar++) {
        if (fLikelihood->GetParameter(0)->IsAtLimit(BestParameters[iPar])) {
          fMinuitStatus = 501;
          fConvergenceStatus |= AtLeastOneFitParameterAtItsLimitMask;
        }
      }
    }
    if (fLikelihood->GetFlagIsNan()== true) {
      fMinuitStatus = 509;
      fConvergenceStatus |= FitAbortedDueToNaNMask;
    } else {
      // check if TF problem
      if (!fLikelihood->NoTFProblem(fLikelihood->GetBestFitParameters())) {
        fMinuitStatus = 510;
        fConvergenceStatus |= InvalidTransferFunctionAtConvergenceMask;
      }
    }

    // calculate integral
    if (fLikelihood->FlagIntegrate()) {
      fLikelihood->SetIntegrationMethod(BCIntegrate::kIntCuba);
      fLikelihood->Normalize();
    }

    // caching parameters
    fLikelihood->SetParametersToCache(index, nperms);
    SetFitStatusToCache(index, nperms);
  }  // end of fitting "else"

  // no error
  return 1;
}

// ---------------------------------------------------------
int KLFitter::Fitter::Fit() {
  // check status
  if (!Status())
    return 0;

  // get number of permutations
  int npermutations = fPermutations->NPermutations();

  // loop over all permutations
  for (int ipermutation = 0; ipermutation < npermutations; ++ipermutation) {
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
    fLikelihood->FindMode(BCIntegrate::kOptMinuit, fLikelihood->GetBestFitParameters());
    fMinuitStatus = fLikelihood->GetMinuitErrorFlag();
  }

  // no error
  return 1;
}

// ---------------------------------------------------------
int KLFitter::Fitter::Status() {
  // check if measured particles exist
  if (!fParticles) {
    std::cout << "KLFitter::Fitter::Status(). Set of measured particles not defined." << std::endl;
    return 0;
  }

  // check if detector exists
  if (!fDetector) {
    std::cout << "KLFitter::Fitter::Status(). No detector defined." << std::endl;
    return 0;
  }

  // check detector
  if (!fDetector->Status()) {
    return 0;
  }

  // no error
  return 1;
}

// ---------------------------------------------------------
int KLFitter::Fitter::GetFitStatusFromCache(int iperm) {
  if ((static_cast<int>(fCachedConvergenceStatusVector.size()) > iperm) && (static_cast<int>(fCachedMinuitStatusVector.size()) > iperm)) {
    fConvergenceStatus = fCachedConvergenceStatusVector.at(iperm);
    fMinuitStatus = fCachedMinuitStatusVector.at(iperm);
  } else {
    std::cout << "KLFitter::Fitter::GetFitStatusFromCache: size of fCachedConvergenceStatusVector or fCachedMinuitStatusVector too small!" << std::endl;
  }

  return 1;
}

// ---------------------------------------------------------
int KLFitter::Fitter::SetFitStatusToCache(int iperm, int nperms) {
  if (iperm == 0) {
    fCachedMinuitStatusVector.clear();
    fCachedMinuitStatusVector.assign(nperms, -1);

    fCachedConvergenceStatusVector.clear();
    fCachedConvergenceStatusVector.assign(nperms, -1);
  }

  if ((iperm > static_cast<int>(fCachedMinuitStatusVector.size())) || (iperm > static_cast<int>(fCachedConvergenceStatusVector.size()))) {
    std::cout << "KLFitter::Fitter::SetFitStatusToCache: iperm > size of fCachedMinuitStatusVector or fCachedConvergenceStatusVector!" << std::endl;
    return 0;
  }
  fCachedMinuitStatusVector.at(iperm) = fMinuitStatus;
  fCachedConvergenceStatusVector.at(iperm) = fConvergenceStatus;

  int dummy;
  int partner = fLikelihood->LHInvariantPermutationPartner(iperm, nperms, &dummy, &dummy);

  if (partner > iperm) {
    if ((static_cast<int>(fCachedMinuitStatusVector.size()) > partner) && (static_cast<int>(fCachedConvergenceStatusVector.size()) > partner)) {
      fCachedMinuitStatusVector.at(partner) = fMinuitStatus;
      fCachedConvergenceStatusVector.at(partner) = fConvergenceStatus;
    } else {
      std::cout << "KLFitter::Fitter::SetFitStatusToCache: size of fCachedMinuitStatusVector or fCachedConvergenceStatusVector too small!" << std::endl;
    }
  }

  GetFitStatusFromCache(iperm);
  return 1;
}

// ---------------------------------------------------------.
int KLFitter::Fitter::ResetCache() {
  fMinuitStatus = -1;
  fConvergenceStatus = -1;

  return 1;
}
