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
 * along with KLFitter. If not, see <http:// www.gnu.org/licenses/>.
 */

#include "KLFitter/Resolutions/ResolutionBase.h"

#include <fstream>
#include <iostream>

// ---------------------------------------------------------
KLFitter::Resolutions::ResolutionBase::ResolutionBase(int npar) {
  if (npar < 0)
    npar = 0;

  fNParameters = npar;

  fParameters.clear();
  fParameters.assign(npar, 0);
}

// ---------------------------------------------------------
KLFitter::Resolutions::ResolutionBase::ResolutionBase(std::vector <double> parameters) {
  fNParameters = parameters.size();

  // clear parameters
  fParameters.clear();

  // copy values
  for (int i = 0; i < fNParameters; ++i)
    fParameters.push_back(parameters[i]);
}

// ---------------------------------------------------------
KLFitter::Resolutions::ResolutionBase::~ResolutionBase() = default;

// ---------------------------------------------------------
int KLFitter::Resolutions::ResolutionBase::Par(int index, double *par) {
  // check parameter range
  if (index < 0 || index >= fNParameters) {
    std::cout << "KLFitter::Resolutions::ResolutionBase::Par(). Index out of range." << std::endl;
    // error
    return 0;
  }

  *par = fParameters[index];

  // no error
  return 1;
}

// ---------------------------------------------------------
int KLFitter::Resolutions::ResolutionBase::SetPar(int index, double value) {
  // check parameter range
  if (index < 0 || index >= fNParameters) {
    std::cout << "KLFitter::Resolutions::ResolutionBase::SetPar(). Index out of range." << std::endl;
    // error
    return 1;
  }

  fParameters[index] = value;

  // no error
  return 1;
}

// ---------------------------------------------------------
int KLFitter::Resolutions::ResolutionBase::SetPar(std::vector <double> parameters) {
  // check vector size
  int npar = parameters.size();

  if (npar != fNParameters) {
    std::cout << "KLFitter::Resolutions::ResolutionBase::SetPar(). Number of parameters is inconsistent." << std::endl;

    // return error code
    return 1;
  }

  // set parameters
  for (int i = 0; i < fNParameters; ++i)
    fParameters[i] = parameters[i];

  // no error
  return 1;
}

// ---------------------------------------------------------
int KLFitter::Resolutions::ResolutionBase::ReadParameters(const char * filename, std::size_t nparameters) {
  // define input file
  std::ifstream inputfile;

  // open file
  inputfile.open(filename);

  // check if file is open
  if (!inputfile.is_open()) {
    fStatus = 0;
    return 0;
  }

  // reset parameters
  fParameters.clear();

  // read parameters
  double par = 0.0;
  while (inputfile >> par) {
    fParameters.push_back(par);
  }

  if (fParameters.size() != nparameters){
    std::cout << "KLFitter::Resolutions::ResolutionBase::ReadParameters(). Expecting " << nparameters
      << ", parameters for Transfer functions but " << fParameters.size() << " parameters found" << std::endl;
    return 0;
  }

  // close file
  inputfile.close();

  // Switch internal status to "ok".
  fStatus = 1;

  // no error
  return 1;
}
