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

#include "KLFitter/ResolutionBase.h"
#include <iostream>
#include <fstream>

// ---------------------------------------------------------
KLFitter::ResolutionBase::ResolutionBase(int npar) {
  if (npar < 0)
    npar = 0;

  fNParameters = npar;

  fParameters.clear();
  fParameters.assign(npar, 0);
}

// ---------------------------------------------------------
KLFitter::ResolutionBase::ResolutionBase(std::vector <double> parameters) {
  fNParameters = int(parameters.size());

  // clear parameters
  fParameters.clear();

  // copy values
  for (int i = 0; i < fNParameters; ++i)
    fParameters.push_back(parameters[i]);
}

// ---------------------------------------------------------
KLFitter::ResolutionBase::~ResolutionBase() {
  fParameters.clear();
}

// ---------------------------------------------------------
int KLFitter::ResolutionBase::Par(int index, double &par) {
  // check parameter range
  if (index < 0 || index >= fNParameters) {
    std::cout << "KLFitter:ResolutionBase::Par(). Index out of range." << std::endl;
    // error
    return 0;
  }

  par = fParameters[index];

  // no error
  return 1;
}

// ---------------------------------------------------------
int KLFitter::ResolutionBase::SetPar(int index, double value) {
  // check parameter range
  if (index < 0 || index >= fNParameters) {
    std::cout << "KLFitter:ResolutionBase::SetPar(). Index out of range." << std::endl;
    // error
    return 1;
  }

  fParameters[index] = value;

  // no error
  return 1;
}

// ---------------------------------------------------------
int KLFitter::ResolutionBase::SetPar(std::vector <double> parameters) {
  // check vector size
  int npar = parameters.size();

  if (npar != fNParameters) {
    std::cout << "KLFitter:ResolutionBase::SetPar(). Number of parameters is inconsistent." << std::endl;

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
int KLFitter::ResolutionBase::ReadParameters(const char * filename, int nparameters) {
  // define input file
  std::ifstream inputfile;

  // open file
  inputfile.open(filename);

  // check if file is open
  if (!inputfile.is_open()) {
    std::cout << "KLFitter::ResolutionBase::ReadParameters(). File \"" << filename << "\" not found." << std::endl;
    return 0;
  }

  // reset parameters
  fParameters.clear();

  // read parameters
  for (int i = 0; i < nparameters; ++i) {
    double par = 0.0;
    inputfile >> par;
    fParameters.push_back(par);
  }

  // close file
  inputfile.close();

  // no error
  return 1;
}

// ---------------------------------------------------------
