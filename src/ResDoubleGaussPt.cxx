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

#include "KLFitter/ResDoubleGaussPt.h" 
#include <iostream>
#include <cmath>

// --------------------------------------------------------- 
KLFitter::ResDoubleGaussPt::ResDoubleGaussPt(const char * filename) : KLFitter::ResDoubleGaussBase(filename)
{ }

// ---------------------------------------------------------
KLFitter::ResDoubleGaussPt::ResDoubleGaussPt(std::vector<double> const& parameters) : KLFitter::ResDoubleGaussBase(parameters)
{ }

// --------------------------------------------------------- 
KLFitter::ResDoubleGaussPt::~ResDoubleGaussPt()
{
  ;
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGaussPt::GetMean1(double x){
  return fParameters[0] + x * fParameters[1];
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGaussPt::GetSigma1(double x){
  return fParameters[2] + x * fParameters[3];
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGaussPt::GetAmplitude2(double x){
  return fParameters[4] + x * fParameters[5];
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGaussPt::GetMean2(double x){
  return fParameters[6] + x * fParameters[7];
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGaussPt::GetSigma2(double x){
  return fParameters[8] + x * fParameters[9];
}

// --------------------------------------------------------- 
