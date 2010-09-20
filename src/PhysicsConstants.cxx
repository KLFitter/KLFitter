#include "PhysicsConstants.h" 
#include <cmath>
#include <iostream> 

// --------------------------------------------------------- 
KLFitter::PhysicsConstants::PhysicsConstants()
{
  fMassBottom =   4.7; // bottom quark mass in GeV/c^{2}
  fMassW      =  80.4; // W boson mass in GeV/c^{2}
  fMassTop    = 170.0; // top quark pole mass in GeV/c^{2}
  fGammaW     =   2.1; // W boson width 
  fGammaTop   =   1.5; // top quark width
  fGF         =   1.16637e-5; // in GeV^{-2} 
  fAlphaS     =   0.118; 
  fPi         =   3.14159265358979312; 
}

// --------------------------------------------------------- 
KLFitter::PhysicsConstants::~PhysicsConstants()
{
}

// --------------------------------------------------------- 
int KLFitter::PhysicsConstants::SetMassBottom(double mass)
{
  // check argument 
  if (mass < 0)
    {
      std::cout << "KLFitter::PhysicsConstants::SetMassBottom(). Mass cannot be negative. Set mass to zero." << std::endl; 
      fMassBottom = 0.; 
      return 0; 
    }

  //set mass 
  fMassBottom = mass; 

  // no error 
  return 1;
}

// --------------------------------------------------------- 
int KLFitter::PhysicsConstants::SetMassTop(double mass)
{
  // check argument 
  if (mass < 0)
    {
      std::cout << "KLFitter::PhysicsConstants::SetMassTop(). Mass cannot be negative. Set mass to zero." << std::endl; 
      fMassTop = 0.; 
      return 0; 
    }

  //set mass 
  fMassTop = mass; 

  // calculate top width
  CalculateGammaTop(); 

  // no error 
  return 1;
}

// --------------------------------------------------------- 
int KLFitter::PhysicsConstants::SetMassW(double mass)
{
  // check argument 
  if (mass < 0)
    {
      std::cout << "KLFitter::PhysicsConstants::SetMassW(). Mass cannot be negative. Set mass to zero." << std::endl; 
      fMassW = 0.; 
      return 0; 
    }

  //set mass 
  fMassW = mass; 

  // no error 
  return 1;
}

// --------------------------------------------------------- 
int KLFitter::PhysicsConstants::SetGammaW(double gamma)
{
  // check argument 
  if (gamma < 0)
    {
      std::cout << "KLFitter::PhysicsConstants::SetGammaW(). Width cannot be negative. Set width to zero." << std::endl; 
      fGammaW = 0.; 
      return 0; 
    }

  //set gamma 
  fGammaW = gamma; 

  // no error 
  return 1;
}

// --------------------------------------------------------- 
int KLFitter::PhysicsConstants::SetGammaTop(double gamma)
{
  // check argument 
  if (gamma < 0)
    {
      std::cout << "KLFitter::PhysicsConstants::SetGammaTop(). Width cannot be negative. Set width to zero." << std::endl; 
      fGammaTop = 0.; 
      return 0; 
    }

  //set gamma 
  fGammaTop = gamma; 

  // no error 
  return 1;
}

// --------------------------------------------------------- 
void KLFitter::PhysicsConstants::CalculateGammaTop()
{
  // use formula used in: 
  // T.~M.~Liss and A.~Quadt,
  // ``The Top Quark,'' in 
  // C.~Amsler {\it et al.}  [Particle Data Group],
  // ``Review of particle physics,''
  // Phys.\ Lett.\  B {\bf 667} (2008) 1.

  double gamma = fGF * fMassTop * fMassTop * fMassTop / 8 / fPi / sqrt(2); 
  gamma *= (1 - fMassW*fMassW/fMassTop/fMassTop) * (1 - fMassW*fMassW/fMassTop/fMassTop); 
  gamma *= (1 + 2 * fMassW * fMassW / fMassTop / fMassTop); 
  gamma *= (1 - 2 * fAlphaS / 3 / fPi * (2*fPi*fPi/3 - 5./2.));

  SetGammaTop(gamma);
}

// --------------------------------------------------------- 
