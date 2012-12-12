#include "PhysicsConstants.h" 
#include <cmath>
#include <iostream> 
#include <cmath>

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
  fMassTopUnc =   1.4; // top quark LHC uncertainty
  //++++++++++++++++//
  fMassHiggs    = 120.0; // Higgs mass in GeV/c^{2}
  fGammaHiggs   = 0.003512; // Higgs width

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
int KLFitter::PhysicsConstants::SetMassHiggs(double mass)
{
  // check argument 
  if (mass < 0)
    {
      std::cout << "KLFitter::PhysicsConstants::SetMassHiggs(). Mass cannot be negative. Set mass to zero." << std::endl; 
      fMassHiggs = 0.; 
      return 0; 
    }

  //set mass 
  fMassHiggs = mass; 

  // calculate top width
  CalculateGammaHiggs(); 

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
int KLFitter::PhysicsConstants::SetGammaHiggs(double gamma)
{
  // check argument 
  if (gamma < 0)
    {
      std::cout << "KLFitter::PhysicsConstants::SetGammaHiggs(). Width cannot be negative. Set width to zero." << std::endl; 
      fGammaHiggs = 0.; 
      return 0; 
    }

  //set gamma 
  fGammaHiggs = gamma; 

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

  double gamma = fGF * fMassTop * fMassTop * fMassTop / 8 / M_PI / sqrt(2); 
  gamma *= (1 - fMassW*fMassW/fMassTop/fMassTop) * (1 - fMassW*fMassW/fMassTop/fMassTop); 
  gamma *= (1 + 2 * fMassW * fMassW / fMassTop / fMassTop); 
  gamma *= (1 - 2 * fAlphaS / 3 / M_PI * (2*M_PI*M_PI/3 - 5./2.));

  SetGammaTop(gamma);
}

// --------------------------------------------------------- 

void KLFitter::PhysicsConstants::CalculateGammaHiggs()
{
  // numbers calculated by HDECAY
  // A.Djouadi, J.Kalinowski and M.Spira, 
  // Comp. Phys. Commun. 108 C (1998) 56, hep-ph/9704448. 
  // 
  //   MHSM(GeV)     WIDTH(GeV/c2)
  //   __________________________
  //
  //    110.000       0.2849E-02	
  //    115.000	      0.3124E-02
  //    120.000       0.3512E-02
  //    125.000       0.4078E-02
  //    130.000       0.4920E-02
  //    140.000	      0.8196E-02
  
  double gamma = 0.0;

  if ( fMassHiggs < 111.0)
    gamma = 0.002849;
  if ( fMassHiggs > 114.0 && fMassHiggs < 116.0)
    gamma = 0.003124;
  if ( fMassHiggs > 119.0 && fMassHiggs < 121.0)
    gamma = 0.003512;
  if ( fMassHiggs > 124.0 && fMassHiggs < 126.0)
    gamma = 0.004078;
  if ( fMassHiggs > 129.0 && fMassHiggs < 131.0)
    gamma = 0.004920;
  if ( fMassHiggs > 139.0  )
    gamma = 0.008196;

  SetGammaHiggs(gamma);
}








// --------------------------------------------------------- 
