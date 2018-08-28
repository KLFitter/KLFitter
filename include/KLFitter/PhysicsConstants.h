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

#ifndef KLFITTER_PHYSICSCONSTANTS_H_
#define KLFITTER_PHYSICSCONSTANTS_H_

// ---------------------------------------------------------

namespace KLFitter {
/**
 * \class KLFitter::PhysicsConstants
 * \brief A class containing physics constants.
 *
 * This class contains physics constants.
 */
class PhysicsConstants final {
 public:
  /**
   * The default constructor.
   */
  PhysicsConstants();

  /**
   * The (defaulted) destructor.
   */
  ~PhysicsConstants();

  /** \name Member functions (Get)  */
  /** @{ */

  /**
   * Return the mass of the bottom quark in GeV/c2.
   * @return The mass of the particle in GeV/c2.
   */
  double MassBottom() { return fMassBottom; }

  /**
   * Return the msas of the W boson in GeV/c2.
   * @return The mass of the particle in GeV/c2.
   */
  double MassW() { return fMassW; }

  /**
   * Return the mass of the Z boson in GeV/c2.
   * @return The mass of the particle in GeV/c2.
   */
  double MassZ() { return fMassZ; }

  /**
   * Return the msas of the top quark in GeV/c2
   * @return The mass of the particle in GeV/c2.
   */
  double MassTop() { return fMassTop; }

  /**
   * Return the msas of the Higgs boson in GeV/c2
   * @return The mass of the particle in GeV/c2.
   */
  double MassHiggs() { return fMassHiggs; }

  /**
   * Return the width of the W boson in GeV/c2.
   * @return The width of the particle in GeV/c2.
   */
  double GammaW() { return fGammaW; }

  /**
   * Return the width of the Z boson in GeV/c2.
   * @return The width of the particle in GeV/c2.
   */
  double GammaZ() { return fGammaZ; }

  /**
   * Return the width of the top quark in GeV/c2
   * @return The width of the particle in GeV/c2.
   */
  double GammaTop() { return fGammaTop; }

  /**
   * Return the width of the Higgs boson in GeV/c2
   * @return The width of the particle in GeV/c2.
   */
  double GammaHiggs() { return fGammaHiggs; }

  double MassTopUnc() { return fMassTopUnc; }

  /** @} */
  /** \name Member functions (Set)  */
  /** @{ */

  /**
   * Set the mass of the bottom quark in GeV/c2.
   * @param mass The mass of the particle in GeV/c2.
   * @return An error code.
   */
  int SetMassBottom(double mass);

  /**
   * Set the mass of the top quark in GeV/c2.
   * @param mass The mass of the particle in GeV/c2.
   * @return An error code.
   */
  int SetMassTop(double mass);

  /**
   * Set the mass of the Higgs boson in GeV/c2.
   * @param mass The mass of the particle in GeV/c2.
   * @return An error code.
   */
  int SetMassHiggs(double mass);

  /**
   * Set the mass of the W boson in GeV/c2.
   * @param mass The mass of the particle in GeV/c2.
   * @return An error code.
   */
  int SetMassW(double mass);

  /**
   * Set the mass of the Z boson in GeV/c2.
   * @param mass The mass of the particle in GeV/c2.
   * @return An error code.
   */
  int SetMassZ(double mass);

  /**
   * Set the width of the W boson in GeV/c2.
   * @param gamma The width of the particle in GeV/c2.
   * @return An error code.
   */
  int SetGammaW(double gamma);

  /**
   * Set the width of the Z boson in GeV/c2.
   * @param gamma The width of the particle in GeV/c2.
   * @return An error code.
   */
  int SetGammaZ(double gamma);

  /**
   * Set the width of the top quark in GeV/c2.
   * @param gamma The width of the particle in GeV/c2.
   * @return An error code.
   */
  int SetGammaTop(double gamma);

  /**
   * Set the width of the Higgs boson in GeV/c2.
   * @param gamma The width of the particle in GeV/c2.
   * @return An error code.
   */
  int SetGammaHiggs(double gamma);

  /** @} */
  /** \name Member functions (misc)  */
  /** @{ */

  /**
   * Calculates the top width at NLO.
   */
  void CalculateGammaTop();

  /**
   * Calculates the Higgs width using HDECAY.
   */
  void CalculateGammaHiggs();

  /** @} */

 private:
  /**
   * The bottom quark pole mass in GeV/c2.
   */
  double fMassBottom;

  /**
   * The W boson pole mass in GeV/c2.
   */
  double fMassW;

  /**
   * The Z boson pole mass in GeV/c2.
   */
  double fMassZ;

  /**
   * The top quark pole mass in GeV/c2.
   */
  double fMassTop;

  /**
   * The Higgs boson mass in GeV/c2.
   */
  double fMassHiggs;

  /**
   * The W boson width in GeV/c2.
   */
  double fGammaW;

  /**
   * The Z boson width in GeV/c2.
   */
  double fGammaZ;

  /**
   * The top quark width in GeV/c2.
   */
  double fGammaTop;

  /**
   * The Higgs boson width in GeV/c2.
   */
  double fGammaHiggs;

  /**
   * The Fermi constant.
   */
  double fGF;

  /**
   * alpha_S at m_Z
   */
  double fAlphaS;

  /**
   * The top quark mass LHC uncertainty in GeV/c2.
   */
  double fMassTopUnc;
};
}  // namespace KLFitter

#endif  // KLFITTER_PHYSICSCONSTANTS_H_
